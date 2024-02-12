from defusedxml import sax
import plyvel
import requests
import pandas as pd
import tables as tb
from rich.progress import track
from requests.adapters import HTTPAdapter, Retry

import re
import os
import csv
import gzip
from pathlib import Path
from random import sample
from itertools import product
from typing import List, Optional
from collections import defaultdict
from xml.sax.handler import ContentHandler, feature_namespaces

from ppi_origami.sources.utils import get_seq, download_file


def download(raw_folder):
    download_file(
        "https://omabrowser.org/All/oma-groups.orthoXML.xml.gz",
        raw_folder / "oma-groups.orthoXML.xml.gz",
        "Downloading Groups",
    )
    download_file(
        "https://omabrowser.org/All/oma-uniprot.txt.gz",
        raw_folder / "oma-uniprot.txt.gz",
        "Downloading OMA/UniProt Maps",
    )


class OMAIDHandler(ContentHandler):
    def __init__(self, limit_taxons: Optional[List[int]]):
        super().__init__()
        self.limit_taxons = limit_taxons
        self.current_species = None
        self.ids = []
        self.gene_ids = []
        self.prot_ids = []
        self.species = []
        self.current_ortholog_group = None
        self.ortholog_groups = defaultdict(lambda: set())

    def startElement(self, tag, attributes):
        if tag == "species":
            # print(attributes.getNames())
            self.current_species = int(attributes.getValue("NCBITaxId"))
        elif tag == "gene":
            if (
                    self.limit_taxons is None
                    or int(self.current_species) in self.limit_taxons
            ):
                self.ids.append(int(attributes.getValue("id")))
                self.gene_ids.append(attributes.getValue("geneId"))
                self.prot_ids.append(attributes.getValue("protId"))
                self.species.append(self.current_species)
        elif tag == "orthologGroup":
            self.current_ortholog_group = int(attributes.getValue("id"))
        elif tag == "geneRef":
            self.ortholog_groups[self.current_ortholog_group].add(
                int(attributes.getValue("id"))
            )

    def endElement(self, tag):
        pass

    def __del__(self):
        pass


def upkb_groups(raw_path, processed_path, limit_taxon):
    parser = sax.make_parser()
    parser.setFeature(feature_namespaces, 0)

    if os.path.isdir(processed_path / "omaid_to_uniprot.leveldb"):
        print("Loading existing omaid_to_uniprot db...")
        db = plyvel.DB(
            str(processed_path / "omaid_to_uniprot.leveldb"), create_if_missing=False
        )
    else:
        print("omaid_to_uniprot db not found, creating...")
        db = plyvel.DB(
            str(processed_path / "omaid_to_uniprot.leveldb"), create_if_missing=True
        )

        df = pd.read_csv(
            raw_path / "oma-uniprot.txt.gz", sep="\t", comment="#", header=None
        )

        for row in track(df.iterrows(), total=len(df)):
            db.put(
                str(row[1].values[0]).encode("utf8"),
                str(row[1].values[1]).encode("utf8"),
            )

    omaid_to_upkb = {}
    omaid_to_species = {}

    handler = OMAIDHandler(limit_taxon)
    parser.setContentHandler(handler)

    print("Parsing XML...")
    with gzip.open(raw_path / "oma-groups.orthoXML.xml.gz") as f:
        parser.parse(f)

    print("Prepping IDs...")
    for oma_id, prot_id, gene_species in track(
        zip(handler.ids, handler.prot_ids, handler.species), total=len(handler.ids)
    ):
        upkb = db.get(prot_id.encode("utf8"))

        omaid_to_upkb[int(oma_id)] = upkb.decode("utf8") if upkb is not None else None
        omaid_to_species[int(oma_id)] = int(gene_species)

    print("Saving OrthologGroups")
    path = processed_path

    os.makedirs(path, exist_ok=True)

    uids = set()

    upkb_ac_pattern = (
        r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
    )

    with gzip.open(path / "ortholog_group.csv.gz", "wt") as f:
        writer = csv.DictWriter(
            f, fieldnames=["orthologGroup_id", "members", "members_species"]
        )
        writer.writeheader()

        for og_id in track(handler.ortholog_groups):
            new_og = list()
            members_species = list()

            for oid in handler.ortholog_groups[og_id]:

                try:
                    uid = omaid_to_upkb[int(oid)]
                except KeyError:
                    print(f"Can't find UniprotKB accessions for the OMA id {oid}. Might be out-of-species. Skipping...")
                    continue

                if (
                    uid is not None
                    and re.match(upkb_ac_pattern, uid) is not None
                    and uid not in new_og
                ):
                    new_og.append(uid)
                    members_species.append(str(omaid_to_species[int(oid)]))
                    uids.add((uid, str(omaid_to_species[int(oid)])))

            if len(new_og) > 1 and len(members_species) > 0:
                writer.writerow(
                    {
                        "orthologGroup_id": og_id,
                        "members": "|".join(new_og),
                        "members_species": "|".join(members_species),
                    }
                )

    print("Saving UPKB IDs and sequences")

    with gzip.open(path / "upkb_ids.txt.gz", "wt") as f:
        for uid, species in uids:
            f.write(f"{uid},{species}\n")


def get_eval_proteins(path: Path) -> List[str]:
    """
    Extract testing and validation proteins from a RAPPPID dataset
    :param path: path to the RAPPPID dataset
    :return:
    """

    pids = set()

    with tb.open_file(path) as dataset:
        for pid in dataset.root.splits.test.iterrows():
            try:
                pids.add(pid[0].decode("utf8"))
            except AttributeError:
                continue

        for pid in dataset.root.splits.val.iterrows():
            try:
                pids.add(pid[0].decode("utf8"))
            except AttributeError:
                continue

    return list(pids)


def build_triplet_dataset(
    processed_path,
    uniref_threshold: int = 90,
    blacklist_species: Optional[List[int]] = None,
    whitelist_species: Optional[List[int]] = None,
    preloaded_split: Optional[Path] = None,
):

    uniref_db_path = str(processed_path / f"uniref{uniref_threshold}_members_upkb.leveldb")

    if not os.path.isdir(uniref_db_path):
        raise IOError(
            f"Can't read the UniRef{uniref_threshold}/UPKB mapping database at {uniref_db_path}. Make sure to run the 'process uniref' task."
        )

    uniref_db = plyvel.DB(uniref_db_path)

    if preloaded_split:
        # Load preloaded_split data
        eval_upkb_acs = get_eval_proteins(preloaded_split)
        eval_unirefs = [
            uniref_db.get(upkb_ac.encode("utf8")).decode("utf8")
            for upkb_ac in eval_upkb_acs
        ]
    else:
        eval_unirefs = []

    retries = Retry(
        total=50, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504]
    )
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    def include_species(taxon):
        if whitelist_species is not None:
            return int(taxon) in whitelist_species
        elif blacklist_species is not None:
            return int(taxon) not in blacklist_species
        else:
            return True

    path = processed_path
    df = pd.read_csv(path / "ortholog_group.csv.gz").sample(frac=1.0)

    db = plyvel.DB(
        str(processed_path / "uniprot_sequences.leveldb"), create_if_missing=True
    )

    upkbs = []

    with gzip.open(path / "upkb_ids.txt.gz", "rt") as f:
        for line in f:
            row = line.strip().split(",")
            upkb_ac = row[0]
            species = row[1]

            if include_species(species):
                upkbs.append(upkb_ac)

    groups = []

    print("Filtering groups...")
    for row_idx, row in track(df.iterrows(), total=len(df)):
        members = str(row.members).split("|")
        member_species = str(row.members_species).split("|")

        filtered_members = []

        for m, ms in zip(members, member_species):
            if ms == "nan":
                continue

            uniprot_id = uniref_db.get(upkb_ac.encode("utf8"))

            if (
                include_species(ms)
                and uniprot_id is not None
                and uniprot_id.decode("utf8") not in eval_unirefs
            ):
                filtered_members.append(m)

        if len(filtered_members) < 2:
            continue

        groups.append(filtered_members)

    print("Making triplets and retrieving sequences...")

    found_seq = set()

    i = 0

    with gzip.open(path / "triplets.txt.gz", "wt") as f:
        with gzip.open(path / "upkb_seqs.csv.gz", "wt") as f2:
            writer_seqs = csv.DictWriter(f2, fieldnames=["upkb_ac", "seq"])
            writer_seqs.writeheader()

            writer_triplets = csv.DictWriter(
                f, fieldnames=["anchor", "positive", "negative"]
            )
            writer_triplets.writeheader()

            for group in track(groups):
                pairs = [
                    _
                    for _ in product(group[: len(group) // 2], group[len(group) // 2 :])
                ]

                for pair in pairs:
                    seq_anchor = get_seq(pair[0], db, session)
                    seq_positive = get_seq(pair[1], db, session)

                    if (
                        "nan" in [seq_anchor, seq_positive]
                        or None in [seq_anchor, seq_positive]
                        or not isinstance(seq_anchor, str)
                        or not isinstance(seq_positive, str)
                        or "" in [seq_anchor, seq_positive]
                    ):
                        print("Skipping NaN")
                        continue

                    seq_negative = "nan"

                    while (
                        seq_negative == "nan"
                        or seq_negative is None
                        or not isinstance(seq_negative, str)
                        or seq_negative == ""
                    ):
                        print("Skipping NaN")
                        negative = sample(upkbs, 1)
                        seq_negative = get_seq(negative[0], db, session)

                    i += 1

                    if i % 100 == 0:
                        print(f"Wrote {i} rows...")

                    writer_triplets.writerow(
                        {
                            "anchor": pair[0],
                            "positive": pair[1],
                            "negative": negative[0],
                        }
                    )

                    if pair[0] not in found_seq:
                        writer_seqs.writerow({"upkb_ac": pair[0], "seq": seq_anchor})

                        found_seq.add(pair[0])

                    if pair[1] not in found_seq:
                        writer_seqs.writerow({"upkb_ac": pair[1], "seq": seq_positive})

                        found_seq.add(pair[1])

                    if negative[0] not in found_seq:
                        writer_seqs.writerow(
                            {"upkb_ac": negative[0], "seq": seq_negative}
                        )

                        found_seq.add(negative[0])
