import os
import gzip
import random
from csv import DictReader
from pathlib import Path
from typing import List, Dict, Set, Optional
from rich.progress import track
import tables as tb
import plyvel
from ppi_origami.sources.utils import upkb_query
import requests
from requests.adapters import HTTPAdapter, Retry


class Protein(tb.IsDescription):
    name = tb.StringCol(20, pos=1)  # UniProt/Ref AC/ID
    sequence = tb.StringCol(3000, pos=2)  # Amino acid sequence


class Pair(tb.IsDescription):
    protein_id1 = tb.StringCol(20, pos=1)
    protein_id2 = tb.StringCol(20, pos=2)
    omid_protein_id = tb.StringCol(20, pos=3)
    omid_id = tb.Int64Col(pos=4)
    label = tb.BoolCol(pos=5)


class OMA(tb.IsDescription):
    ortholog_group_id = tb.Int64Col(pos=1)
    protein_id = tb.StringCol(20, pos=2)


def write_intrepppid(
    processed_folder: Path,
    path: Path,
    c_types: List[int],
    interactions,
    proteins: Dict[str, List[str]],
    scramble_interactions: bool,
    oma: Dict[int, List[str]],
    oma_ids: Set[int],
):
    db = plyvel.DB(
        str(processed_folder / "uniprot_sequences.leveldb"), create_if_missing=False
    )

    filters = tb.Filters(
        complevel=9, complib="blosc:zstd", fletcher32=True, bitshuffle=True
    )
    intrepppid_dataset = tb.open_file(str(path), mode="w", filters=filters)

    group_interactions = intrepppid_dataset.create_group(
        "/", "interactions", "Protein Interactions"
    )

    group = {
        "C3": intrepppid_dataset.create_group(group_interactions, "c3", "C3 Datasets"),
        "C2": intrepppid_dataset.create_group(group_interactions, "c2", "C2 Datasets"),
        "C1": intrepppid_dataset.create_group(group_interactions, "c1", "C1 Datasets"),
    }

    table_pairs = {f"C{c}": {} for c in c_types}

    for c in c_types:
        for split in ["train", "test", "val"]:
            table_pairs[f"C{c}"][split] = intrepppid_dataset.create_table(
                group[f"C{c}"],
                f"c{c}_{split}",
                Pair,
                f"C{c} {split.title()} Dataset",
                filters=filters,
            )

    group_splits = intrepppid_dataset.create_group("/", "splits", "Protein Splits")

    train_proteins = intrepppid_dataset.create_carray(
        group_splits,
        "train",
        tb.StringCol(20),
        (len(proteins["train"]),),
        filters=filters,
    )
    val_proteins = intrepppid_dataset.create_carray(
        group_splits, "val", tb.StringCol(20), (len(proteins["val"]),), filters=filters
    )
    test_proteins = intrepppid_dataset.create_carray(
        group_splits,
        "test",
        tb.StringCol(20),
        (len(proteins["test"]),),
        filters=filters,
    )

    table_sequences = intrepppid_dataset.create_table(
        "/", "sequences", Protein, "Protein Sequences", filters=filters
    )

    table_orthologs = intrepppid_dataset.create_table(
        "/", "orthologs", OMA, "Orthologs", filters=filters
    )

    for c_type in c_types:
        for split in ["train", "val", "test"]:
            for interaction in track(
                interactions[split][f"C{c_type}"], "Saving pairs..."
            ):
                pair = table_pairs[f"C{c_type}"][split].row
                pair["protein_id1"] = interaction["pid1"]
                pair["protein_id2"] = interaction["pid2"]
                pair["omid_id"] = interaction["oma_id"]
                pair["omid_protein_id"] = interaction["oma_pid"]
                pair["label"] = interaction["label"]
                pair.append()

    train_proteins[:] = list(proteins["train"])
    val_proteins[:] = list(proteins["val"])
    test_proteins[:] = list(proteins["test"])

    all_proteins = list(set(proteins["train"] + proteins["val"] + proteins["test"]))

    if scramble_interactions:
        scrambled_proteins = all_proteins.copy()
        random.shuffle(scrambled_proteins)

        protein_scramble_map = {
            protein: scrambled_protein
            for protein, scrambled_protein in zip(all_proteins, scrambled_proteins)
        }

    for oma_id in track(oma_ids, description="Writing OMA IDs"):

        for prot in oma[oma_id]:
            oma_row = table_orthologs.row
            oma_row["ortholog_group_id"] = oma_id
            oma_row["protein_id"] = prot
            all_proteins.append(prot)
            oma_row.append()

    missing_proteins = []

    for target_protein in track(all_proteins):
        sequence = table_sequences.row

        sequence["name"] = target_protein

        if scramble_interactions:
            target_protein = protein_scramble_map[target_protein]

        retries = Retry(
            total=50, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504]
        )
        session = requests.Session()
        session.mount("https://", HTTPAdapter(max_retries=retries))

        seq = db.get(target_protein.encode("utf8"))

        if seq is None:
            missing_proteins.append(target_protein.encode("utf8"))
        else:
            sequence["sequence"] = seq
            sequence.append()

    print(f"{len(missing_proteins)} missing proteins")
    missing_proteins = [f"accession:{missing_protein.decode('utf8')}" for missing_protein in missing_proteins]

    num_batches = len(missing_proteins) // 1000
    num_batches = max(num_batches, 1)

    for batch_num in track(range(num_batches), description="Finding Missing Sequences"):
        min_index = batch_num * 1000
        max_index = min_index + 1000
        batch = missing_proteins[min_index:max_index]
        query = "+OR+".join(batch)
        query_url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&format=tsv&fields=accession,sequence"

        for upkb_ac, seq in upkb_query(query_url):
            sequence = table_sequences.row
            sequence["name"] = upkb_ac
            sequence["sequence"] = seq
            sequence.append()
            db.put(upkb_ac.encode('utf8'), seq.encode('utf8'))

        intrepppid_dataset.flush()

    return intrepppid_dataset


def filter_members(
    members: List[str], member_species: List[int], allowlist_taxon: Optional[List[int]],
    denylist_taxon: Optional[List[int]]
):
    new_members = []

    # Whitelist is ignored if both blacklist and whitelist are specified.
    if allowlist_taxon and denylist_taxon:
        denylist_taxon = None

    if allowlist_taxon:
        for member, member_s in zip(members, member_species):
            if member_s in allowlist_taxon:
                new_members.append(member)

    if denylist_taxon:
        for member, member_s in zip(members, member_species):
            if member_s not in denylist_taxon:
                new_members.append(member)

    return new_members


def get_interactions(table, upkb_to_oma_id, oma, uniref_db, eval_unirefs):
    interactions = []
    oma_ids = set()

    for row_idx, row in enumerate(table):

        if row_idx % 1000 == 0:
            print(f"{row_idx}/{len(table)} ({row_idx*100/len(table):.4}%) rows.")

        pid1 = row["protein_id1"].decode("utf8")
        pid2 = row["protein_id2"].decode("utf8")
        label = row["label"]

        if pid1 in upkb_to_oma_id:
            oma_pid = pid1
            oma_id = upkb_to_oma_id[pid1]
        elif pid2 in upkb_to_oma_id:
            oma_pid = pid2
            oma_id = upkb_to_oma_id[pid2]
        else:
            while True:
                oma_id = random.sample(oma.keys(), 1)[0]
                members = oma[oma_id]
                oma_pid = random.sample(members, 1)[0]
                oma_uniref_pid = uniref_db.get(oma_pid.encode("utf8"))

                if oma_uniref_pid is not None and oma_uniref_pid not in eval_unirefs:
                    break

        interactions.append(
            {
                "pid1": pid1,
                "pid2": pid2,
                "label": label,
                "oma_pid": oma_pid,
                "oma_id": oma_id,
            }
        )

        oma_ids.add(oma_id)

    return interactions, oma_ids


def load_delac(raw_path: Path):
    delac_sp_path = str(raw_path / "delac_sp.txt")
    if not os.path.isfile(delac_sp_path):
        raise IOError(
            f"Can't read the UPKB delac database at {delac_sp_path}. Make sure to run the 'download uniprot_delac' task."
        )

    delac_tr_path = str(raw_path / "delac_tr.txt.gz")
    if not os.path.isfile(delac_sp_path):
        raise IOError(
            f"Can't read the UPKB delac database at {delac_tr_path}. Make sure to run the 'download uniprot_delac' task."
        )

    delac = set()

    for f in [open(delac_sp_path), gzip.open(delac_tr_path)]:
        start = False
        for line in f:
            if line.strip() == "________________":
                start = True
                continue

            if start is False:
                continue

            delac.add(line.strip())

        f.close()

    return delac


def rapppid_to_intrepppid(
    processed_path: Path,
    rapppid_path: Path,
    intrepppid_path: Path,
    c_types: List[int],
    allowlist_taxon: Optional[List[int]],
    denylist_taxon: Optional[List[int]],
    scramble_interactions: bool,
    scramble_orthologs: bool,
    uniref_threshold: int = 90,
    max_members: Optional[int] = 10
):

    if allowlist_taxon and denylist_taxon:
        print("Allowlist ignored since Denylist is specified.")

    rapppid_dataset = tb.open_file(str(rapppid_path))
    uniref_db_path = str(processed_path / f"uniref{uniref_threshold}_members_upkb.leveldb")
    if not os.path.isdir(uniref_db_path):
        raise IOError(
            f"Can't read the UniRef{uniref_threshold}/UPKB mapping database at {uniref_db_path}. Make sure to run the 'process uniref' task."
        )

    uniref_db = plyvel.DB(uniref_db_path)

    proteins_upkb = {"train": list(), "val": list(), "test": list()}

    for split_name, split in zip(
        ["train", "val", "test"],
        [
            rapppid_dataset.root.splits.train,
            rapppid_dataset.root.splits.val,
            rapppid_dataset.root.splits.test,
        ],
    ):
        for row in split.iterrows():
            name = row

            if not isinstance(name, int):
                name = name.decode("utf8")

            proteins_upkb[split_name].append(name)

    proteins_uniref = dict()

    for split in ["train", "val", "test"]:
        proteins_uniref[split] = set(
            [
                uniref_db.get(upkb_ac.encode("utf8")).decode("utf8")
                for upkb_ac in proteins_upkb[split]
                if not isinstance(upkb_ac, int)
            ]
        )

    oma = dict()
    upkb_to_oma_id = dict()

    num_oma = 0

    if scramble_orthologs:

        print("!!!! WARNING SCRAMBLING ORTHOLOG IDS")

        all_oma_ids = set()

        with gzip.open(
                processed_path / "ortholog_group.csv.gz", "rt", newline=""
        ) as f:
            oma_reader = DictReader(f)

            for row in oma_reader:
                og_id = int(row["orthologGroup_id"])
                all_oma_ids.add(og_id)

    if scramble_orthologs:
        all_oma_members = set()

    with gzip.open(
        processed_path / "ortholog_group.csv.gz", "rt", newline=""
    ) as f:
        oma_reader = DictReader(f)

        for row in oma_reader:
            # print(row)
            og_id = int(row["orthologGroup_id"])
            members = row["members"].split("|")
            members_species = [
                int(x) for x in row["members_species"].split("|") if x != ""
            ]

            if allowlist_taxon or denylist_taxon:
                members = filter_members(members, members_species, allowlist_taxon, denylist_taxon)

            if max_members:
                members = members[:max_members]

            if len(members) > 2:
                oma[og_id] = members

                for member in members:
                    if scramble_orthologs:
                        all_oma_members.add(member)
                    upkb_to_oma_id[member] = og_id

            num_oma += 1

        if scramble_orthologs:
            # reset the oma dictionary
            oma = dict()

            # convert to list so that we can sample
            all_oma_members = list(all_oma_members)
            all_oma_ids = list(all_oma_ids)

            for og_id in all_oma_ids:
                # assign random members to each ortholog group
                oma[og_id] = random.sample(all_oma_members, 10)

            for member in all_oma_members:
                # assign a random oma id to each member
                upkb_to_oma_id[member] = random.sample(all_oma_ids, 1)[0]


    print(f"Num OMA: {num_oma}")

    interactions = dict()

    for split in ["train", "val", "test"]:
        interactions[split] = dict()
        for c in ["C1", "C2", "C3"]:
            interactions[split][c] = []

    eval_proteins = proteins_uniref["val"].union(proteins_uniref["test"])

    oma_ids = set()

    for split, c, table in track([
        ("train", "C1", rapppid_dataset.root.interactions.c1.c1_train),
        ("train", "C2", rapppid_dataset.root.interactions.c2.c2_train),
        ("train", "C3", rapppid_dataset.root.interactions.c3.c3_train),
        ("val", "C1", rapppid_dataset.root.interactions.c1.c1_val),
        ("val", "C2", rapppid_dataset.root.interactions.c2.c2_val),
        ("val", "C3", rapppid_dataset.root.interactions.c3.c3_val),
        ("test", "C1", rapppid_dataset.root.interactions.c1.c1_test),
        ("test", "C2", rapppid_dataset.root.interactions.c2.c2_test),
        ("test", "C3", rapppid_dataset.root.interactions.c3.c3_test),
    ]):
        interactions[split][c], split_oma_ids = get_interactions(
            table, upkb_to_oma_id, oma, uniref_db, eval_proteins
        )

        oma_ids = oma_ids.union(split_oma_ids)

    print("Writing INTREPPPID...")
    h5file = write_intrepppid(
        processed_path,
        intrepppid_path,
        c_types,
        interactions,
        proteins_upkb,
        scramble_interactions,
        oma,
        oma_ids,
    )

    class Metadata(tb.IsDescription):
        rapppid_path = tb.StringCol(560, pos=1)
        c1 = tb.BoolCol(pos=2)
        c2 = tb.BoolCol(pos=3)
        c3 = tb.BoolCol(pos=4)
        scramble_interactions = tb.BoolCol(pos=5)
        scramble_orthologs = tb.BoolCol(pos=6)
        uniref_threshold = tb.IntCol(pos=7)
        max_members = tb.IntCol(pos=8)
        allowlist_taxon = tb.StringCol(560, pos=9)
        denylist_taxon = tb.StringCol(560, pos=10)

    table_metadata = h5file.create_table("/", "metadata", Metadata, "Metadata")

    if allowlist_taxon:
        allowlist_taxon = [str(taxon) for taxon in allowlist_taxon]
    if denylist_taxon:
        denylist_taxon = [str(taxon) for taxon in denylist_taxon]

    metadata = table_metadata.row
    metadata["rapppid_path"] = rapppid_path.stem
    metadata["c1"] = 1 in c_types
    metadata["c2"] = 2 in c_types
    metadata["c3"] = 3 in c_types
    metadata["scramble_interactions"] = scramble_interactions
    metadata["scramble_orthologs"] = scramble_orthologs
    metadata["uniref_threshold"] = uniref_threshold
    metadata["max_members"] = max_members
    metadata["allowlist_taxon"] = ", ".join(allowlist_taxon) if allowlist_taxon else ""
    metadata["denylist_taxon"] = ", ".join(denylist_taxon) if denylist_taxon else ""
    metadata.append()
