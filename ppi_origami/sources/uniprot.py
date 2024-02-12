import gzip
from csv import reader, DictWriter, DictReader
from pathlib import Path
from typing import Optional

import plyvel
import requests
from rich.progress import Progress

from . import utils


def build_ensembl_pro_to_upkb_db(processed_folder: Path):
    db = plyvel.DB(
        str(processed_folder / "ensembl_pro_upkb.leveldb"), create_if_missing=True
    )

    with gzip.open(
        processed_folder / "uniprot_idmapping_Ensembl_PRO.csv.gz", "rt"
    ) as f:
        r = DictReader(f)

        for row_idx, row in enumerate(r):
            # remove version numbers
            ensembl = row["Ensembl_PRO"].split(".")[0].encode("utf8")
            upkb_ac = row["upkb_ac"].split("-")[0].encode("utf8")
            db.put(ensembl, upkb_ac)


def make_seq_db(processed_folder, taxon: Optional[int] = None):

    if taxon is None:
        download_path = "https://rest.uniprot.org/uniprotkb/search?fields=accession%2Creviewed%2Cid%2Csequence%2Corganism_id&format=tsv&query=%28%2A%29&size=500"
    else:
        download_path = f"https://rest.uniprot.org/uniprotkb/search?fields=accession%2Creviewed%2Cid%2Csequence%2Corganism_id&format=tsv&query=%28%2A%29%20AND%20%28taxonomy_id%3A{taxon}%29&size=500"

    db_all = plyvel.DB(
        str(processed_folder / "uniprot_sequences.leveldb"), create_if_missing=True
    )

    if taxon is not None:
        db_taxon = plyvel.DB(
            str(processed_folder / f"uniprot_sequences_{taxon}.leveldb"),
            create_if_missing=True,
        )

    num_rows = 0

    with Progress() as progress:
        r = requests.request("HEAD", download_path)
        total = float(r.headers["x-total-results"])

        task = progress.add_task("[red]Downloading...", total=total)

        for row in utils.upkb_query(download_path):
            entry, reviewed, entry_name, sequence, organism = row

            db_all.put(entry.encode("utf8"), sequence.encode("utf8"))
            if taxon is not None:
                db_taxon.put(entry.encode("utf8"), sequence.encode("utf8"))

            num_rows += 1
            progress.update(task, advance=1)


def download(dataset_type: str, taxon: Optional[int], path: Path):
    if dataset_type == "sec_ac" and taxon is None:
        url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/sec_ac.txt"
        utils.download_file(
            url, path / "sec_ac.txt", "Downloading Uniprot Secondary Accession Numbers"
        )
    elif dataset_type == "id_mapping" and taxon is None:
        url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz"
        utils.download_file(
            url, path / "uniprot_idmapping.dat.gz", "Downloading Uniprot ID mapping"
        )
    elif dataset_type == "delac" and taxon is None:
        url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/delac_sp.txt"
        utils.download_file(
            url, path / "delac_sp.txt", "Downloading Deleted Swiss-Prot IDs"
        )
        url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/delac_tr.txt.gz"
        utils.download_file(
            url, path / "delac_tr.txt.gz", "Downloading Deleted TREmbl IDs"
        )

    elif dataset_type == "uniprot_seq_db":
        make_seq_db(path, taxon)
    elif dataset_type == "uniprot_seq_fasta":
        sprot_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
        trembl_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"

        utils.download_file(
            sprot_url, path / "uniprot_sprot.fasta.gz", "Downloading Swiss-Prot Sequences"
        )

        utils.download_file(
            trembl_url, path / "uniprot_trembl.fasta.gz", "Downloading TrEMBL Sequences"
        )
    else:
        raise ValueError


def process(dataset_type: str, raw_folder: Path, processed_folder: Path):
    if dataset_type == "id_mapping":
        seen_id_types = set()

        with gzip.open(raw_folder / "uniprot_idmapping.dat.gz", "rt") as f:
            csv_reader = reader(f, delimiter="\t")

            for upkb_ac, id_type, sec_id in csv_reader:
                if id_type in seen_id_types:
                    write_mode = "at"
                    write_header = False
                else:
                    write_mode = "wt"
                    write_header = True

                with gzip.open(
                    processed_folder / f"uniprot_idmapping_{id_type}.csv.gz", write_mode
                ) as f_out:
                    csv_writer = DictWriter(f_out, fieldnames=["upkb_ac", id_type])

                    if write_header:
                        csv_writer.writeheader()

                    csv_writer.writerow({"upkb_ac": upkb_ac, id_type: sec_id})

                seen_id_types.add(id_type)

