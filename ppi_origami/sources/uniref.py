import csv
import gzip
from pathlib import Path
from xml.sax.handler import ContentHandler, feature_namespaces

import humanize
import plyvel
from defusedxml import sax

from . import utils


class UniRefHandler(ContentHandler):
    def __init__(self, save_folder: Path, threshold: int, console):
        super().__init__()
        self.console = console
        self.threshold = threshold
        self.in_entry = False
        self.in_representativeMember = False
        self.in_sequence = False
        self.in_member = False
        self.in_dbReference = False
        self.num_records = 0
        self.entry = {
            "entry_id": None,
            "updated": None,
            "sequence": "",
            "members_upkb": [],
            "members_uniparc": [],
        }

        self.db_members_upkb = plyvel.DB(
            str(save_folder / f"uniref{threshold}_members_upkb.leveldb"), create_if_missing=True
        )

        self.db_members_uniparc = plyvel.DB(
            str(save_folder / f"uniref{threshold}_members_uniparc.leveldb"), create_if_missing=True
        )

        self.db_sequences = plyvel.DB(
            str(save_folder / f"uniref{threshold}_sequences.leveldb"), create_if_missing=True
        )

        self.fh_members_upkb = gzip.open(
            save_folder / f"uniref{threshold}_members_upkb.csv.gz", "wt"
        )
        self.fh_members_uniparc = gzip.open(
            save_folder / f"uniref{threshold}_members_uniparc.csv.gz", "wt"
        )
        self.fh_sequences = gzip.open(save_folder / f"uniref{threshold}_sequences.csv.gz", "wt")

        fieldnames = ["member_upkb_ac", f"uniref{threshold}"]
        self.writer_members_upkb = csv.DictWriter(
            self.fh_members_upkb, fieldnames=fieldnames
        )
        self.writer_members_upkb.writeheader()

        fieldnames = ["member_uniparc", f"uniref{threshold}"]
        self.writer_members_uniparc = csv.DictWriter(
            self.fh_members_uniparc, fieldnames=fieldnames
        )
        self.writer_members_uniparc.writeheader()

        fieldnames = [f"uniref{threshold}", "sequence"]
        self.writer_sequences = csv.DictWriter(self.fh_sequences, fieldnames=fieldnames)
        self.writer_sequences.writeheader()

    def startElement(self, tag, attributes):
        if tag == "entry":
            self.in_entry = True
            self.entry["entry_id"] = attributes["id"]
            self.entry["updated"] = attributes["updated"]

        elif tag == "representativeMember" and self.in_entry:
            self.in_representativeMember = True

        elif tag == "sequence" and self.in_entry and self.in_representativeMember:
            self.in_sequence = True

        elif tag == "member":
            self.in_member = True

        elif tag == "dbReference" and (self.in_member or self.in_representativeMember):
            self.in_dbReference = True
            if attributes["type"] == "UniParc ID":
                self.entry["members_uniparc"].append(attributes["id"])

        elif (
                tag == "property"
                and (self.in_member or self.in_representativeMember)
                and self.in_dbReference
        ):
            if attributes["type"] == "UniProtKB accession":
                self.entry["members_upkb"].append(attributes["value"])

    def characters(self, content):
        if self.in_sequence:
            self.entry["sequence"] += content

    def endElement(self, tag):
        if tag == "entry":
            self.in_entry = False

            for member in self.entry["members_upkb"]:
                self.writer_members_upkb.writerow(
                    {"member_upkb_ac": member, f"uniref{self.threshold}": self.entry["entry_id"]}
                )

                self.db_members_upkb.put(
                    member.encode("utf8"), self.entry["entry_id"].encode("utf8")
                )

            for member in self.entry["members_uniparc"]:
                self.writer_members_uniparc.writerow(
                    {"member_uniparc": member, f"uniref{self.threshold}": self.entry["entry_id"]}
                )

                self.db_members_uniparc.put(
                    member.encode("utf8"), self.entry["entry_id"].encode("utf8")
                )

            self.writer_sequences.writerow(
                {f"uniref{self.threshold}": self.entry["entry_id"], "sequence": self.entry["sequence"]}
            )

            self.db_sequences.put(
                self.entry["entry_id"].encode("utf8"),
                self.entry["sequence"].encode("utf8"),
            )

            self.num_records += 1

            if self.num_records % 100000 == 0:
                self.console.log(
                    f"Processed {humanize.intword(self.num_records)} records."
                )

            self.entry = {
                "entry_id": None,
                "updated": None,
                "sequence": "",
                "members_upkb": [],
                "members_uniparc": [],
            }

        elif tag == "representativeMember":
            self.in_representativeMember = False

        elif tag == "sequence":
            self.in_sequence = False

        elif tag == "member":
            self.in_member = False

        elif tag == "dbReference":
            self.in_dbReference = False

    def __del__(self):
        self.fh_members_upkb.close()
        self.fh_members_uniparc.close()
        self.fh_sequences.close()


def process(raw_path: Path, processed_folder: Path, threshold: int, console):
    parser = sax.make_parser()
    parser.setFeature(feature_namespaces, 0)

    Handler = UniRefHandler(processed_folder, threshold, console)
    parser.setContentHandler(Handler)

    with gzip.open(raw_path / f"uniref_uniref{threshold}.xml.gz") as f:
        parser.parse(f)

    build_db(processed_folder, threshold)


def build_db(processed_folder: Path, threshold: int):
    db = plyvel.DB(
        str(processed_folder / f"uniref{threshold}_sequences.leveldb"), create_if_missing=True
    )

    with gzip.open(processed_folder / f"uniref{threshold}_sequences.csv.gz", "rt") as f:
        reader = csv.DictReader(f)

        for row_idx, row in enumerate(reader):
            name = row[f"uniref{threshold}"].encode("utf8")
            sequence = row["sequence"].encode("utf8")
            db.put(name, sequence)


def download(path: Path, threshold: int):
    if threshold not in [50, 90, 100]:
        raise ValueError(f"Invalid threshold value. Got {threshold}, expected one of 50, 90, or 100.")

    root = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/"
    xml_url = f"{root}/uniref{threshold}/uniref{threshold}.xml.gz"
    rel_url = f"{root}/uniref{threshold}/uniref{threshold}.release_note"

    utils.download_file(xml_url, path / f"uniref_uniref{threshold}.xml.gz", f"Downloading Uniref{threshold}")
    utils.download_file(rel_url, path / f"uniref_uniref{threshold}.release_note",
                        f"Downloading Uniref{threshold} Release Notes")
