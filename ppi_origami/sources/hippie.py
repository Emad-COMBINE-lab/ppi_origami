import gzip
import re
from csv import DictReader
from pathlib import Path

from . import utils


def download(path: Path, version: str = "current"):
    url = f"http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/HIPPIE-{current}.mitab.txt"
    utils.download_file(url, path / "HIPPIE-current.mitab.txt", "Downloading HIPPIE")


def join_upkb(raw_folder: Path, processed_folder: Path, console):
    with open(raw_folder / "HIPPIE-current.mitab.txt", "r") as f:
        reader = DictReader(f, delimiter="\t")

        id_map = {}

        for idx, row in enumerate(reader):
            ids_a = utils.mitab_deserialize(
                row["ID Interactor A"]
            ) | utils.mitab_deserialize(row["Alt IDs Interactor A"])
            ids_b = utils.mitab_deserialize(
                row["ID Interactor B"]
            ) | utils.mitab_deserialize(row["Alt IDs Interactor B"])

            if idx == 0:
                # https://www.uniprot.org/help/accession_numbers
                upkb_ac_re = re.compile(
                    "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
                )
                # https://www.uniprot.org/help/entry_name
                upkb_name_re = re.compile("[A-Za-z\d]{1,5}_[A-Za-z\d]{1,5}")

                if "uniprotkb" in ids_a and "uniprotkb" in ids_b:
                    id_key = "uniprotkb"
                    a_is_ac = upkb_ac_re.match(ids_a["uniprotkb"])
                    b_is_ac = upkb_ac_re.match(ids_b["uniprotkb"])
                    if a_is_ac and b_is_ac:
                        console.log(
                            "Found UPKB accession ids in ambiguous ID tag 'uniprotkb'"
                        )
                        id_is_ac = True
                    else:
                        id_is_ac = False

                    if id_is_ac is False:
                        a_is_upkb_name = upkb_name_re.match(ids_a["uniprotkb"])
                        b_is_upkb_name = upkb_name_re.match(ids_b["uniprotkb"])
                        if a_is_upkb_name and b_is_upkb_name:
                            console.log(
                                "Found UPKB entry name in ambiguous ID tag 'uniprotkb'"
                            )
                            id_is_upkb_name = True
                        else:
                            id_is_upkb_name = False

                    if id_is_upkb_name is False and id_is_ac is False:
                        console.log(
                            "Ambiguous tag 'uniprotkb' appears to be neither entry name or accession. Breaking."
                        )
                        return None
                else:
                    console.log(
                        "Did not find a tag we understand to be UniProt related. Breaking"
                    )
                    return None

            if id_is_ac:
                ac_id_a = ids_a[id_key]
                ac_id_b = ids_b[id_key]
            elif id_is_upkb_name:
                if len(id_map) == 0:
                    console.log("Loading ID UniProtKB name to ac map into memory...")
                    with gzip.open(
                            processed_folder / "uniprot_idmapping_UniProtKB-ID.csv.gz",
                        "rt",
                    ) as f:
                        map_reader = DictReader(f)

                        for row in map_reader:
                            id_map[row["UniProtKB-ID"]] = row["upkb_ac"]

                ac_id_a = id_map[ids_a[id_key]]
                ac_id_b = id_map[ids_b[id_key]]

                print(ac_id_a, ac_id_b)

            if idx < 3:
                break
