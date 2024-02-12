import csv
import gzip
import json
import os
from pathlib import Path
from typing import Optional

import plyvel
import pandas as pd
from rich.progress import track

from . import utils


def load_upkb_aliases(raw_folder: Path, processed_folder: Path, version: str, console) -> dict:
    if os.path.isfile(processed_folder / "string_upkb_aliases.json"):
        console.log("Aliases file 'string_upkb_aliases.json' found. Loading them.")
        with open(processed_folder / "string_upkb_aliases.json") as f:
            aliases = json.load(f)
    else:
        console.log("Creating alias file 'string_upkb_aliases.json'")
        aliases = {}
        with gzip.open(
            raw_folder / f"string_protein.aliases.v{version}.txt.gz", "rt"
        ) as f_in:
            reader = csv.DictReader(f_in, delimiter="\t")
            for row in reader:
                if "UniProt_AC" in row["source"]:
                    aliases[row["#string_protein_id"]] = row["alias"]

        with open(processed_folder / "string_upkb_aliases.json", "w") as f_out:
            json.dump(aliases, f_out)

    return aliases


def update_secondaries(
    raw_folder: Path, processed_folder: Path, aliases: dict, console
) -> dict:
    if not os.path.isfile(raw_folder / "sec_ac.txt"):
        console.log(
            "WARNING: Did not find a secondary ac index",
            style="red",
        )
        console.log(
            'run "ppi_origami download uniprot_sec_ac ./raw" to download the secondary ac index.'
        )
        raise Exception("Did not find a secondary ac index. run \"ppi_orgiami download uniprot_sec_ac RAW_FOLDER\" to download the secondary ac index.")
    else:
        if os.path.isfile(processed_folder / "sec_ac.json"):
            console.log("Using cached sec_ac.json file")
            with open(processed_folder / "sec_ac.json") as f:
                sec_ac_map = json.load(f)
        else:
            console.log("Loading and transforming sec_ac.txt")
            sec_ac_map = pd.read_csv(
                raw_folder / "sec_ac.txt",
                skiprows=30,
                sep="\s+",
                names=["secondary", "primary"],
                header=0,
            )
            sec_ac_map = sec_ac_map.set_index("secondary").to_dict()
            console.log("Saving cached sec_ac.json")
            with open(processed_folder / "sec_ac.json", "w") as f:
                json.dump(sec_ac_map, f)

        for key in track(aliases, description="Updating aliases..."):
            try:
                aliases[key] = sec_ac_map["primary"][aliases[key]]
            except KeyError:
                continue

    return aliases


def process_string_upkb(
    raw_folder: Path, processed_path: Path, taxon: Optional[int], version: str, console
):
    links_path = raw_folder / f"string_{taxon}.protein.links.detailed.v{version}.txt.gz"

    if not os.path.exists(links_path):
        console.log(f"Didn't find STRING links TSV path, expected at {links_path}")
        return False

    console.log("Loading links as Pandas DataFrame")

    aliases = load_upkb_aliases(raw_folder, processed_path, version, console)

    aliases = update_secondaries(raw_folder, processed_path, aliases, console)

    console.log("Joining UPKB aliases...")
    missed_ids = 0
    total_ids = 0
    with gzip.open(
        raw_folder / f"string_{taxon}.protein.links.detailed.v{version}.txt.gz", "rt"
    ) as f_in:
        with gzip.open(
                processed_path
                / f"string_{taxon}.protein.links.detailed.v{version}_upkb.csv.gz",
            "wt",
        ) as f_out:
            reader = csv.DictReader(f_in, delimiter=" ")
            writer = csv.DictWriter(
                f_out,
                [
                    "protein1",
                    "protein2",
                    "protein1_upkb",
                    "protein2_upkb",
                    "experimental",
                    "database",
                    "textmining",
                    "combined_score",
                ],
            )
            writer.writeheader()

            for row_idx, row in enumerate(reader):
                try:
                    upkb1 = aliases[row["protein1"]]
                except KeyError:
                    console.log(
                        f"Couldn't find {row['protein1']} @ row {row_idx}", style="red"
                    )
                    missed_ids += 1
                    continue
                try:
                    upkb2 = aliases[row["protein2"]]
                except KeyError:
                    console.log(
                        f"Couldn't find {row['protein2']} @ row {row_idx}", style="red"
                    )
                    missed_ids += 1
                    continue

                writer.writerow(
                    {
                        "protein1": row["protein1"],
                        "protein2": row["protein2"],
                        "protein1_upkb": upkb1,
                        "protein2_upkb": upkb2,
                        "experimental": row["experimental"],
                        "textmining": row["textmining"],
                        "combined_score": row["combined_score"],
                    }
                )

                total_ids += 1

    if missed_ids > 0:
        console.log(
            f"WARNING: {missed_ids}/{total_ids} ids don't have corresponding UPKB ids"
        )


def process_string_uniref(
    processed_path: Path, threshold: int, taxon: Optional[int], version: str, console
):

    uniref_db_path = str(processed_path / f"uniref{threshold}_members_upkb.leveldb")

    if not os.path.isdir(uniref_db_path):
        raise IOError(
            f"Can't read the UniRef90/UPKB mapping database at {uniref_db_path}. Make sure to run the 'process uniref' task."
        )

    uniref_db = plyvel.DB(uniref_db_path)

    console.log(f"Looking up UniRef{threshold} aliases in STRING UPKB")
    links_path = (
            processed_path / f"string_{taxon}.protein.links.detailed.v{version}_upkb.csv.gz"
    )
    missed_ids = 0
    total_ids = 0
    with gzip.open(
        processed_path / f"string_{taxon}.protein.links.detailed.v{version}_uniref{threshold}.csv.gz",
        "wt",
    ) as f_out:
        writer = csv.DictWriter(
            f_out,
            [
                "protein1",
                "protein2",
                "protein1_upkb",
                "protein2_upkb",
                f"protein1_uniref{threshold}",
                f"protein2_uniref{threshold}",
                "experimental",
                "database",
                "textmining",
                "combined_score",
            ],
        )
        writer.writeheader()
        with gzip.open(links_path, "rt") as f_in:
            reader = csv.DictReader(f_in)
            # missed = False
            for row_idx, row in enumerate(reader):
                if (
                    row_idx == 0
                    and ("protein1_upkb" not in row.keys()
                    or "protein2_upkb" not in row.keys())
                ):
                    raise KeyError(
                        "Please run string_upkb before running string_uniref"
                    )

                uniref_1 = uniref_db.get(str(row["protein1_upkb"]).encode("utf8"))
                uniref_2 = uniref_db.get(str(row["protein2_upkb"]).encode("utf8"))

                if uniref_1 is None or uniref_2 is None:
                    missed_ids += 1

                if uniref_1 is None:
                    uniref_1 = "n.a"
                else:
                    uniref_1 = uniref_1.decode("utf8")

                if uniref_2 is None:
                    uniref_2 = "n.a"
                else:
                    uniref_2 = uniref_2.decode("utf8")

                writer.writerow(
                    {
                        "protein1": row["protein1"],
                        "protein2": row["protein2"],
                        "protein1_upkb": row["protein1_upkb"],
                        "protein2_upkb": row["protein2_upkb"],
                        f"protein1_uniref{threshold}": uniref_1,
                        f"protein2_uniref{threshold}": uniref_2,
                        "experimental": row["experimental"],
                        "database": row["database"],
                        "textmining": row["textmining"],
                        "combined_score": row["combined_score"],
                    }
                )

                total_ids += 1
                if row_idx % 100000 == 0:
                    console.log("Wrote 100K lines...")

    console.log(f"{missed_ids / total_ids * 100:.4}% of IDs are missing uniref{threshold} ids")


def parse_string(source_path: Path, identifier: str):
    """
    Parses STRING rows and output the fields required to construct a commmon file.
    """

    first_row = True

    with gzip.open(source_path, "rt") as f_in:
        reader = csv.DictReader(f_in)

        for row in reader:
            if first_row:
                first_row = False
                keys = list(row.keys())
                if (
                    f"protein1_{identifier}" not in keys
                    and f"protein2_{identifier}" not in keys
                ):
                    raise ValueError(
                        f"Expected to find headers 'protein1_{identifier}' and 'protein2_{identifier}' \
                                            'in {source_path}, but did not find it."
                    )

            score = utils.serialize_dict(
                {
                    "string_combined_score": row["combined_score"],
                    "string_experimental": row["experimental"],
                    "string_database": row["database"],
                    "string_textmining": row["textmining"],
                }
            )

            id_a = row[f"protein1_{identifier}"]
            id_b = row[f"protein2_{identifier}"]
            interaction_id = utils.canonical_interaction_id(id_a, id_b)

            yield {
                "interaction_id": interaction_id,
                "protein1": id_a,
                "protein2": id_b,
                "score": score,
            }


def download(dataset_type: str, taxon: Optional[int], path: Path, version: str):

    root_link = "https://stringdb-downloads.org/download/"

    if dataset_type == "links" and taxon is None:
        url = f"{root_link}protein.links.v{version}/protein.links.v{version}.txt.gz"
        utils.download_file(
            url,
            path / f"string_protein.links.v{version}.txt.gz",
            "STRING Protein Links (Taxon: all)",
        )
    elif dataset_type == "links":

        if taxon is None:
            raise ValueError("Taxon must not be None.")

        url = f"{root_link}/protein.links.v{version}/{taxon}.protein.links.v{version}.txt.gz"
        utils.download_file(
            url,
            path / f"string_{taxon}.string_protein.links.v{version}.txt.gz",
            f"STRING Protein Links (Taxon: {taxon})",
        )
    elif dataset_type == "links_detailed" and taxon is None:
        url = f"{root_link}/protein.links.detailed.v{version}.txt.gz"
        utils.download_file(
            url,
            path / f"string_protein.links.detailed.v{version}.txt.gz",
            "STRING Protein Links Detailed (Taxon: all)",
        )
    elif dataset_type == "links_detailed":

        if taxon is None:
            raise ValueError("Taxon must not be None.")

        url = f"{root_link}protein.links.detailed.v{version}/{taxon}.protein.links.detailed.v{version}.txt.gz"
        utils.download_file(
            url,
            path / f"string_{taxon}.protein.links.detailed.v{version}.txt.gz",
            f"STRING Protein Links Detailed (Taxon: {taxon})",
        )
    elif dataset_type == "physical_links" and taxon is None:

        if taxon is None:
            raise ValueError("Taxon must not be None.")

        url = f"{root_link}protein.physical.links.v{version}/protein.physical.links.v{version}.txt.gz"
        utils.download_file(
            url,
            path / f"string_protein.links.v{version}.txt.gz",
            f"STRING Protein Physical Links (Taxon: all)",
        )
    elif dataset_type == "physical_links":

        if taxon is None:
            raise ValueError("Taxon must not be None.")

        url = f"{root_link}protein.physical.links.v{version}/{taxon}.protein.physical.links.v{version}.txt.gz"
        utils.download_file(
            url,
            path / f"string_{taxon}.protein.links.v{version}.txt.gz",
            f"STRING Protein Physical Links (Taxon: {taxon})",
        )
    elif dataset_type == "physical_detailed_links" and taxon is None:

        if taxon is None:
            raise ValueError("Taxon must not be None.")

        url = f"{root_link}protein.physical.links.detailed.v{version}/protein.physical.links.detailed.v{version}.txt.gz"
        utils.download_file(
            url,
            path / f"string_protein.links.detailed.v{version}.txt.gz",
            f"STRING Protein Physical Links Detailed (Taxon: all)",
        )
    elif dataset_type == "physical_detailed_links":

        if taxon is None:
            raise ValueError("Taxon must not be None.")

        url = f"{root_link}protein.physical.links.detailed.v{version}/{taxon}.protein.physical.links.detailed.v{version}.txt.gz"
        utils.download_file(
            url,
            path / f"string_{taxon}.protein.links.detailed.v{version}.txt.gz",
            f"STRING Protein Physical Links Detailed (Taxon: {taxon})",
        )
    elif dataset_type == "aliases" and taxon is None:
        url = f"{root_link}/protein.aliases.v{version}.txt.gz"
        utils.download_file(
            url,
            path / f"string_protein.aliases.v{version}.txt.gz",
            "STRING Protein Aliases (Taxon: all)",
        )
    elif dataset_type == "aliases":
        url = f"{root_link}/protein.aliases.v{version}/{taxon}.protein.aliases.v{version}.txt.gz"
        utils.download_file(
            url,
            path / f"string_{taxon}.protein.aliases.v{version}.txt.gz",
            f"STRING Protein Aliases (Taxon: {taxon})",
        )
