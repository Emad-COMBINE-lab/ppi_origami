import csv
import gzip
import json
import os
from csv import DictWriter, DictReader
from csv import reader as csv_reader
from pathlib import Path
from typing import Iterable
from typing import Optional, List

import plyvel

from . import utils, uniprot
from .string import load_upkb_aliases, update_secondaries


def parse_dscript(
    source_paths: Iterable[Path],
    identifier: Optional[str],
    interaction: Optional[int] = 1,
):
    if interaction not in [0, 1, None]:
        raise ValueError("Invalid value for `interaction`. Must be 0, 1, or None.")

    for source_path in source_paths:
        first_row = True

        with gzip.open(source_path, "rt") as f_in:
            reader = DictReader(f_in)

            for row in reader:
                if first_row:
                    first_row = False
                    keys = list(row.keys())
                    if (
                        f"protein1_{identifier}"
                        if identifier
                        else "protein1" not in keys and f"protein2_{identifier}"
                        if identifier
                        else "protein2" not in keys
                    ):
                        raise ValueError(
                            f"Expected to find headers 'protein1_{identifier}' and 'protein2_{identifier}' \
                                                'in {source_path}, but did not find it."
                            if identifier
                            else f"Expected to find headers 'protein1' and 'protein2' \
                                                'in {source_path}, but did not find it."
                        )

                if interaction is None or int(row["interaction"]) == interaction:
                    id_a = row[f"protein1_{identifier}"]
                    id_b = row[f"protein2_{identifier}"]
                    interaction_id = utils.canonical_interaction_id(id_a, id_b)

                    yield {
                        "interaction_id": interaction_id,
                        "protein1": id_a,
                        "protein2": id_b,
                        "score": "",
                    }


def process_dscript_uniref(processed_folder: Path, taxon: int, threshold: int, console):
    train_paths = (
        processed_folder / f"dscript_{taxon}_train_upkb.csv.gz",
        processed_folder / f"dscript_{taxon}_train_uniref{threshold}.csv.gz",
    )
    test_paths = (
        processed_folder / f"dscript_{taxon}_test_upkb.csv.gz",
        processed_folder / f"dscript_{taxon}_test_uniref{threshold}.csv.gz",
    )

    uniref_db = plyvel.DB(str(processed_folder / f"uniref{threshold}_members_upkb.leveldb"))

    fieldnames = [
        "protein1_string",
        "protein1_upkb",
        f"protein1_uniref{threshold}",
        "protein2_string",
        "protein2_upkb",
        f"protein2_uniref{threshold}",
        "interaction",
    ]

    missed_ids = 0
    hit_ids = 0

    for path_in, path_out in [train_paths, test_paths]:
        f_in = gzip.open(path_in, "rt")
        f_out = gzip.open(path_out, "wt")

        csv_in, csv_out = (
            DictReader(f_in),
            DictWriter(f_out, fieldnames),
        )
        csv_out.writeheader()

        for row in csv_in:
            if len(row) < 5:
                console.log(f'Non-conforming line "{row}"')
                continue

            protein1_upkb = row["protein1_upkb"]
            protein2_upkb = row["protein2_upkb"]
            interaction = row["interaction"]

            protein1_uniref = uniref_db.get(protein1_upkb.encode("utf8"))
            protein2_uniref = uniref_db.get(protein2_upkb.encode("utf8"))

            if protein1_uniref is None or protein2_uniref is None:
                missed_ids += 1
                continue
            else:
                hit_ids += 1

            protein1_uniref = protein1_uniref.decode("utf8")
            protein2_uniref = protein2_uniref.decode("utf8")

            csv_out.writerow(
                {
                    "protein1_string": row["protein1_string"],
                    "protein2_string": row["protein1_string"],
                    "protein1_upkb": row["protein1_upkb"],
                    "protein2_upkb": row["protein2_upkb"],
                    f"protein1_uniref{threshold}": protein1_uniref,
                    f"protein2_uniref{threshold}": protein2_uniref,
                    "interaction": interaction,
                }
            )

    if missed_ids > 0:
        console.log(
            f"WARNING: {missed_ids}/{(hit_ids + missed_ids)} ({100 * missed_ids / (hit_ids + missed_ids):.3}%) ids don't"
            "have corresponding UniRef ids"
        )


def process_dscript_upkb(
    raw_folder: Path, processed_folder: Path, taxon: int, console
):

    aliases = update_secondaries(raw_folder, processed_folder, {}, console)

    train_paths = (
        raw_folder / f"dscript_{taxon}_train.tsv",
        processed_folder / f"dscript_{taxon}_train_upkb.csv.gz",
    )
    test_paths = (
        raw_folder / f"dscript_{taxon}_test.tsv",
        processed_folder / f"dscript_{taxon}_test_upkb.csv.gz",
    )

    if not os.path.isfile(processed_folder / "ensembl_pro_upkb.leveldb"):
        console.log("Building ensembl/upkb database...")
        uniprot.build_ensembl_pro_to_upkb_db(processed_folder)

    db = plyvel.DB(str(processed_folder / "ensembl_pro_upkb.leveldb"))

    fieldnames = [
        "protein1_string",
        "protein1_upkb",
        "protein2_string",
        "protein2_upkb",
        "interaction",
    ]

    missed_ids = 0
    hit_ids = 0

    for path_in, path_out in [train_paths, test_paths]:
        f_in = open(path_in)
        f_out = gzip.open(path_out, "wt")

        csv_in, csv_out = (
            csv_reader(f_in, delimiter="\t"),
            DictWriter(f_out, fieldnames),
        )
        csv_out.writeheader()

        for row in csv_in:
            if len(row) < 3:
                console.log(f'Non-conforming line "{row}"')
                continue

            protein1_string = row[0]
            protein2_string = row[1]
            interaction = row[2]

            try:
                protein1_upkb = aliases[protein1_string]
            except KeyError:
                try:
                    protein1_upkb = db.get(protein1_string.split(".")[0].encode("utf8"))

                    if protein1_upkb is None:
                        missed_ids += 1
                        continue
                    else:
                        protein1_upkb = protein1_upkb.decode("utf8")

                except KeyError:
                    missed_ids += 1
                    continue

            try:
                protein2_upkb = aliases[protein2_string]
            except KeyError:
                try:
                    protein2_upkb = db.get(protein2_string.split(".")[0].encode("utf8"))

                    if protein2_upkb is None:
                        missed_ids += 1
                        continue
                    else:
                        protein2_upkb = protein2_upkb.decode("utf8")

                except KeyError:
                    missed_ids += 1
                    continue

            hit_ids += 1

            csv_out.writerow(
                {
                    "protein1_string": protein1_string,
                    "protein2_string": protein2_string,
                    "protein1_upkb": protein1_upkb,
                    "protein2_upkb": protein2_upkb,
                    "interaction": interaction,
                }
            )

    if missed_ids > 0:
        console.log(
            f"WARNING: {missed_ids}/{(hit_ids + missed_ids)} ({100 * missed_ids / (hit_ids + missed_ids):.3}%) ids don't have corresponding UPKB ids"
        )


def common_to_dscript(
    console,
    processed_folder: Path,
    common_path: Path,
    c_types: List[int],
    train_proportion: float,
    val_proportion: float,
    test_proportion: float,
    neg_proportion: float,
    uniref_threshold: int,
    score_key: Optional[str] = None,
    score_threshold: Optional[float] = None,
    preloaded_protein_splits_path: Optional[Path] = None,
    seed: int = 8675309,
    trim_unseen_proteins: bool = True,
    negatives_path: Optional[Path] = None,
    taxon: Optional[int] = None,
    weighted_random: bool = False,
):
    """
    Convert a dataset in the common format to the D-SCRIPT format.

    :param processed_folder: The processed folder, where you will deposit the new dataset.
    :param common_path: The path to the common file.
    :param c_types: The different Park & Marcotte C-type levels to generate. Takes a list. e.g.: [1,2]
    :param train_proportion: The proportion of interactions to assign to the training fold.
    :param val_proportion: The proportion of interactions to assign to the validation fold.
    :param test_proportion: The proportion of interactions to assign to the testing fold.
    :param neg_proportion: The proportion of interactions that will be negative interactions.
    :param uniref_threshold: The UniRef threshold rate to use to ensure proteins between splits are not too similar.
    :param score_key: The scoring key to theshold by, if any.
    :param score_threshold: The value to threshold the score key by. Values below this value will be filtered out.
    :param preloaded_protein_splits_path: Load protein splits from another RAPPPID dataset.
    :param seed: An integer that will serve as the random seed for datasets.
    :param trim_unseen_proteins: If true, when a protein loaded from the preloaded_protein_splits is not found in the common dataset, it is not included in the dataset.
    :param negatives_path: Optional, path to file with negative interactions.
    :param taxon: Optional, restrict a dataset to a certain organism.
    :param weighted_random: If true, negative samples will be sampled in such a way as to maintain the same protein degree as the positive samples.
    :return: None
    """
    interactions, proteins = utils.c_ify(
        console,
        common_path,
        processed_folder,
        c_types,
        train_proportion,
        val_proportion,
        test_proportion,
        neg_proportion,
        uniref_threshold,
        score_key,
        score_threshold,
        preloaded_protein_splits_path,
        seed,
        None,
        trim_unseen_proteins,
        negatives_path,
        taxon,
        weighted_random,
    )


    metadata = {
        "common_path": common_path.stem,
        "c_types": c_types,
        "train_proportion": train_proportion,
        "val_proportion": val_proportion,
        "test_proportion": test_proportion,
        "neg_proportion": neg_proportion,
        "score_key": score_key,
        "score_threshold": score_threshold,
        "preloaded_protein_splits_path": None,
        "seed": seed,
        "trim_unseen_proteins": trim_unseen_proteins,
        "negatives_path": negatives_path,
    }

    dataset_hash = utils.hash_dict(metadata)

    dataset_folder = Path(
        processed_folder / f"dscript_[{common_path.stem}]_{dataset_hash}"
    )

    os.makedirs(dataset_folder, exist_ok=True)

    with open(dataset_folder / "metadata.json", "w") as f:
        json.dump(metadata, f, indent=4)

    for split in ["train", "test", "val"]:
        for c_type in ["C1", "C2", "C3"]:
            os.makedirs(dataset_folder / f"{c_type}", exist_ok=True)

            with open(dataset_folder / f"{c_type}/{split}.tsv", "w") as f:
                writer = csv.writer(f, delimiter="\t", quoting=csv.QUOTE_MINIMAL)

                for protein1, protein2, label in interactions[c_type][split]:
                    writer.writerow([protein1, protein2, 1 if label else 0])


def download(path: Path, taxon: int):
    if taxon == 9606:
        species = "human"
    elif taxon == 10090:
        species = "mouse"
    elif taxon == 7227:
        species = "fly"
    elif taxon == 4932:
        species = "yeast"
    elif taxon == 6239:
        species = "worm"
    elif taxon == 511145:
        species = "ecoli"
    else:
        raise ValueError(f"Invalid taxon value {taxon}.")

    url = f"https://raw.githubusercontent.com/samsledje/D-SCRIPT/main/data/pairs/{species}_train.tsv"
    utils.download_file(
        url,
        path / f"dscript_{taxon}_train.tsv",
        f"D-SCRIPT Training Set (Taxon: {taxon})",
    )

    url = f"https://raw.githubusercontent.com/samsledje/D-SCRIPT/main/data/pairs/{species}_test.tsv"
    utils.download_file(
        url,
        path / f"dscript_{taxon}_test.tsv",
        f"D-SCRIPT Testing Set (Taxon: {taxon})",
    )

    url = f"https://raw.githubusercontent.com/samsledje/D-SCRIPT/main/data/seqs/{species}.fasta"
    utils.download_file(
        url,
        path / f"dscript_{taxon}.fasta",
        f"D-SCRIPT Sequences (Taxon: {taxon})",
    )
