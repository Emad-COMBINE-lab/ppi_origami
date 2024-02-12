import gzip
import json
import math
import tempfile
from base64 import urlsafe_b64encode
from collections import defaultdict
from hashlib import blake2b
from pathlib import Path
from random import shuffle
from typing import List, Optional

import plyvel
import sentencepiece as sp
import tables as tb
from rich.console import Console
from rich.progress import track, Progress

from .utils import hash_dict, c_ify, canonical_id_sort, get_random_negative, MaxIterReached


class Protein(tb.IsDescription):
    name = tb.StringCol(20, pos=1)  # UniProt/Ref AC/ID
    sequence = tb.StringCol(3000, pos=2)  # Amino acid sequence


class Pair(tb.IsDescription):
    protein_id1 = tb.StringCol(20, pos=1)  # UniProt/Ref AC/ID
    protein_id2 = tb.StringCol(20, pos=2)  # Amino acid sequence
    label = tb.BoolCol(pos=3)


def write_rapppid2(
    console, interactions, proteins, processed_folder, outpath, c_types, scramble=False
):
    if scramble:
        console.log(
            "[red bold]WARNING: Scrambling activated!!! Proteins are randomly associated with sequences of other proteins."
        )

    filters = tb.Filters(
        complevel=9, complib="blosc:zstd", fletcher32=True, bitshuffle=True
    )
    h5file = tb.open_file(outpath, mode="w", filters=filters)

    group_interactions = h5file.create_group(
        "/", "interactions", "Protein Interactions"
    )

    group = {}
    group["C3"] = h5file.create_group(group_interactions, "c3", "C3 Datasets")
    group["C2"] = h5file.create_group(group_interactions, "c2", "C2 Datasets")
    group["C1"] = h5file.create_group(group_interactions, "c1", "C1 Datasets")

    table_pairs = {f"C{c}": {} for c in c_types}

    for c in c_types:
        for split in ["train", "test", "val"]:
            table_pairs[f"C{c}"][split] = h5file.create_table(
                group[f"C{c}"],
                f"c{c}_{split}",
                Pair,
                f"C{c} {split.title()} Dataset",
                filters=filters,
            )

    group_splits = h5file.create_group("/", "splits", "Protein Splits")

    train_proteins = h5file.create_carray(
        group_splits,
        "train",
        tb.StringCol(20),
        (len(proteins["train"]),),
        filters=filters,
    )
    val_proteins = h5file.create_carray(
        group_splits, "val", tb.StringCol(20), (len(proteins["val"]),), filters=filters
    )
    test_proteins = h5file.create_carray(
        group_splits,
        "test",
        tb.StringCol(20),
        (len(proteins["test"]),),
        filters=filters,
    )

    table_sequences = h5file.create_table(
        "/", "sequences", Protein, "Protein Sequences", filters=filters
    )

    for c_type in c_types:
        for split in ["train", "val", "test"]:
            for p1, p2, label in track(
                interactions[f"C{c_type}"][split], "Saving pairs..."
            ):
                pair = table_pairs[f"C{c_type}"][split].row
                pair["protein_id1"] = p1
                pair["protein_id2"] = p2
                pair["label"] = label
                pair.append()

    train_proteins[:] = list(proteins["train"])
    val_proteins[:] = list(proteins["val"])
    test_proteins[:] = list(proteins["test"])

    all_proteins = list(set(proteins["train"] + proteins["val"] + proteins["test"]))

    console.log("Saving sequences")

    db = plyvel.DB(
        str(processed_folder / "uniprot_sequences.leveldb"), create_if_missing=False
    )

    if scramble:
        scrambled_proteins = all_proteins.copy()
        shuffle(scrambled_proteins)

        protein_scramble_map = {
            protein: scrambled_protein
            for protein, scrambled_protein in zip(all_proteins, scrambled_proteins)
        }

    for target_protein in track(all_proteins):
        sequence = table_sequences.row

        sequence["name"] = target_protein

        if scramble:
            target_protein = protein_scramble_map[target_protein]
        try:
            return_seq = db.get(target_protein.encode("utf8"))
        except Exception:
            print(f"Error finding sequence {target_protein}")
            return_seq = None
        if return_seq is None:
            print(f"Could not find sequence for {target_protein}")
        else:
            sequence["sequence"] = return_seq.decode("utf8")

        sequence.append()

    return h5file


def rapppid_to_sprint(
    dataset_path,
    c_type,
    seq_path,
    train_pos_path,
    train_neg_path,
    train_seq_path,
    val_pos_path,
    val_neg_path,
    test_pos_path,
    test_neg_path,
    metadata_path,
):
    rapppid_dataset = tb.open_file(dataset_path)

    metadata = rapppid_dataset.root.metadata
    cols = metadata.colnames
    rows = [row[:] for row in metadata.iterrows()]

    metadata_dict = defaultdict(lambda: [])

    for col_idx, col in enumerate(cols):
        for row in rows:
            try:
                value = row[col_idx].decode("utf8")
            except AttributeError:
                value = row[col_idx]

            metadata_dict[col].append(value)

    metadata_dict["rapppid_file"] = dataset_path

    with open(metadata_path, "w") as f:
        json.dump(metadata_dict, f, indent=4)

    seqs = {}

    with open(seq_path, "w") as f:
        for row in rapppid_dataset.root.sequences:
            name = row["name"].decode("utf8")
            sequence = row["sequence"].decode("utf8")

            f.write(f">{name}\n")
            f.write(f"{sequence}\n")
            f.write("\n")

            seqs[name] = sequence

    with open(train_seq_path, "w") as f_seq:
        with open(train_pos_path, "w") as f_pos:
            with open(train_neg_path, "w") as f_neg:
                if c_type == 3:
                    c_table = rapppid_dataset.root.interactions.c3.c3_train
                elif c_type == 2:
                    c_table = rapppid_dataset.root.interactions.c2.c2_train
                elif c_type == 1:
                    c_table = rapppid_dataset.root.interactions.c1.c1_train
                else:
                    raise ValueError("Unexpected value for c_type")

                for row in c_table:
                    p1 = row["protein_id1"].decode("utf8")
                    p2 = row["protein_id2"].decode("utf8")

                    if row["label"] is True:
                        f_pos.write(f"{p1} {p2}\n")
                    else:
                        f_neg.write(f"{p1} {p2}\n")

                    for name in [p1, p2]:
                        f_seq.write(f">{name}\n")
                        f_seq.write(f"{seqs[name]}\n")
                        f_seq.write("\n")

    with open(test_pos_path, "w") as f_pos:
        with open(test_neg_path, "w") as f_neg:
            if c_type == 3:
                c_table = rapppid_dataset.root.interactions.c3.c3_test
            elif c_type == 2:
                c_table = rapppid_dataset.root.interactions.c2.c2_test
            elif c_type == 1:
                c_table = rapppid_dataset.root.interactions.c1.c1_test
            else:
                raise ValueError("Unexpected value for c_type")

            for row in c_table:
                if row["label"] is True:
                    f_pos.write(
                        f"{row['protein_id1'].decode('utf8')} {row['protein_id2'].decode('utf8')}\n"
                    )
                else:
                    f_neg.write(
                        f"{row['protein_id1'].decode('utf8')} {row['protein_id2'].decode('utf8')}\n"
                    )

    with open(val_pos_path, "w") as f_pos:
        with open(val_neg_path, "w") as f_neg:
            if c_type == 3:
                c_table = rapppid_dataset.root.interactions.c3.c3_val
            elif c_type == 2:
                c_table = rapppid_dataset.root.interactions.c2.c2_val
            elif c_type == 1:
                c_table = rapppid_dataset.root.interactions.c1.c1_val
            else:
                raise ValueError("Unexpected value for c_type")

            for row in c_table:
                if row["label"] is True:
                    f_pos.write(
                        f"{row['protein_id1'].decode('utf8')} {row['protein_id2'].decode('utf8')}\n"
                    )
                else:
                    f_neg.write(
                        f"{row['protein_id1'].decode('utf8')} {row['protein_id2'].decode('utf8')}\n"
                    )


def rapppid_to_dscript(
    dataset_path, c_type, seq_path, train_path, test_path, val_path, metadata_path, trunc_len
):
    rapppid_dataset = tb.open_file(dataset_path)

    metadata = rapppid_dataset.root.metadata
    cols = metadata.colnames
    rows = [row[:] for row in metadata.iterrows()]

    metadata_dict = defaultdict(lambda: [])

    for col_idx, col in enumerate(cols):
        for row in rows:
            try:
                value = row[col_idx].decode("utf8")
            except AttributeError:
                value = row[col_idx]

            metadata_dict[col].append(value)

    with open(metadata_path, "w") as f:
        json.dump(metadata_dict, f, indent=4)

    with open(seq_path, "w") as f:
        for row in rapppid_dataset.root.sequences:
            f.write(f">{row['name'].decode('utf8')}\n")
            f.write(f"{row['sequence'].decode('utf8')[:trunc_len]}\n")
            f.write("\n")

    with open(train_path, "w") as f:
        if c_type == 3:
            c_table = rapppid_dataset.root.interactions.c3.c3_train
        elif c_type == 2:
            c_table = rapppid_dataset.root.interactions.c2.c2_train
        elif c_type == 1:
            c_table = rapppid_dataset.root.interactions.c1.c1_train
        else:
            raise ValueError("Unexpected value for c_type")

        for row in c_table:
            if row["label"] is True:
                label = 1
            else:
                label = 0

            f.write(
                f"{row['protein_id1'].decode('utf8')}\t{row['protein_id2'].decode('utf8')}\t{label}\n"
            )

    with open(test_path, "w") as f:
        if c_type == 3:
            c_table = rapppid_dataset.root.interactions.c3.c3_test
        elif c_type == 2:
            c_table = rapppid_dataset.root.interactions.c2.c2_test
        elif c_type == 1:
            c_table = rapppid_dataset.root.interactions.c1.c1_test
        else:
            raise ValueError("Unexpected value for c_type")

        for row in c_table:
            if row["label"] is True:
                label = 1
            else:
                label = 0

            f.write(
                f"{row['protein_id1'].decode('utf8')}\t{row['protein_id2'].decode('utf8')}\t{label}\n"
            )

    with open(val_path, "w") as f:
        if c_type == 3:
            c_table = rapppid_dataset.root.interactions.c3.c3_val
        elif c_type == 2:
            c_table = rapppid_dataset.root.interactions.c2.c2_val
        elif c_type == 1:
            c_table = rapppid_dataset.root.interactions.c1.c1_val
        else:
            raise ValueError("Unexpected value for c_type")

        for row in c_table:
            if row["label"] is True:
                label = 1
            else:
                label = 0
            f.write(
                f"{row['protein_id1'].decode('utf8')}\t{row['protein_id2'].decode('utf8')}\t{label}\n"
            )


def common_to_rapppid(
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
    scramble_proteins: bool = False,
    exclude_preloaded_from_neg: bool = True
):
    """
    Convert a dataset in the common format to the RAPPPID format.

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
    :param scramble_proteins: Scramble the association between protein ids and their sequences.
    :param exclude_preloaded_from_neg: Set this to true if you don't want preloaded proteins to leak into negative.
    :return: None
    """
    class Metadata(tb.IsDescription):
        common_path = tb.StringCol(560, pos=1)
        seed = tb.Int64Col(pos=2)
        c1 = tb.BoolCol(pos=3)
        c2 = tb.BoolCol(pos=4)
        c3 = tb.BoolCol(pos=5)
        train_proportion = tb.FloatCol(pos=6)
        val_proportion = tb.FloatCol(pos=7)
        test_proportion = tb.FloatCol(pos=8)
        neg_proportion = tb.FloatCol(pos=9)
        score_key = tb.StringCol(280, pos=10)
        score_threshold = tb.Float64Col(pos=11)
        preloaded_protein_splits_path = tb.StringCol(560, pos=12)
        trim_unseen_proteins = tb.BoolCol(pos=13)
        negatives_path = tb.StringCol(560, pos=14)
        taxon = tb.Int64Col(pos=15)
        weighted_random = tb.BoolCol(pos=16)
        scramble = tb.BoolCol(pos=17)
        exclude_preloaded_from_neg = tb.BoolCol(pos=18)

    if preloaded_protein_splits_path is not None:
        console.log("Preloading protein splits...")

    if (
        preloaded_protein_splits_path is not None
        and str(preloaded_protein_splits_path).endswith(".json.gz")
    ):
        with gzip.open(preloaded_protein_splits_path, "rt") as f:
            preloaded_protein_splits = json.load(f)

    elif (
        preloaded_protein_splits_path is not None
        and str(preloaded_protein_splits_path).endswith(".h5")
    ):
        dataset = tb.open_file(str(preloaded_protein_splits_path))

        table = {
            "val_proteins": dataset.root["splits"].val,
            "train_proteins": dataset.root["splits"].train,
            "test_proteins": dataset.root["splits"].test,
        }

        preloaded_protein_splits = {
            "val_proteins": [],
            "train_proteins": [],
            "test_proteins": [],
        }

        for split in ["val_proteins", "train_proteins", "test_proteins"]:
            for row in table[split].iterrows():
                preloaded_protein_splits[split].append(row.decode("utf8"))
    else:
        preloaded_protein_splits = None

    interactions, proteins = c_ify(
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
        preloaded_protein_splits,
        seed,
        trim_unseen_proteins=trim_unseen_proteins,
        negatives_path=negatives_path,
        taxon=taxon,
        weighted_random=weighted_random,
        exclude_preloaded_from_neg=exclude_preloaded_from_neg
    )

    for c in [1, 2, 3]:
        for split in ['train', 'val', 'test']:
            console.log(f"C{c} {split} edges: {len(interactions['C'+str(c)][split])}")

    for split in ['train', 'val', 'test']:
        console.log(f"{split} proteins: {len(proteins[split])}")

    dataset_hash = hash_dict(
        {
            "common_path": common_path.stem,
            "c_types": c_types,
            "train_proportion": train_proportion,
            "val_proportion": val_proportion,
            "test_proportion": test_proportion,
            "neg_proportion": neg_proportion,
            "score_key": score_key,
            "score_threshold": score_threshold,
            "preloaded_protein_splits_path": preloaded_protein_splits_path,
            "seed": seed,
            "trim_unseen_proteins": trim_unseen_proteins,
            "negatives_path": negatives_path,
            "taxon": taxon,
            "weighted_random": weighted_random,
            "scramble": scramble_proteins,
            "exclude_preloaded_from_neg": exclude_preloaded_from_neg
        }
    )

    outpath = processed_folder / f"rapppid_[{common_path.stem}]_{dataset_hash}.h5"

    h5file = write_rapppid2(
        console, interactions, proteins, processed_folder, outpath, c_types, scramble_proteins
    )

    table_metadata = h5file.create_table("/", "metadata", Metadata, "Metadata")

    metadata = table_metadata.row
    metadata["common_path"] = common_path.stem
    metadata["seed"] = seed
    metadata["c1"] = 1 in c_types
    metadata["c2"] = 2 in c_types
    metadata["c3"] = 3 in c_types
    metadata["train_proportion"] = train_proportion
    metadata["val_proportion"] = val_proportion
    metadata["test_proportion"] = test_proportion
    metadata["neg_proportion"] = neg_proportion
    metadata["score_key"] = score_key
    metadata["score_threshold"] = score_threshold
    metadata["preloaded_protein_splits_path"] = preloaded_protein_splits_path
    metadata["trim_unseen_proteins"] = trim_unseen_proteins
    metadata["negatives_path"] = negatives_path if negatives_path is not None else ""
    metadata["taxon"] = taxon if taxon is not None else -1
    metadata["weighted_random"] = weighted_random
    metadata["scramble"] = scramble_proteins
    metadata["exclude_preloaded_from_neg"] = exclude_preloaded_from_neg
    metadata.append()


def merge_rapppid(
    console,
    processed_folder: Path,
    dataset_path1: Path,
    dataset_path2: Path,
    scramble: bool = False,
):
    class Metadata(tb.IsDescription):
        dataset_num = tb.Int64Col(pos=1)
        dataset_path = tb.StringCol(560, pos=2)
        common_path = tb.StringCol(560, pos=3)
        seed = tb.Int64Col(pos=4)
        c1 = tb.BoolCol(pos=5)
        c2 = tb.BoolCol(pos=6)
        c3 = tb.BoolCol(pos=7)
        train_proportion = tb.FloatCol(pos=8)
        val_proportion = tb.FloatCol(pos=9)
        test_proportion = tb.FloatCol(pos=10)
        neg_proportion = tb.FloatCol(pos=11)
        score_key = tb.StringCol(280, pos=12)
        score_threshold = tb.Float64Col(pos=13)
        preloaded_protein_splits_path = tb.StringCol(560, pos=14)

    dataset1 = tb.open_file(str(dataset_path1))
    dataset2 = tb.open_file(str(dataset_path2))

    c1_1 = [x["c1"] for x in dataset1.root["metadata"].iterrows()][0]
    c1_2 = [x["c1"] for x in dataset2.root["metadata"].iterrows()][0]
    c2_1 = [x["c2"] for x in dataset1.root["metadata"].iterrows()][0]
    c2_2 = [x["c2"] for x in dataset2.root["metadata"].iterrows()][0]
    c3_1 = [x["c3"] for x in dataset1.root["metadata"].iterrows()][0]
    c3_2 = [x["c3"] for x in dataset2.root["metadata"].iterrows()][0]

    if c1_1 != c1_2 or c2_1 != c2_2 or c3_1 != c3_2:
        console.log(
            f'Dataset 1: {"C1 " if c1_1 else ""}{"C2 " if c2_1 else ""}{"C3" if c3_1 else ""}'
        )
        console.log(
            f'Dataset 2: {"C1 " if c1_2 else ""}{"C2 " if c2_2 else ""}{"C3" if c3_2 else ""}'
        )
        raise ValueError("Both datasets must have the same C values")

    c_types = []

    if c1_1:
        c_types.append(1)

    if c2_1:
        c_types.append(2)

    if c3_1:
        c_types.append(3)

    interactions = {
        f"C{c}": {split: [] for split in [f"train", f"val", f"test"]} for c in [1, 2, 3]
    }
    proteins = {split: [] for split in ["train", "val", "test"]}
    labels = {f"C{c}": [] for c in [1, 2, 3]}

    for dataset_num, dataset in enumerate([dataset1, dataset2]):
        for split in ["train", "val", "test"]:
            # interactions
            for c in [1, 2, 3]:
                intxn_table = dataset.root["interactions"][f"c{c}"][f"c{c}_{split}"]

                for row in intxn_table.iterrows():
                    p1 = row[0].decode("utf8")
                    p2 = row[1].decode("utf8")
                    label = row[2]
                    interactions[f"C{c}"][split].append((p1, p2, label))
                    if split == "test":
                        labels[f"C{c}"].append(dataset_num)

            # sequences
            protein_tables = dataset.root["splits"][split]
            for row in protein_tables.iterrows():
                proteins[split].append(row.decode("utf8"))

    for c in [1, 2, 3]:
        shuffle(interactions[f"C{c}"]["train"])

    merged_filename = f"rapppid_merged_[{dataset_path1.stem}]X[{dataset_path2.stem}].h5"

    merged_dataset = write_rapppid2(
        console,
        interactions,
        proteins,
        processed_folder,
        processed_folder / merged_filename,
        c_types,
        scramble,
    )

    table_metadata = merged_dataset.create_table("/", "metadata", Metadata, "Metadata")

    for dataset_num, (dataset, dataset_path) in enumerate(
        [(dataset1, dataset_path1), (dataset2, dataset_path2)]
    ):
        metadata = table_metadata.row
        metadata["dataset_num"] = dataset_num
        metadata["dataset_path"] = dataset_path
        metadata["common_path"] = [
            x["common_path"] for x in dataset.root["metadata"].iterrows()
        ][0]
        metadata["seed"] = [x["seed"] for x in dataset.root["metadata"].iterrows()][0]
        metadata["c1"] = [x["c1"] for x in dataset.root["metadata"].iterrows()][0]
        metadata["c2"] = [x["c2"] for x in dataset.root["metadata"].iterrows()][0]
        metadata["c3"] = [x["c3"] for x in dataset.root["metadata"].iterrows()][0]
        metadata["train_proportion"] = [
            x["train_proportion"] for x in dataset.root["metadata"].iterrows()
        ][0]
        metadata["test_proportion"] = [
            x["test_proportion"] for x in dataset.root["metadata"].iterrows()
        ][0]
        metadata["val_proportion"] = [
            x["val_proportion"] for x in dataset.root["metadata"].iterrows()
        ][0]
        metadata["neg_proportion"] = [
            x["neg_proportion"] for x in dataset.root["metadata"].iterrows()
        ][0]


def multimerge_rapppid(
    console, processed_folder: Path, dataset_paths: List[Path], scramble: bool = False
):
    class Metadata(tb.IsDescription):
        dataset_num = tb.Int64Col(pos=1)
        dataset_path = tb.StringCol(560, pos=2)
        common_path = tb.StringCol(560, pos=3)
        seed = tb.Int64Col(pos=4)
        c1 = tb.BoolCol(pos=5)
        c2 = tb.BoolCol(pos=6)
        c3 = tb.BoolCol(pos=7)
        train_proportion = tb.FloatCol(pos=8)
        val_proportion = tb.FloatCol(pos=9)
        test_proportion = tb.FloatCol(pos=10)
        neg_proportion = tb.FloatCol(pos=11)
        score_key = tb.StringCol(280, pos=12)
        score_threshold = tb.Float64Col(pos=13)
        preloaded_protein_splits_path = tb.StringCol(560, pos=14)

    class TestDatasetLabels(tb.IsDescription):
        dataset_num = tb.Int64Col(pos=1)

    datasets = []

    cx = {"c1": [], "c2": [], "c3": []}

    for dataset_path in dataset_paths:
        dataset = tb.open_file(str(dataset_path))
        datasets.append(dataset)
        for c in ["c1", "c2", "c3"]:
            cx[c].append([x[c] for x in dataset.root["metadata"].iterrows()][0])

    c_types = []

    for c in [1, 2, 3]:
        # All datasets must have the same value (all true, or all false)
        if sum(cx[f"c{c}"]) != len(cx[f"c{c}"]) and sum(cx[f"c{c}"]) != 0:
            console.log(cx)
            raise ValueError("Both datasets must have the same C values")

        if cx[f"c{c}"][0] is True:
            c_types.append(c)

    interactions = {
        f"C{c}": {split: [] for split in [f"train", f"val", f"test"]} for c in [1, 2, 3]
    }
    proteins = {split: [] for split in ["train", "val", "test"]}
    labels = {f"C{c}": [] for c in [1, 2, 3]}

    for dataset_num, dataset in enumerate(datasets):
        for split in ["train", "val", "test"]:
            # interactions
            for c in [1, 2, 3]:
                intxn_table = dataset.root["interactions"][f"c{c}"][f"c{c}_{split}"]

                for row in intxn_table.iterrows():
                    p1 = row[0].decode("utf8")
                    p2 = row[1].decode("utf8")
                    label = row[2]
                    interactions[f"C{c}"][split].append((p1, p2, label))
                    if split == "test":
                        labels[f"C{c}"].append(dataset_num)

            # sequences
            protein_tables = dataset.root["splits"][split]
            for row in protein_tables.iterrows():
                proteins[split].append(row.decode("utf8"))

    for c in [1, 2, 3]:
        shuffle(interactions[f"C{c}"]["train"])

    hasher = blake2b(digest_size=40)

    for seq_file in sorted(dataset_paths):
        hasher.update(str(seq_file).encode("utf8"))

    merged_filename = (
        f"rapppid_merged_{urlsafe_b64encode(hasher.digest()).decode('utf8')}.h5"
    )

    merged_dataset = write_rapppid2(
        console,
        interactions,
        proteins,
        processed_folder,
        processed_folder / merged_filename,
        c_types,
        scramble,
    )

    table_metadata = merged_dataset.create_table("/", "metadata", Metadata, "Metadata")

    for dataset_num, (dataset, dataset_path) in enumerate(zip(datasets, dataset_paths)):
        metadata = table_metadata.row
        metadata["dataset_num"] = dataset_num
        metadata["dataset_path"] = dataset_path
        metadata["common_path"] = [
            x["common_path"] for x in dataset.root["metadata"].iterrows()
        ][0]
        metadata["seed"] = [x["seed"] for x in dataset.root["metadata"].iterrows()][0]
        metadata["c1"] = [x["c1"] for x in dataset.root["metadata"].iterrows()][0]
        metadata["c2"] = [x["c2"] for x in dataset.root["metadata"].iterrows()][0]
        metadata["c3"] = [x["c3"] for x in dataset.root["metadata"].iterrows()][0]
        metadata["train_proportion"] = [
            x["train_proportion"] for x in dataset.root["metadata"].iterrows()
        ][0]
        metadata["test_proportion"] = [
            x["test_proportion"] for x in dataset.root["metadata"].iterrows()
        ][0]
        metadata["val_proportion"] = [
            x["val_proportion"] for x in dataset.root["metadata"].iterrows()
        ][0]
        metadata["neg_proportion"] = [
            x["neg_proportion"] for x in dataset.root["metadata"].iterrows()
        ][0]
        metadata["score_key"] = [
            x["score_key"] for x in dataset.root["metadata"].iterrows()
        ][0]
        metadata["score_threshold"] = [
            x["score_threshold"] for x in dataset.root["metadata"].iterrows()
        ][0]
        try:
            metadata["preloaded_protein_splits_path"] = [
                x["preloaded_protein_splits_path"]
                for x in dataset.root["metadata"].iterrows()
            ][0]
        except KeyError:
            metadata["preloaded_protein_splits_path"] = ""
        metadata.append()

    for c in [1, 2, 3]:
        table_labels = merged_dataset.create_table(
            merged_dataset.root["interactions"][f"c{c}"],
            "test_dataset_labels",
            TestDatasetLabels,
            "Test Dataset Labels",
        )
        for label in labels[f"C{c}"]:
            row = table_labels.row
            row["dataset_num"] = label
            row.append()


def train_sentencepiece_model(
    processed_folder: Path, dataset_path: Path, seed: int, vocab_size: int
):
    console = Console()

    console.log("Opening dataset file")
    sp.set_random_generator_seed(seed)
    dataset = tb.open_file(str(dataset_path))

    with tempfile.TemporaryDirectory() as tmpdirname:
        with open(Path(tmpdirname) / "seqs.txt", "w") as f:
            for row in track(
                dataset.root.splits.train.iterrows(),
                total=len(dataset.root.splits.train),
            ):
                sequence = dataset.root.sequences.read_where(
                    f'name=="{row.decode("utf8")}"'
                )[0][1].decode("utf8")
                f.write(sequence + "\n")

        sp.SentencePieceTrainer.train(
            input=Path(tmpdirname) / "seqs.txt",
            model_prefix=Path(dataset_path).stem,
            vocab_size=vocab_size,
            character_coverage=1.0,
            bos_id=1,
            eos_id=2,
            pad_id=0,
            unk_id=3,
        )


def inflate_eval_negatives(processed_folder: Path, dataset_path: Path, negative_ratio: float, c_types: List[int]):

    #c_types = ['C1', 'C2', 'C3']
    #c_types = [3]

    class Metadata(tb.IsDescription):
        common_path = tb.StringCol(560, pos=1)
        seed = tb.Int64Col(pos=2)
        c1 = tb.BoolCol(pos=3)
        c2 = tb.BoolCol(pos=4)
        c3 = tb.BoolCol(pos=5)
        train_proportion = tb.FloatCol(pos=6)
        val_proportion = tb.FloatCol(pos=7)
        test_proportion = tb.FloatCol(pos=8)
        neg_proportion = tb.FloatCol(pos=9)
        score_key = tb.StringCol(280, pos=10)
        score_threshold = tb.Float64Col(pos=11)
        preloaded_protein_splits_path = tb.StringCol(560, pos=12)
        trim_unseen_proteins = tb.BoolCol(pos=13)
        negatives_path = tb.StringCol(560, pos=14)
        taxon = tb.Int64Col(pos=15)
        weighted_random = tb.BoolCol(pos=16)
        scramble = tb.BoolCol(pos=17)
        exclude_preloaded_from_neg = tb.BoolCol(pos=18)

    console = Console()

    dataset = tb.open_file(str(dataset_path))

    metadata_table = dataset.root.metadata

    protein_table = {
        'train': dataset.root.splits.train,
        'test': dataset.root.splits.test,
        'val': dataset.root.splits.val
    }

    interaction_table = {
        'C3': {
            'train': dataset.root.interactions.c3.c3_train,
            'test': dataset.root.interactions.c3.c3_test,
            'val': dataset.root.interactions.c3.c3_val,
            },
        'C2': {
            'train': dataset.root.interactions.c2.c2_train,
            'test': dataset.root.interactions.c2.c2_test,
            'val': dataset.root.interactions.c2.c2_val,
        },
        'C1': {
            'train': dataset.root.interactions.c1.c1_train,
            'test': dataset.root.interactions.c1.c1_test,
            'val': dataset.root.interactions.c1.c1_val,
        }
    }

    proteins = {split: set() for split in ['train', 'test', 'val']}
    interactions = {c: {split: {label: set() for label in [1, 0]} for split in ['train', 'test', 'val']} for c in ['C1', 'C2', 'C3']}

    for split in ['train', 'test', 'val']:
        for _, row in enumerate(protein_table[split].iterrows()):
            proteins[split].add(row.decode('utf8'))

        for c in [f'C{c}' for c in c_types]:
            for _, row in enumerate(interaction_table[c][split].iterrows()):

                p1 = row['protein_id1'].decode('utf8')
                p2 = row['protein_id2'].decode('utf8')
                p1, p2 = canonical_id_sort(p1, p2)

                label = int(row['label'])

                interactions[c][split][label].add((p1, p2, label))

    total_negatives = 0

    for split in ['test', 'val']:
        for c in [f'C{c}' for c in c_types]:
            current_nr = len(interactions[c][split][0]) / len(interactions[c][split][1])

            if current_nr > negative_ratio:
                console.log(f"ERROR! - Current negative ratio ({current_nr}) is larger than target negative ratio ({negative_ratio})")

            diff_nr = negative_ratio - current_nr

            total_negatives += diff_nr * len(interactions[c][split][1])

    with Progress() as progress:
        task = progress.add_task("Draw Negatives...", total=total_negatives)

        for split in ['test', 'val']:
            for c in [f'C{c}' for c in c_types]:
                while len(interactions[c][split][0]) / len(interactions[c][split][1]) < negative_ratio:

                    try:
                        p1, p2 = get_random_negative(
                            proteins[split],
                            [],
                            interactions[c]['train'][1],
                            interactions[c]['val'][1],
                            interactions[c]['test'][1],
                            interactions[c]['train'][0],
                            interactions[c]['val'][0],
                            interactions[c]['test'][0],
                            False
                        )
                    except MaxIterReached:
                        console.log(f"[red][bold]MaxIterReached: {c}, {split}")
                        continue

                    p1, p2 = canonical_id_sort(p1, p2)
                    interactions[c][split][0].add((p1, p2, False))
                    progress.update(task, advance=1)

    print(len(interactions['C3']['test'][0]), len(interactions['C3']['test'][1]))

    for c in [f'C{c}' for c in c_types]:
        for split in ['train', 'test', 'val']:
            temp_interactions = list(interactions[c][split][0].union(interactions[c][split][1]))
            shuffle(temp_interactions)
            interactions[c][split] = temp_interactions

    for row in metadata_table.iterrows():
        metadata_row = row

    metadata_dict = {}

    fields = [
        ('common_path', 'path'),
        ('c_types', 'list'),
        ('train_proportion', 'float'),
        ('val_proportion', 'float'),
        ('test_proportion', 'float'),
        ('neg_proportion', 'float'),
        ('score_key', 'str'),
        ('score_threshold', 'int'),
        ('preloaded_protein_splits_path', 'path'),
        ('seed', 'int'),
        ('trim_unseen_proteins', 'bool'),
        ('negatives_path', 'path'),
        ('taxon', 'int'),
        ('weighted_random', 'bool'),
        ('scramble', 'bool'),
        ('exclude_preloaded_from_neg', 'bool')
    ]

    for field, field_type in fields:
        try:
            if field == 'neg_proportion':
                metadata_dict[field] = negative_ratio
            elif field == 'c_types':
                metadata_dict[field] = c_types
            else:
                metadata_dict[field] = metadata_row[field]
        except Exception:
            if field_type == 'path':
                metadata_dict[field] = ''
            elif field_type == 'list':
                metadata_dict[field] = []
            elif field_type == 'float':
                metadata_dict[field] = 0.0
            elif field_type == 'int':
                metadata_dict[field] = 0
            elif field_type == 'bool':
                metadata_dict[field] = False

    for split in ['train', 'val', 'test']:
        proteins[split] = list(proteins[split])

    dataset_hash = hash_dict(metadata_dict)

    outpath = processed_folder / f"rapppid_[{metadata_dict['common_path'].decode('utf8')}]_{dataset_hash}.h5"

    h5file = write_rapppid2(
            console, interactions, proteins, processed_folder, outpath, c_types, scramble=False
    )

    table_metadata = h5file.create_table("/", "metadata", Metadata, "Metadata")

    metadata = table_metadata.row

    for field, _ in fields:
        if field == 'c_types':
            metadata['c1'] = True if 1 in metadata_dict['c_types'] else False
            metadata['c2'] = True if 2 in metadata_dict['c_types'] else False
            metadata['c3'] = True if 3 in metadata_dict['c_types'] else False
        else:
            metadata[field] = metadata_dict[field]

    metadata.append()

