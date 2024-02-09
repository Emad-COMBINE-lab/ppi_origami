import csv
import gzip
import os
import time
from base64 import urlsafe_b64encode
from collections import defaultdict
from hashlib import sha1
from os.path import isfile
from pathlib import Path
from random import randint, random, sample, choices, shuffle
from random import seed as set_seed
from typing import Dict, Tuple, List, Optional, Generator, Set
import re
from urllib.parse import urlparse, parse_qs

import plyvel
import requests
from requests.adapters import HTTPAdapter, Retry
from rich.progress import Progress, track


def get_aa_code(aa):
    # Codes based on IUPAC-IUB
    # https://web.expasy.org/docs/userman.html#AA_codes

    aas = [
        "PAD",
        "A",
        "R",
        "N",
        "D",
        "C",
        "Q",
        "E",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
        "O",
        "U",
    ]
    wobble_aas = {
        "B": ["D", "N"],
        "Z": ["Q", "E"],
        "X": [
            "A",
            "R",
            "N",
            "D",
            "C",
            "Q",
            "E",
            "G",
            "H",
            "I",
            "L",
            "K",
            "M",
            "F",
            "P",
            "S",
            "T",
            "W",
            "Y",
            "V",
        ],
    }

    if aa in aas:
        return aas.index(aa)

    elif aa in ["B", "Z", "X"]:
        # Wobble
        idx = randint(0, len(wobble_aas[aa]) - 1)
        return aas.index(wobble_aas[aa][idx])


def encode_seq(seq):
    return [get_aa_code(aa) for aa in seq]


def mitab_deserialize(seialized_string: str) -> Dict[str, str]:
    return {
        tagged_value.split(":")[0]: tagged_value.split(":")[1]
        for tagged_value in seialized_string.split("|")
    }


def canonical_id_sort(id_a: str, id_b: str) -> Tuple[str, str]:
    """
    Deterministically sorts protein ids according to their SHA1 hash.

    :param id_a: A protein in the interaction.
    :param id_b: Another protein in the interaciton.
    :return: Sorted ids.
    """

    if int(sha1(str(id_a).encode("utf8")).hexdigest(), 16) > int(
        sha1(str(id_b).encode("utf8")).hexdigest(), 16
    ):
        return id_a, id_b
    else:
        return id_b, id_a


def canonical_interaction_id(id_a: str, id_b: str) -> str:
    """
    Provides a canonical id for an interaction between to proteins with identifiers `id_a` and `id_b`.

    :param id_a: A protein in the interaction.
    :param id_b: Another protein in the interaciton.
    :return: A canonical id for the interaction between these two proteins.
    """

    id_a = id_a.replace("><", "")
    id_b = id_b.replace("><", "")

    return "><".join(canonical_id_sort(id_a, id_b)) + ".1"


def serialize_dict(dictionary: dict) -> str:
    """
    Serialize a dictionary in a "key1:value1|key2:value2" format.
    :param dictionary: The dictionary to serialize.
    :return: A serialized string.
    """
    return "|".join([f"{key}:{dictionary[key]}" for key in dictionary])


def download_file(url: str, outfile: Path, description: str = "Downloading"):
    if isfile(outfile):
        raise IOError

    # Streaming, so we can iterate over the response.
    response = requests.get(url, stream=True)
    total_size_in_bytes = int(response.headers.get("content-length", 0))
    block_size = 1024  # 1 Kibibyte
    with Progress() as progress:
        progress_bar = progress.add_task(
            f"[green]{description}...", total=total_size_in_bytes
        )
        with open(outfile, "wb") as file:
            for data in response.iter_content(block_size):
                progress.update(progress_bar, advance=len(data))
                file.write(data)


def make_hashable(o):
    # from https://stackoverflow.com/questions/5884066/hashing-a-dictionary
    if isinstance(o, (tuple, list)):
        return tuple((make_hashable(e) for e in o))

    if isinstance(o, dict):
        return tuple(sorted((k, make_hashable(v)) for k, v in o.items()))

    if isinstance(o, (set, frozenset)):
        return tuple(sorted(make_hashable(e) for e in o))

    return o


def hash_dict(d: dict) -> str:
    """
    Create a portable, deterministic hash of a dictionary d
    :param d: dictionary to hash
    :return: The hash of the dict
    """
    return urlsafe_b64encode(sha1(repr(make_hashable(d)).encode()).digest()).decode()


def threeway_split(rate1: float, rate2: float) -> int:
    if rate1 < 0 or rate1 > 1 or rate2 < 0 or rate2 > 1:
        raise ValueError(
            "train_proportion and val_proportion need to be between 0 and 1"
        )

    if rate1 + rate2 > 1:
        raise ValueError(
            "The sum of train_proportion and val_proportion must be less than 1"
        )

    dice = random()

    if dice < rate1:
        return 0
    elif rate1 <= dice < rate1 + rate2:
        return 1
    else:
        return 2

def get_split(
    protein: str,
    train_proteins: set,
    val_proteins: set,
    test_proteins: set,
    train_interactions: set,
    val_interactions: set,
    test_interactions: set,
    train_proportion: float,
    val_proportion: float,
    test_proportion: float,
    uniref_db,
):
    # Convert `protein` to UniRef90
    protein = uniref_db.get(protein.encode("utf8"))

    if protein is None:
        return None, None
    else:
        protein = protein.decode("utf8")

    if protein in train_proteins:
        return "train", protein
    elif protein in val_proteins:
        return "val", protein
    elif protein in test_proteins:
        return "test", protein

    num_interactions = len(train_interactions) + len(val_interactions) + len(test_interactions)

    num_zero = 0

    if num_interactions > 0:
        current_train_prop = len(train_interactions) / num_interactions
        current_val_prop = len(val_interactions) / num_interactions
        current_test_prop = len(test_interactions) / num_interactions

        if current_test_prop > test_proportion:
            test_proportion = 0
            num_zero += 1

        if current_train_prop > train_proportion:
            train_proportion = 0
            num_zero += 1

        if current_val_prop > val_proportion:
            val_proportion = 0
            num_zero += 1

        if train_proportion == 0 and test_proportion == 0 and val_proportion == 0:
            train_proportion = 1 / 3
            val_proportion = 1 / 3
        else:
            remaining_proportion = 1 - (
                train_proportion + test_proportion + val_proportion
            )

            if train_proportion != 0:
                train_proportion += remaining_proportion / (3 - num_zero)

            if val_proportion != 0:
                val_proportion += remaining_proportion / (3 - num_zero)

            if test_proportion != 0:
                test_proportion += remaining_proportion / (3 - num_zero)

    split = ["train", "val", "test"][threeway_split(train_proportion, val_proportion)]

    # Return UniRef90 id

    if split == "train":
        return "train", protein
    elif split == "val":
        return "val", protein
    elif split == "test":
        return "test", protein


class MaxIterReached(Exception):
    pass


def get_random_negative(
    proteins: List,
    counts: List[int],
    train_interactions: List,
    val_interactions: List,
    test_interactions: List,
    train_neg_interactions: List,
    val_neg_interactions: List,
    test_neg_interactions: List,
    weighted_random: bool
 ) -> Tuple[str, str]:

    MAX_ITER = 1000

    iters = 0

    while True:
        if weighted_random:
            protein1, protein2 = sample(population=proteins, k=2, counts=counts)
        else:
            protein1, protein2 = sample(population=proteins, k=2)

        protein1, protein2 = canonical_id_sort(protein1, protein2)

        if (
            (protein1, protein2, True) not in train_interactions
            and (protein1, protein2, True) not in val_interactions
            and (protein1, protein2, True) not in test_interactions
            and (protein1, protein2, False) not in train_neg_interactions
            and (protein1, protein2, False) not in val_neg_interactions
            and (protein1, protein2, False) not in test_neg_interactions
        ):
            break

        iters += 1

        if iters > MAX_ITER:
            raise MaxIterReached

    return protein1, protein2


def parse_common(common_path: Path) -> Generator[Tuple[str, str, str], None, None]:
    interaction_num_rows = 0

    with gzip.open(common_path, "rt") as f:
        for _ in f:
            interaction_num_rows += 1

    with gzip.open(common_path, "rt") as f:
        reader = csv.DictReader(f)

        for row in track(
                reader,
                total=interaction_num_rows,
                description=f"Loading edges",
        ):

            if row["protein1"] == "n.a" or row["protein2"] == "n.a":
                continue

            yield row["protein1"], row["protein2"], row["score"]


def c_ify(
    console,
    common_path: Path,
    processed_folder: Path,
    c_types: List[int],
    train_proportion: float,
    val_proportion: float,
    test_proportion: float,
    neg_proportion: float,
    uniref_threshold: int,
    score_key: Optional[str] = None,
    score_threshold: Optional[float] = None,
    preloaded_protein_splits: Optional[dict] = None,
    seed: int = 8675309,
    negative_pools: Optional[List[dict]] = None,
    trim_unseen_proteins: bool = False,
    negatives_path: Optional[Path] = None,
    taxon: Optional[int] = None,
    weighted_random: bool = False,
    exclude_preloaded_from_neg: bool = True
) -> Tuple[Dict[str, Dict[str, Tuple[str, str, bool]]], Dict[str, List[str]]]:

    set_seed(seed)

    uniref_db_path = str(processed_folder / f"uniref{uniref_threshold}_members_upkb.leveldb")

    if not os.path.isdir(uniref_db_path):
        raise IOError(
            f"Can't read the UniRef90/UPKB mapping database at {uniref_db_path}. Make sure to run the 'process uniref' task."
        )

    uniref_db = plyvel.DB(uniref_db_path)

    # START validating arguments
    for c_type in c_types:
        if c_type not in [1, 2, 3]:
            raise ValueError(f"c_type must be one of 1, 2, or 3, not {c_type}")

    if preloaded_protein_splits is not None and (
        "train_proteins" not in preloaded_protein_splits.keys()
        or "val_proteins" not in preloaded_protein_splits.keys()
        or "test_proteins" not in preloaded_protein_splits.keys()
    ):
        raise ValueError(
            '"preloaded_protein_splits" must have the keys "train_proteins", "val_proteins", and'
            '"test_proteins."'
        )
    # END validating arguments

    # Seed the split proteins with preloaded_protein_splits
    if preloaded_protein_splits is not None:
        train_proteins_upkb = set(preloaded_protein_splits["train_proteins"])
        val_proteins_upkb = set(preloaded_protein_splits["val_proteins"])
        test_proteins_upkb = set(preloaded_protein_splits["test_proteins"])

        train_proteins_uniref = set()
        val_proteins_uniref = set()
        test_proteins_uniref = set()

        for protein_upkb in train_proteins_upkb:
            protein_uniref = uniref_db.get(str(protein_upkb).encode("utf8"))

            if protein_uniref is None:
                console.log(
                    f"WARNING: Preloaded protein from train split with id {protein_upkb} has missing UniRef id. Omitting."
                )
            else:

                protein_uniref = protein_uniref.decode("utf8")

                # make sure no overlap slips through
                in_val = protein_uniref in val_proteins_uniref
                in_test = protein_uniref in test_proteins_uniref

                if in_val:
                    train_proteins_uniref = train_proteins_uniref - set(protein_uniref)
                    train_proteins_upkb = train_proteins_upkb - set(protein_upkb)

                if in_test:
                    train_proteins_uniref = train_proteins_uniref - set(protein_uniref)
                    train_proteins_upkb = train_proteins_upkb - set(protein_upkb)

                if not in_val and not in_test:
                    train_proteins_uniref.add(protein_uniref)

        for protein_upkb in val_proteins_upkb:
            protein_uniref = uniref_db.get(str(protein_upkb).encode("utf8"))

            if protein_uniref is None:
                console.log(
                    f"WARNING: Preloaded protein from val split with id {protein_upkb} has missing UniRef id. Omitting."
                )
            else:
                protein_uniref = protein_uniref.decode("utf8")

                # make sure no overlap slips through
                in_train = protein_uniref in train_proteins_uniref
                in_test = protein_uniref in test_proteins_uniref

                if in_train:
                    val_proteins_uniref = val_proteins_uniref - set(protein_uniref)
                    val_proteins_upkb = val_proteins_upkb - set(protein_upkb)

                if in_test:
                    val_proteins_uniref = val_proteins_uniref - set(protein_uniref)
                    val_proteins_upkb = val_proteins_upkb - set(protein_upkb)

                if not in_train and not in_test:
                    val_proteins_uniref.add(protein_uniref)

        for protein_upkb in test_proteins_upkb:
            protein_uniref = uniref_db.get(str(protein_upkb).encode("utf8"))

            if protein_uniref is None:
                console.log(
                    f"WARNING: Preloaded protein from test split with id {protein_upkb} has missing UniRef id. Omitting."
                )
            else:
                protein_uniref = protein_uniref.decode("utf8")

                # make sure no overlap slips through
                in_train = protein_uniref in train_proteins_uniref
                in_val = protein_uniref in val_proteins_uniref

                if in_train:
                    test_proteins_uniref = test_proteins_uniref - set(protein_uniref)
                    test_proteins_upkb = test_proteins_upkb - set(protein_upkb)

                if in_val:
                    test_proteins_uniref = test_proteins_uniref - set(protein_uniref)
                    test_proteins_upkb = test_proteins_upkb - set(protein_upkb)

                if not in_val and not in_train:
                    test_proteins_uniref.add(protein_uniref)

        train_val_overlap = train_proteins_uniref.intersection(val_proteins_uniref)
        train_test_overlap = train_proteins_uniref.intersection(test_proteins_uniref)
        test_val_overlap = test_proteins_uniref.intersection(val_proteins_uniref)

        if len(train_val_overlap) > 0:
            console.log(
                f"WARNING: {len(train_val_overlap)} proteins are common in training and val splits!"
            )

        if len(train_test_overlap) > 0:
            console.log(
                f"WARNING: {len(train_test_overlap)} proteins are common in training and test splits!"
            )

        if len(test_val_overlap) > 0:
            console.log(
                f"WARNING: {len(test_val_overlap)} proteins are common in test and val splits!"
            )

        console.log(
            f"Preloaded {len(train_proteins_uniref)} train proteins, {len(val_proteins_uniref)} val proteins, {len(test_proteins_uniref)} test proteins."
        )

    else:
        train_proteins_uniref = set()
        val_proteins_uniref = set()
        test_proteins_uniref = set()

        train_proteins_upkb = set()
        val_proteins_upkb = set()
        test_proteins_upkb = set()

    # Set-up interaction sets
    train_pos_interactions = {c_type: set() for c_type in c_types}
    val_pos_interactions = {c_type: set() for c_type in c_types}
    test_pos_interactions = {c_type: set() for c_type in c_types}

    train_neg_interactions = {c_type: set() for c_type in c_types}
    val_neg_interactions = {c_type: set() for c_type in c_types}
    test_neg_interactions = {c_type: set() for c_type in c_types}

    # protein counts are important for preserving degree when sampling negatives
    protein_counts = defaultdict(lambda: 0)

    # Load edges
    if negatives_path is None:
        paths = [("positive", common_path)]
    else:
        paths = [("positive", common_path), ("negative", negatives_path)]

    if taxon is not None:
        db_taxon = plyvel.DB(
            str(processed_folder / f"uniprot_sequences_{taxon}.leveldb"),
            create_if_missing=True,
        )
        out_of_taxon = 0

    interaction_num_rows = 0

    for edge_type, path in paths:
        console.log(f"Making edges from {path}")

        interaction_num_rows += 1

        for protein1_upkb, protein2_upkb, score in parse_common(common_path):

            protein_counts[protein1_upkb] += 1
            protein_counts[protein2_upkb] += 1

            # Filter edges on their score, if score_key specified
            if score_key:
                row_score = float(mitab_deserialize(score)[score_key])

                if row_score < score_threshold:
                    continue

            if protein1_upkb == "n.a" or protein2_upkb == "n.a":
                continue

            if taxon is not None and (
                    db_taxon.get(protein1_upkb.encode("utf8")) is None
                    or db_taxon.get(protein2_upkb.encode("utf8")) is None
            ):
                out_of_taxon += 1
                continue

            if edge_type == "positive":
                if 3 in c_types:
                    train_interactions = train_pos_interactions[3]
                elif 2 in c_types:
                    train_interactions = train_pos_interactions[2]
                else:
                    train_interactions = train_pos_interactions[1]

                if 3 in c_types:
                    val_interactions = val_pos_interactions[3]
                elif 2 in c_types:
                    val_interactions = val_pos_interactions[2]
                else:
                    val_interactions = val_pos_interactions[1]

                if 3 in c_types:
                    test_interactions = test_pos_interactions[3]
                elif 2 in c_types:
                    test_interactions = test_pos_interactions[2]
                else:
                    test_interactions = test_pos_interactions[1]
            else:
                if 3 in c_types:
                    train_interactions = train_neg_interactions[3]
                elif 2 in c_types:
                    train_interactions = train_neg_interactions[2]
                else:
                    train_interactions = train_neg_interactions[1]

                if 3 in c_types:
                    val_interactions = val_neg_interactions[3]
                elif 2 in c_types:
                    val_interactions = val_neg_interactions[2]
                else:
                    val_interactions = val_neg_interactions[1]

                if 3 in c_types:
                    test_interactions = test_neg_interactions[3]
                elif 2 in c_types:
                    test_interactions = test_neg_interactions[2]
                else:
                    test_interactions = test_neg_interactions[1]

            protein1_split, protein1_uniref = get_split(
                protein1_upkb,
                train_proteins_uniref,
                val_proteins_uniref,
                test_proteins_uniref,
                train_interactions,
                val_interactions,
                test_interactions,
                train_proportion,
                val_proportion,
                test_proportion,
                uniref_db
            )

            protein2_split, protein2_uniref = get_split(
                protein2_upkb,
                train_proteins_uniref,
                val_proteins_uniref,
                test_proteins_uniref,
                train_interactions,
                val_interactions,
                test_interactions,
                train_proportion,
                val_proportion,
                test_proportion,
                uniref_db
            )


            if protein1_split == "train":
                train_proteins_uniref.add(protein1_uniref)
                train_proteins_upkb.add(protein1_upkb)
            elif protein1_split == "val":
                val_proteins_uniref.add(protein1_uniref)
                val_proteins_upkb.add(protein1_upkb)
            elif protein1_split == "test":
                test_proteins_uniref.add(protein1_uniref)
                test_proteins_upkb.add(protein1_upkb)

            if protein2_split == "train":
                train_proteins_uniref.add(protein2_uniref)
                train_proteins_upkb.add(protein2_upkb)
            elif protein2_split == "val":
                val_proteins_uniref.add(protein2_uniref)
                val_proteins_upkb.add(protein2_upkb)
            elif protein2_split == "test":
                test_proteins_uniref.add(protein2_uniref)
                test_proteins_upkb.add(protein2_upkb)

            protein1_upkb, protein2_upkb = canonical_id_sort(
                protein1_upkb, protein2_upkb
            )

            if 1 in c_types:
                dice = random()

                if dice < train_proportion:

                    if edge_type == "positive":
                        train_pos_interactions[1].add((protein1_upkb, protein2_upkb, True))
                    else:
                        train_neg_interactions[1].add((protein1_upkb, protein2_upkb, False))

                elif train_proportion <= dice < train_proportion + val_proportion:
                    if edge_type == "positive":
                        val_pos_interactions[1].add((protein1_upkb, protein2_upkb, True))
                    else:
                        val_neg_interactions[1].add((protein1_upkb, protein2_upkb, False))
                else:
                    if edge_type == "positive":
                        test_pos_interactions[1].add((protein1_upkb, protein2_upkb, True))
                    else:
                        test_neg_interactions[1].add((protein1_upkb, protein2_upkb, False))

            if 2 in c_types:

                if {protein1_split, protein2_split} == {"train", "train"}:
                    if edge_type == "positive":
                        train_pos_interactions[2].add((protein1_upkb, protein2_upkb, True))
                    else:
                        train_neg_interactions[2].add((protein1_upkb, protein2_upkb, False))

                elif {protein1_split, protein2_split} == {"train", "test"} or {
                    protein1_split,
                    protein2_split,
                } == {"test", "test"}:

                    if edge_type == "positive":
                        test_pos_interactions[2].add((protein1_upkb, protein2_upkb, True))
                    else:
                        test_neg_interactions[2].add((protein1_upkb, protein2_upkb, False))

                elif {protein1_split, protein2_split} == {"train", "val"} or {
                    protein1_split,
                    protein2_split,
                } == {"val", "val"}:

                    if edge_type == "positive":
                        val_pos_interactions[2].add((protein1_upkb, protein2_upkb, True))
                    else:
                        val_neg_interactions[2].add((protein1_upkb, protein2_upkb, False))

            if 3 in c_types:

                if {protein1_split, protein2_split} == {"train", "train"}:

                    if edge_type == "positive":
                        train_pos_interactions[3].add((protein1_upkb, protein2_upkb, True))
                    else:
                        train_neg_interactions[3].add((protein1_upkb, protein2_upkb, False))

                elif {protein1_split, protein2_split} == {"test", "test"}:

                    if edge_type == "positive":
                        test_pos_interactions[3].add((protein1_upkb, protein2_upkb, True))
                    else:
                        test_neg_interactions[3].add((protein1_upkb, protein2_upkb, False))

                elif {protein1_split, protein2_split} == {"val", "val"}:

                    if edge_type == "positive":
                        val_pos_interactions[3].add((protein1_upkb, protein2_upkb, True))
                    else:
                        val_neg_interactions[3].add((protein1_upkb, protein2_upkb, False))

    console.log(f"# train proteins: {len(train_proteins_upkb)}")
    console.log(f"# test proteins: {len(test_proteins_upkb)}")
    console.log(f"# val proteins: {len(val_proteins_upkb)}")

    if taxon is not None:
        console.log(
            f"{out_of_taxon} ({out_of_taxon * 100 / interaction_num_rows:.3})% edges are out-of-taxon. Skipping them."
        )

    for c in [1, 2, 3]:
        console.log(f"+ C{c}")
        console.log(f"+--- # train pos interactions: {len(train_pos_interactions[c])}")
        console.log(f"+--- # val pos interactions: {len(val_pos_interactions[c])}")
        console.log(f"+--- # test pos interactions: {len(test_pos_interactions[c])}")

    if negatives_path is None:
        # making negative edges

        # so what's up with the exclusion sets?
        # they exist to remove proteins that belong to the preloaded_protein_split from the population of negative proteins.
        # this is useful when you include preloaded_protein_splits from other species.
        # if however your preloaded_protein_splits are from the same species, then this can eliminate too many (possibly all)
        # proteins from the pool of proteins from which negatives are drawn. 
        if preloaded_protein_splits is None or exclude_preloaded_from_neg is False:
            train_exclusion_set = set()
            val_exclusion_set = set()
            test_exclusion_set = set()
        else:
            train_exclusion_set = set(preloaded_protein_splits["train_proteins"])
            val_exclusion_set = set(preloaded_protein_splits["val_proteins"])
            test_exclusion_set = set(preloaded_protein_splits["test_proteins"])

        all_proteins = list(
            (train_proteins_upkb).union(val_proteins_upkb).union(test_proteins_upkb)
        )
        train_proteins_upkb = list(train_proteins_upkb)
        train_protein_upkb_counts = [
            protein_counts[protein] for protein in train_proteins_upkb
        ]

        val_proteins_upkb = list(val_proteins_upkb)
        val_protein_upkb_counts = [
            protein_counts[protein] for protein in val_proteins_upkb
        ]

        test_proteins_upkb = list(test_proteins_upkb)
        test_protein_upkb_counts = [
            protein_counts[protein] for protein in test_proteins_upkb
        ]

        console.log("Generating negatives...")

        if negative_pools is not None:
            console.log("Found negative pools.")
            negative_pool_size = sum(
                [
                    len(p["train"]) + len(p["test"]) + len(p["val"])
                    for p in negative_pools
                ]
            )
            negative_ratios = [
                (len(p["train"]) + len(p["test"]) + len(p["val"])) / negative_pool_size
                for p in negative_pools
            ]
            console.log("The pool relative sizes are: ", negative_ratios)

            for idx, negative_pool in enumerate(negative_pools):
                console.log(f"Negative pool {idx} stats")
                console.log(f'\t# Train: {len(negative_pool["train"])}')
                console.log(f'\t# Val: {len(negative_pool["val"])}')
                console.log(f'\t# Test: {len(negative_pool["test"])}')

        all_intxns = [
            ("train", train_pos_interactions),
            ("val", val_pos_interactions),
            ("test", test_pos_interactions),
        ]
        for intxns_split, intxns in all_intxns:
            for c_type in c_types:
                for _ in track(
                    range(int(len(intxns[c_type]) * neg_proportion)),
                    description=f"{intxns_split} C{c_type}",
                ):
                    if c_type == 1:
                        proteins = set(all_proteins) - (
                            train_exclusion_set.union(val_exclusion_set).union(
                                test_exclusion_set
                            )
                        )
                        counts = [protein_counts[protein] for protein in all_proteins]
                    else:
                        if negative_pools is not None:
                            negative_pool = choices(negative_pools, negative_ratios)[0]
                            if intxns_split == "train":
                                proteins = list(negative_pool["train"])
                                counts = [
                                    protein_counts[protein]
                                    for protein in negative_pool["train"]
                                ]
                            elif intxns_split == "val":
                                proteins = list(negative_pool["val"])
                                counts = [
                                    protein_counts[protein]
                                    for protein in negative_pool["val"]
                                ]
                            else:
                                proteins = list(negative_pool["test"])
                                counts = [
                                    protein_counts[protein]
                                    for protein in negative_pool["test"]
                                ]
                        else:
                            if intxns_split == "train":
                                proteins = set(train_proteins_upkb) - train_exclusion_set
                                counts = train_protein_upkb_counts
                            elif intxns_split == "val":
                                proteins = set(val_proteins_upkb) - val_exclusion_set
                                counts = val_protein_upkb_counts
                            else:
                                proteins = set(test_proteins_upkb) - test_exclusion_set
                                counts = test_protein_upkb_counts

                    if len(proteins) == 0:
                        console.log(f"[b][red]ERROR![/red][b] No proteins are valid candidates to draw negative samples from. This can happen when you preload protein splits.")

                    protein1, protein2 = get_random_negative(
                        proteins,
                        counts,
                        train_pos_interactions[c_type],
                        val_pos_interactions[c_type],
                        test_pos_interactions[c_type],
                        train_neg_interactions[c_type],
                        val_neg_interactions[c_type],
                        test_neg_interactions[c_type],
                        weighted_random,
                    )

                    neg_intxn = (protein1, protein2, False)
                    if intxns_split == "train":
                        train_neg_interactions[c_type].add(neg_intxn)
                    elif intxns_split == "val":
                        val_neg_interactions[c_type].add(neg_intxn)
                    elif intxns_split == "test":
                        test_neg_interactions[c_type].add(neg_intxn)

    interactions = {f"C{c}": None for c in c_types}

    for c_type in track(c_types, description="Saving files..."):
        train_interactions = list(
            train_pos_interactions[c_type].union(train_neg_interactions[c_type])
        )
        val_interactions = list(
            val_pos_interactions[c_type].union(val_neg_interactions[c_type])
        )
        test_interactions = list(
            test_pos_interactions[c_type].union(test_neg_interactions[c_type])
        )

        shuffle(train_interactions)
        shuffle(val_interactions)
        shuffle(test_interactions)

        interactions[f"C{c_type}"] = {
            "train": train_interactions,
            "val": val_interactions,
            "test": test_interactions,
        }

    proteins = {
        "train": list(train_proteins_upkb),
        "test": list(test_proteins_upkb),
        "val": list(val_proteins_upkb),
    }

    return interactions, proteins


def get_seq(upkb_ac: str, db, session):
    seq = db.get(upkb_ac.encode("utf8"))

    if seq is None:
        print(f"Missing seq for {upkb_ac}, retrieving from unirpot...")
        time.sleep(1)
        response = session.get(f"https://rest.uniprot.org/uniprotkb/{upkb_ac}.fasta")
        response.raise_for_status()
        seq = "".join(response.text.splitlines()[1:])
        db.put(upkb_ac.encode("utf8"), seq.encode("utf8"))
    else:
        seq = seq.decode("utf8")

    return seq


def upkb_query(query_url: str) -> Generator[List[str], None, None]:
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    def get_next_link(headers):
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def get_batch(batch_url):
        while batch_url:
            response = session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = get_next_link(response.headers)

    args = urlparse(query_url)
    queries = parse_qs(args.query)

    # if format not specified, specify it as tsv
    if 'format' not in queries:
        query_url += "&format=tsv"
    # if format specified as anything but tsv, raise error
    elif 'tsv' not in queries['format']:
        raise ValueError("Query URL must be TSV")

    for batch, total in get_batch(query_url):
        for line in batch.text.splitlines()[1:]:
            yield line.split('\t')


def verify_taxon(upkb_data: List[str], taxon_ids: Set[int]) -> Tuple[bool, Set[int], Dict[str, int]]:

    observed_taxon_ids = set()

    upkb_taxa = {}

    # upkb_ac, taxon_id, seq
    for upkb_ac, taxon_id, _ in upkb_data:
        observed_taxon_ids.add(int(taxon_id))
        upkb_taxa[upkb_ac] = taxon_id

    if observed_taxon_ids != taxon_ids:
        return False, observed_taxon_ids, upkb_taxa
    else:
        return True, observed_taxon_ids, upkb_taxa


def verify_observed_splits(observed_proteins_split: Dict[int, Dict[str, Set[str]]]) -> Tuple[bool, Dict[int, Dict[str, Optional[Set]]]]:

    observed_overlap = {c: {split_pair: None for split_pair in ["train/val", "train/test", "test/val"]} for c in [1, 2, 3]}

    for c in [1, 2, 3]:
        observed_overlap[c]['train/val'] = observed_proteins_split[c]["train"].intersection(observed_proteins_split[c]["val"])
        observed_overlap[c]['train/test'] = observed_proteins_split[c]["train"].intersection(observed_proteins_split[c]["test"])
        observed_overlap[c]['test/val'] = observed_proteins_split[c]["test"].intersection(observed_proteins_split[c]["val"])

    if len(observed_overlap[3]['train/val']) == 0 and len(observed_overlap[3]['train/test']) == 0 and len(observed_overlap[3]['test/val']) == 0:
        verify_pass = True
    else:
        verify_pass = False

    return verify_pass, observed_overlap
