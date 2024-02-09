import csv
import time
import urllib.parse
from pathlib import Path
from typing import Dict, Iterable, Set, List, Any

import tables as tb
from ppi_origami.sources.utils import canonical_interaction_id, upkb_query
from rich.progress import track


def dscript_rapppid_equivalence(
    dscript_path: Path,
    rapppid_path: Path,
    c_types: Iterable[int] = (1, 2, 3),
) -> Dict[str, Dict[str, set]]:
    """
    Test the equivalence of a D-SCRIPT dataset and a RAPPPID dataset.

    :param dscript_path: The path to the D-SCRIPT dataset
    :param rapppid_path: The path to the RAPPPID dataset
    :param c_types: The C-types to check equivalence for.
    :return: A dictionary where keys are the C-types, and the value are a set of mismatched interactions
    """

    interactions = dict()

    for c_type in c_types:
        interactions[f"C{c_type}"] = {"dscript": set(), "rapppid": set()}

        if c_type not in [1, 2, 3]:
            raise ValueError(
                "Invalid value for c_types, integers must be between 1 and 3 inclusively."
            )

        # Get D-SCRIPT interactions
        paths = [
            dscript_path / f"C{c_type}/{split}.tsv"
            for split in ["train", "test", "val"]
        ]

        for path in paths:
            with open(path) as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    p1 = row[0]
                    p2 = row[1]
                    label = int(row[2])

                    pair = canonical_interaction_id(p1, p2)

                    interactions[f"C{c_type}"]["dscript"].add((pair, label))

        # Get RAPPPID interactions
        with tb.open_file(str(rapppid_path)) as dataset:
            for split in ["train", "test", "val"]:
                for row in dataset.root["interactions"][f"c{c_type}"][
                    f"c{c_type}_{split}"
                ].iterrows():
                    p1 = row["protein_id1"].decode("utf8")
                    p2 = row["protein_id2"].decode("utf8")
                    label = 1 if row["label"] else 0

                    pair = canonical_interaction_id(p1, p2)

                    interactions[f"C{c_type}"]["rapppid"].add((pair, label))

        interactions[f"C{c_type}"] = {
            "rapppid": interactions[f"C{c_type}"]["rapppid"].difference(
                interactions[f"C{c_type}"]["dscript"]
            ),
            "dscript": interactions[f"C{c_type}"]["dscript"].difference(
                interactions[f"C{c_type}"]["rapppid"]
            ),
        }

    return interactions


def rapppid_stats(console, path: Path, c_types: List[int]):
    dataset = tb.open_file(str(path))

    tables = []

    for c in c_types:
        if c == 3:
            tables += [(3, "test", dataset.root.interactions.c3.c3_test),
                (3, "val", dataset.root.interactions.c3.c3_val),
                (3, "train", dataset.root.interactions.c3.c3_train)]
        elif c == 2:
            tables += [(2, "test", dataset.root.interactions.c2.c2_test),
                        (2, "val", dataset.root.interactions.c2.c2_val),
                        (2, "train", dataset.root.interactions.c2.c2_train)]
        elif c == 1:
            tables += [(1, "test", dataset.root.interactions.c1.c1_test),
                        (1, "val", dataset.root.interactions.c1.c1_val),
                        (1, "train", dataset.root.interactions.c1.c1_train)]

    interactions = {
        c: {split: set() for split in ["test", "val", "train"]} for c in c_types
    }
    proteins = {
        c: {split: { label: set() for label in ['pos', 'neg']} for split in ["test", "val", "train"]} for c in c_types
    }

    negative_interactions = {c: {split: 0 for split in ['test', 'val', 'train']} for c in c_types}
    positive_interactions = {c: {split: 0 for split in ['test', 'val', 'train']} for c in c_types}

    for c, split, table in tables:
        console.log(f"Computing stats for C{c} {split}...")
        for idx, row in enumerate(track(table.iterrows(), total=len(table))):
            interaction_id = canonical_interaction_id(
                row["protein_id1"].decode('utf8'), row["protein_id2"].decode('utf8')
            )

            interactions[c][split].add(interaction_id)

            if row['label']:
                label = 'pos'
                positive_interactions[c][split] += 1
            else:
                label = 'neg'
                negative_interactions[c][split] += 1

            proteins[c][split][label].add(row["protein_id1"])
            proteins[c][split][label].add(row["protein_id2"])

    total_interactions = []
    total_proteins = []
    neg_proteins = []
    pos_proteins = []

    for c in c_types:
        interaction_test_overlap = interactions[c]["train"].intersection(
            interactions[c]["test"]
        )
        interaction_val_overlap = interactions[c]["train"].intersection(
            interactions[c]["val"]
        )
        console.log(
            f"# of C{c} Test Interactions in Train:", len(interaction_test_overlap)
        )
        console.log(
            f"# of C{c} Val Interactions in Train:", len(interaction_val_overlap)
        )

        total_interactions = (
            total_interactions
            + list(interactions[c]["train"])
            + list(interactions[c]["test"])
            + list(interactions[c]["val"])
        )

        train_proteins = proteins[c]["train"]['pos'].union(proteins[c]["train"]['neg'])
        test_proteins = proteins[c]["test"]['pos'].union(proteins[c]["test"]['neg'])
        val_proteins = proteins[c]["val"]['pos'].union(proteins[c]["val"]['neg'])

        protein_test_interaction = train_proteins.intersection(test_proteins)
        protein_val_interaction = train_proteins.intersection(val_proteins)

        total_proteins = total_proteins + list(train_proteins) + list(test_proteins) + list(val_proteins)
        neg_proteins = neg_proteins + list(proteins[c]["train"]['neg']) + list(proteins[c]["test"]['neg']) + list(proteins[c]["val"]['neg'])
        pos_proteins = pos_proteins + list(proteins[c]["train"]['pos']) + list(proteins[c]["test"]['pos']) + list(proteins[c]["val"]['pos'])

        console.log(f"# of C{c} Test Proteins in Train", len(protein_test_interaction))
        console.log(f"# of C{c} Val Proteins in Train", len(protein_val_interaction))

    console.log("")
    console.log(f"# of interactions:", len(total_interactions))
    console.log(f"# of unique interactions:", len(set(total_interactions)))
    console.log("")
    console.log(f"# of proteins:", len(total_proteins))
    console.log(f"# of unique proteins:", len(set(total_proteins)))
    console.log("")
    console.log(f"# of proteins in positive interactions:", len(pos_proteins))
    console.log(f"# of unique proteins in positive interactions:", len(set(pos_proteins)))
    console.log("")
    console.log(f"# of proteins in negative interactions:", len(neg_proteins))
    console.log(f"# of unique proteins in negative interactions:", len(set(neg_proteins)))

    for c in c_types:
        for split in ['train', 'test', 'val']:
            npr = negative_interactions[c][split] / positive_interactions[c][split]
            console.log(f"C{3} {split} Negative/Positive Interactions: {npr}")


def get_uniprot_data(proteins: Set[str]) -> List[str]:

    def rq(batch: List[str], return_fields: str) -> Any:
        time.sleep(2)

        batch = urllib.parse.quote(" OR ".join(batch))

        for line in upkb_query(f'https://rest.uniprot.org/uniprotkb/search?&query={batch}&fields={return_fields}'):
            yield line

    batch = []
    return_fields = "accession,organism_id,sequence"

    results = []

    for idx, protein in enumerate(track(proteins)):

        batch.append(f"accession_id:{protein}")

        if (idx + 1) % 100 == 0:

            for line in rq(batch, return_fields):
                results.append(line)

            batch = []

    if len(batch) > 0:
        for line in rq(batch, return_fields):
            results.append(line)

    return results
