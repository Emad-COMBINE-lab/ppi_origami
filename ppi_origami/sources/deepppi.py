import json
from collections import defaultdict

import tables as tb


def rapppid_to_deepppi(
    dataset_path, c_type, train_path, test_path, val_path, metadata_path
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

    seqs = dict()
    for row in rapppid_dataset.root.sequences:
        seqs[row["name"].decode("utf8")] = row["sequence"].decode("utf8")

    with open(train_path, "w") as f:
        if c_type == 3:
            c_table = rapppid_dataset.root.interactions.c3.c3_train
        elif c_type == 2:
            c_table = rapppid_dataset.root.interactions.c2.c2_train
        elif c_type == 1:
            c_table = rapppid_dataset.root.interactions.c1.c1_train
        else:
            raise ValueError("Unexpected value for c_type")

        f.write("1166\n")

        for row in c_table:
            if row["label"] is True:
                label = 1
            else:
                label = 0

            pid1 = row["protein_id1"].decode("utf8")
            pid2 = row["protein_id2"].decode("utf8")
            seq1 = seqs[pid1]
            seq2 = seqs[pid2]

            f.write(f"{pid1} {pid2} {seq1} {seq2} {label}\n")

    with open(test_path, "w") as f:
        if c_type == 3:
            c_table = rapppid_dataset.root.interactions.c3.c3_test
        elif c_type == 2:
            c_table = rapppid_dataset.root.interactions.c2.c2_test
        elif c_type == 1:
            c_table = rapppid_dataset.root.interactions.c1.c1_test
        else:
            raise ValueError("Unexpected value for c_type")

        f.write("1166\n")

        for row in c_table:
            if row["label"] is True:
                label = 1
            else:
                label = 0

            pid1 = row["protein_id1"].decode("utf8")
            pid2 = row["protein_id2"].decode("utf8")
            seq1 = seqs[pid1]
            seq2 = seqs[pid2]

            f.write(f"{pid1} {pid2} {seq1} {seq2} {label}\n")

    with open(val_path, "w") as f:
        if c_type == 3:
            c_table = rapppid_dataset.root.interactions.c3.c3_val
        elif c_type == 2:
            c_table = rapppid_dataset.root.interactions.c2.c2_val
        elif c_type == 1:
            c_table = rapppid_dataset.root.interactions.c1.c1_val
        else:
            raise ValueError("Unexpected value for c_type")

        f.write("1166\n")

        for row in c_table:
            if row["label"] is True:
                label = 1
            else:
                label = 0

            pid1 = row["protein_id1"].decode("utf8")
            pid2 = row["protein_id2"].decode("utf8")
            seq1 = seqs[pid1]
            seq2 = seqs[pid2]

            f.write(f"{pid1} {pid2} {seq1} {seq2} {label}\n")
