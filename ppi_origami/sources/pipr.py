import json
from collections import defaultdict

import tables as tb


def rapppid_to_pipr(
    dataset_path, c_type, seq_path, train_path, test_path, val_path, metadata_path
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
            f.write(f"{row['name'].decode('utf8')}\t{row['sequence'].decode('utf8')}\n")

    with open(train_path, "w") as f:
        if c_type == 3:
            c_table = rapppid_dataset.root.interactions.c3.c3_train
        elif c_type == 2:
            c_table = rapppid_dataset.root.interactions.c2.c2_train
        elif c_type == 1:
            c_table = rapppid_dataset.root.interactions.c1.c1_train
        else:
            raise ValueError("Unexpected value for c_type")

        f.write("v1\tv2\tlabel\n")

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

        f.write("v1\tv2\tlabel\n")

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

        f.write("v1\tv2\tlabel\n")

        for row in c_table:
            if row["label"] is True:
                label = 1
            else:
                label = 0
            f.write(
                f"{row['protein_id1'].decode('utf8')}\t{row['protein_id2'].decode('utf8')}\t{label}\n"
            )
