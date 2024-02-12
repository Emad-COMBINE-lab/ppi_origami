import csv
import gzip
from pathlib import Path
from os.path import isfile
from typing import Optional

from .dscript import parse_dscript
from .string import parse_string


def to_common(processed_folder: Path, source: str, taxon: Optional[int], identifier: str, version: Optional[str] = None) -> None:
    """
    Convert PPI datasets from multiple sources into a common format that is used by
    this library.

    :param version: The version of the source database to process.
    :param processed_folder: The folder wherein we keep processed raw files.
    :param source: The source dataset. Accepted values currently are "string_links_detailed" and "dscript".
    :param taxon: The NCBI taxon number of the dataset to convert.
    :param identifier: Which protein identifier to use. Must be one of 'upkb' or 'uniref50', 'uniref90', or 'uniref100'.
    :return: None
    """

    valid_identifiers = ["upkb", "uniref50", "uniref90", "uniref100"]
    valid_identifiers_str = ", ".join(valid_identifiers)

    if identifier not in valid_identifiers:
        raise ValueError(f'Identifier must be one of {valid_identifiers_str}')

    valid_sources = ["string_links_detailed", "dscript"]
    valid_sources_str = ", ".join(valid_sources)

    if source not in valid_sources:
        raise ValueError(f'Source must be one of {valid_sources_str}')

    if source == "string_links_detailed":

        if taxon:
            source_path = (
                    processed_folder
                    / f"string_{taxon}.protein.links.detailed.v{version}_{identifier}.csv.gz"
            )
        else:
            source_path = (
                    processed_folder
                    / f"string.protein.links.detailed.v{version}_{identifier}.csv.gz"
            )

        if not isfile(source_path):
            raise IOError(f"Couldn't find a file for STRING v{version} dataset at {source_path}")

        parser = parse_string(source_path, identifier)

        out_path = (
                processed_folder
                / f"common_string_{taxon}.protein.links.detailed.v{version}_{identifier}.csv.gz"
        )
    elif source == "dscript":

        if taxon is None:
            raise ValueError("taxon must not be None if source is dscript.")

        source_path = (
            processed_folder / f"dscript_{taxon}_train_{identifier}.csv.gz",
            processed_folder / f"dscript_{taxon}_test_{identifier}.csv.gz",
        )

        if not isfile(source_path[0]) or not isfile(source_path[1]):
            raise IOError(f"Couldn't find a file for the D-SCRIPT dataset.")

        parser = parse_dscript(source_path, identifier, 1)

        out_path = processed_folder / f"common_dscript_{taxon}_{identifier}.csv.gz"

    with gzip.open(out_path, "wt") as f_out:
        writer = csv.DictWriter(
            f_out, fieldnames=["interaction_id", "protein1", "protein2", "score"]
        )
        writer.writeheader()

        for row in parser:
            writer.writerow(row)
