from pathlib import Path

from . import utils


def download(path: Path):
    url = "https://irefindex.vib.be/download/irefindex/data/archive/release_18.0/psi_mitab/MITAB2.6/All.mitab.06-11-2021.txt.zip"
    utils.download_file(
        url, path / "All.mitab.06-11-2021.txt.zip", "Downloading iRefWeb"
    )
