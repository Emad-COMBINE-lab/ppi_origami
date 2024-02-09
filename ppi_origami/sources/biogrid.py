from pathlib import Path
from shutil import unpack_archive

from . import utils


def download(path: Path, version: str = "4.4.224"):
    url = f"https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-{version}/BIOGRID-ORGANISM-{version}.mitab.zip"
    utils.download_file(
        url, path / f"BIOGRID-ORGANISM-{version}.mitab.zip", f"Downloading BioGRID v{version}"
    )
    unpack_archive(path / f"BIOGRID-ORGANISM-{version}.mitab.zip", path)
