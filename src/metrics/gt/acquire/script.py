import os
import sys
import io
import zipfile
import tempfile
import requests
from pathlib import Path


sys.path.append("src/utils")
from util import download_and_uncompress_zip


if __name__ == "__main__":

    download_and_uncompress_zip(
        "https://zenodo.org/records/3701939/files/BEELINE-Networks.zip?download=1",
        "resources/grn_benchmark/ground_truth/beeline"
    )
