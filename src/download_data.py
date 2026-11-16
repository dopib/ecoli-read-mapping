"""Download the two input files. Skip reads if a local copy already exists."""

import argparse
import gzip
from pathlib import Path
import shutil
import urllib.request


REFERENCE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/"
    "GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
)
READS_URL = "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022075/ERR022075_1.fastq.gz"


def download(url, destination):
    if destination.exists():
        print(f"Already exists: {destination}")
        return
    temporary = destination.with_suffix(destination.suffix + ".part")
    print(f"Downloading {url}", flush=True)
    with urllib.request.urlopen(url, timeout=120) as source, temporary.open("wb") as output:
        shutil.copyfileobj(source, output, length=1024 * 1024)
    temporary.replace(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path("data"))
    parser.add_argument("--reference-only", action="store_true")
    args = parser.parse_args()
    args.directory.mkdir(parents=True, exist_ok=True)
    compressed = args.directory / "reference.fna.gz"
    download(REFERENCE_URL, compressed)
    reference = args.directory / "reference.fna"
    if not reference.exists():
        with gzip.open(compressed, "rb") as source, reference.open("wb") as output:
            shutil.copyfileobj(source, output)
    if not args.reference_only:
        download(READS_URL, args.directory / "ERR022075_1.fastq.gz")


if __name__ == "__main__":
    main()
