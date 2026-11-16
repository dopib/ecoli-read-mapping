"""Check the report inputs against Bowtie's log and coverage tables."""

import csv
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


def main():
    directory = ROOT / "results"
    m = json.loads((directory / "metrics.json").read_text())
    log = (directory / "bowtie2.log").read_text()
    expected = {
        "total_reads": int(re.search(r"^(\d+) reads;", log).group(1)),
        "unmapped_reads": int(re.search(r"(\d+) \([^\n]+\) aligned 0 times", log).group(1)),
        "one_placement_reads": int(re.search(r"(\d+) \([^\n]+\) aligned exactly 1 time", log).group(1)),
        "multiple_placement_reads": int(re.search(r"(\d+) \([^\n]+\) aligned >1 times", log).group(1)),
    }
    for name, count in expected.items():
        assert m[name] == count, (name, m[name], count)
    assert m["total_reads"] == m["fastq_validated_reads"]
    assert m["input_bases"] == m["fastq_validated_bases"]
    assert m["mapped_reads"] + m["unmapped_reads"] == m["total_reads"]
    assert m["one_placement_reads"] + m["multiple_placement_reads"] == m["mapped_reads"]
    histogram = {}
    with (directory / "depth_histogram.tsv").open() as source:
        for row in csv.DictReader(source, delimiter="\t"):
            histogram[int(row["depth"])] = int(row["genome_positions"])
    assert sum(histogram.values()) == m["genome_length"]
    c = m["coverage"]
    assert sum(depth * count for depth, count in histogram.items()) == c["depth_sum"]
    assert sum(count for depth, count in histogram.items() if depth >= 1) == c["covered_positions"]
    for threshold in (10, 30):
        assert sum(count for depth, count in histogram.items() if depth >= threshold) == c[f"positions_at_least_{threshold}x"]
    bases, weighted_depth = 0, 0.0
    with (directory / "coverage_10kb.tsv").open() as source:
        for row in csv.DictReader(source, delimiter="\t"):
            length = int(row["end_exclusive"]) - int(row["start_0based"])
            bases += length
            weighted_depth += length * float(row["mean_depth"])
    assert bases == m["genome_length"]
    assert abs(weighted_depth - c["depth_sum"]) <= bases * 0.000001
    quality = m["quality"]
    edits = {int(k): v for k, v in quality["edit_distance_histogram"].items()}
    assert sum(edits.values()) == m["mapped_reads"]
    assert sum(n * count for n, count in edits.items()) == quality["edit_distance_sum"]
    assert quality["aligned_read_bases"] + quality["soft_clipped_bases"] == quality["mapped_read_bases"]
    print("Checks passed: FASTQ counts, Bowtie 2 log, depth histogram, coverage windows and edit distances.")


if __name__ == "__main__":
    main()
