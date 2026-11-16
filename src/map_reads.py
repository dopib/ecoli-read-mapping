"""Map a complete FASTQ file with Bowtie 2 and summarize its SAM stream.

Only one primary alignment per read contributes to quality and coverage.
Secondary alignments are used to detect reads with alternative placements.
"""

import argparse
from array import array
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import threading
import time

import numpy as np


CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def read_reference(path):
    sequences = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="ascii") as source:
        for line in source:
            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = 0
            else:
                sequences[name] += len(line.strip())
    return sequences


def feed_reads(path, destination, errors, counts):
    """Validate four-line FASTQ records while feeding batches to the aligner."""
    opener = gzip.open if str(path).endswith(".gz") else open
    batch = bytearray()
    try:
        with opener(path, "rb") as source:
            while True:
                header = source.readline()
                if not header:
                    break
                sequence = source.readline()
                separator = source.readline()
                quality = source.readline()
                if (not header.startswith(b"@") or not separator.startswith(b"+")
                        or not quality or len(sequence.rstrip()) != len(quality.rstrip())):
                    raise ValueError(f"Invalid FASTQ record {counts['reads'] + 1}")
                counts["reads"] += 1
                counts["bases"] += len(sequence.rstrip())
                batch.extend(header + sequence + separator + quality)
                if len(batch) >= 1024 * 1024:
                    destination.write(batch)
                    batch.clear()
            if batch:
                destination.write(batch)
    except Exception as error:
        errors.append(error)
    finally:
        destination.close()


class AlignmentStats:
    def __init__(self, reference):
        self.deltas = {name: array("i", [0]) * (length + 1)
                       for name, length in reference.items()}
        self.total = self.mapped = self.multiple = 0
        self.edits = self.columns = self.aligned_bases = 0
        self.input_bases = self.mapped_read_bases = self.soft_clipped = 0
        self.score_sum = self.retained_80 = 0
        self.edit_histogram = Counter()
        self.length_histogram = Counter()
        self.last_primary = None
        self.alternative_seen = False

    def add(self, line):
        fields = line.rstrip(b"\r\n").split(b"\t")
        flag = int(fields[1])
        if flag & 2048:
            raise ValueError("Unexpected supplementary alignment")
        if flag & 256:
            # Bowtie 2 emits each read's primary and secondary records together.
            if fields[0] != self.last_primary:
                raise ValueError("Secondary alignment is not adjacent to its primary")
            if not self.alternative_seen:
                self.multiple += 1
                self.alternative_seen = True
            return
        self.last_primary = fields[0]
        self.alternative_seen = False
        self.total += 1
        length = len(fields[9])
        self.input_bases += length
        self.length_histogram[length] += 1
        if flag & 4:
            return
        self.mapped += 1
        self.mapped_read_bases += length
        tags = {}
        for item in fields[11:]:
            name, value_type, value = item.split(b":", 2)
            tags[name] = value
        edits = int(tags[b"NM"])
        self.edits += edits
        self.edit_histogram[edits] += 1
        self.score_sum += int(tags[b"AS"])
        position = int(fields[3]) - 1
        delta = self.deltas[fields[2].decode("ascii")]
        retained = 0
        for size_text, operation in CIGAR.findall(fields[5].decode("ascii")):
            size = int(size_text)
            if operation in "M=X":
                delta[position] += 1
                delta[position + size] -= 1
                position += size
                retained += size
                self.columns += size
            elif operation == "I":
                retained += size
                self.columns += size
            elif operation == "D":
                position += size
                self.columns += size
            elif operation == "N":
                position += size
            elif operation == "S":
                self.soft_clipped += size
        self.aligned_bases += retained
        self.retained_80 += retained >= 0.8 * length

    def finish(self, output):
        genome_size = sum(len(delta) - 1 for delta in self.deltas.values())
        depth_sum = covered = depth_10 = depth_30 = 0
        hist = Counter()
        with (output / "coverage_10kb.tsv").open("w", encoding="utf-8") as table:
            table.write("reference\tstart_0based\tend_exclusive\tmean_depth\tcovered_fraction\n")
            for name, delta in self.deltas.items():
                depth = np.cumsum(np.asarray(delta, dtype=np.int32), dtype=np.int64)[:-1]
                if np.any(depth < 0):
                    raise ValueError("Negative coverage")
                depth_sum += int(depth.sum())
                covered += int(np.count_nonzero(depth))
                depth_10 += int(np.count_nonzero(depth >= 10))
                depth_30 += int(np.count_nonzero(depth >= 30))
                values, counts = np.unique(depth, return_counts=True)
                hist.update(dict(zip(map(int, values), map(int, counts))))
                for start in range(0, len(depth), 10000):
                    window = depth[start:start + 10000]
                    table.write(f"{name}\t{start}\t{start + len(window)}\t"
                                f"{window.mean():.6f}\t{np.count_nonzero(window)/len(window):.6f}\n")
        with (output / "depth_histogram.tsv").open("w") as table:
            table.write("depth\tgenome_positions\n")
            for depth, count in sorted(hist.items()):
                table.write(f"{depth}\t{count}\n")
        pct = lambda n, d: 100 * n / d if d else 0.0
        return {
            "total_reads": self.total, "mapped_reads": self.mapped,
            "unmapped_reads": self.total - self.mapped,
            "one_placement_reads": self.mapped - self.multiple,
            "multiple_placement_reads": self.multiple,
            "mapping_percent": pct(self.mapped, self.total),
            "one_placement_percent_all": pct(self.mapped - self.multiple, self.total),
            "multiple_placement_percent_all": pct(self.multiple, self.total),
            "one_placement_percent_mapped": pct(self.mapped - self.multiple, self.mapped),
            "multiple_placement_percent_mapped": pct(self.multiple, self.mapped),
            "read_length_histogram": dict(sorted(self.length_histogram.items())),
            "input_bases": self.input_bases, "genome_length": genome_size,
            "quality": {
                "edit_distance_sum": self.edits, "alignment_columns": self.columns,
                "pooled_identity_percent": pct(self.columns - self.edits, self.columns),
                "mean_edit_distance": self.edits / self.mapped if self.mapped else 0,
                "mean_alignment_score": self.score_sum / self.mapped if self.mapped else 0,
                "aligned_read_bases": self.aligned_bases,
                "mapped_read_bases": self.mapped_read_bases,
                "soft_clipped_bases": self.soft_clipped,
                "retained_bases_percent": pct(self.aligned_bases, self.mapped_read_bases),
                "reads_retaining_at_least_80_percent": self.retained_80,
                "reads_retaining_at_least_80_percent_of_mapped": pct(self.retained_80, self.mapped),
                "edit_distance_histogram": dict(sorted(self.edit_histogram.items())),
            },
            "coverage": {
                "covered_positions": covered, "covered_percent": pct(covered, genome_size),
                "depth_sum": depth_sum, "mean_depth": depth_sum / genome_size,
                "positions_at_least_10x": depth_10, "at_least_10x_percent": pct(depth_10, genome_size),
                "positions_at_least_30x": depth_30, "at_least_30x_percent": pct(depth_30, genome_size),
            },
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--aligner", default="bowtie2")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    command = [args.aligner, "--very-sensitive-local", "-k", "2", "--seed", "42",
               "--reorder", "-p", str(args.threads), "-x", str(args.index), "-U", "-"]
    started = datetime.now(timezone.utc).isoformat()
    start = time.perf_counter()
    stats = AlignmentStats(read_reference(args.reference))
    counts = Counter()
    errors = []
    print("Starting:", subprocess.list2cmdline(command), flush=True)
    with (args.output / "bowtie2.log").open("wb") as log:
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=log)
        feeder = threading.Thread(target=feed_reads, args=(args.reads, process.stdin, errors, counts))
        feeder.start()
        try:
            with (args.output / "alignment_sample.sam").open("wb") as sample:
                for line in process.stdout:
                    if line.startswith(b"@"):
                        sample.write(line)
                        continue
                    stats.add(line)
                    if stats.total <= 1000:
                        sample.write(line)
                    if stats.total and stats.total % 1000000 == 0 and not int(line.split(b"\t", 2)[1]) & 256:
                        print(f"Processed {stats.total:,} reads ({time.perf_counter() - start:.0f} s)", flush=True)
            feeder.join()
            returncode = process.wait()
            if errors:
                raise errors[0]
            if returncode:
                raise RuntimeError(f"Bowtie 2 exited with {returncode}; see bowtie2.log")
            if stats.total != counts["reads"] or stats.input_bases != counts["bases"]:
                raise ValueError(f"FASTQ/SAM counts differ: {counts}, {stats.total}, {stats.input_bases}")
        except BaseException:
            process.kill()
            feeder.join()
            process.wait()
            raise
    result = stats.finish(args.output)
    result.update({"started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
                   "elapsed_seconds": time.perf_counter() - start, "command": command,
                   "fastq_validated_reads": counts["reads"], "fastq_validated_bases": counts["bases"],
                   "reference_sha256": hashlib.sha256(args.reference.read_bytes()).hexdigest(),
                   "reads_file": args.reads.name, "reads_compressed_bytes": args.reads.stat().st_size,
                   "aligner_version": subprocess.run([args.aligner, "--version"], capture_output=True,
                                                      text=True, check=True).stdout.splitlines()[0]})
    (args.output / "metrics.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("total_reads", "mapping_percent", "elapsed_seconds")}), flush=True)


if __name__ == "__main__":
    main()
