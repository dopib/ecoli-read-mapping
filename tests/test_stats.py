"""Small hand-calculated alignments check the statistics independently of Bowtie."""

import gzip
import io
from pathlib import Path
import sys
import tempfile
import unittest
from collections import Counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from map_reads import AlignmentStats, feed_reads


def sam(name, flag, position, cigar, sequence, edits=0, score=10):
    return (f"{name}\t{flag}\tchr\t{position}\t20\t{cigar}\t*\t0\t0\t"
            f"{sequence}\t{'I' * len(sequence)}\tNM:i:{edits}\tAS:i:{score}\n").encode()


class StatsTests(unittest.TestCase):
    def finish(self, stats):
        with tempfile.TemporaryDirectory() as directory:
            return stats.finish(Path(directory))

    def test_alternatives_do_not_inflate_coverage(self):
        stats = AlignmentStats({"chr": 10})
        stats.add(sam("a", 0, 1, "4M", "ACGT"))
        stats.add(sam("a", 256, 5, "4M", "ACGT"))
        stats.add(sam("b", 16, 3, "4M", "ACGT", edits=1))
        stats.add(sam("c", 4, 0, "*", "NNNN"))
        result = self.finish(stats)
        self.assertEqual(result["total_reads"], 3)
        self.assertEqual(result["mapped_reads"], 2)
        self.assertEqual(result["one_placement_reads"], 1)
        self.assertEqual(result["multiple_placement_reads"], 1)
        self.assertEqual(result["coverage"]["covered_positions"], 6)
        self.assertEqual(result["coverage"]["depth_sum"], 8)
        self.assertAlmostEqual(result["coverage"]["mean_depth"], 0.8)
        self.assertAlmostEqual(result["quality"]["pooled_identity_percent"], 87.5)

    def test_clipping_insertions_and_deletions(self):
        stats = AlignmentStats({"chr": 10})
        # Two soft clips, five M bases, one inserted base, one deleted base.
        stats.add(sam("a", 0, 2, "1S2M1I1M1D2M1S", "AACCGGTT", edits=3))
        result = self.finish(stats)
        self.assertEqual(result["coverage"]["depth_sum"], 5)
        self.assertEqual(result["coverage"]["covered_positions"], 5)
        self.assertEqual(result["quality"]["aligned_read_bases"], 6)
        self.assertEqual(result["quality"]["soft_clipped_bases"], 2)
        self.assertEqual(result["quality"]["alignment_columns"], 7)
        self.assertAlmostEqual(result["quality"]["pooled_identity_percent"], 400 / 7)

    def test_unmapped_input_keeps_zero_depth_positions(self):
        stats = AlignmentStats({"chr": 10})
        stats.add(sam("a", 4, 0, "*", "NNNN"))
        result = self.finish(stats)
        self.assertEqual(result["coverage"]["mean_depth"], 0)
        self.assertEqual(result["mapping_percent"], 0)

    def test_reference_last_base_is_counted(self):
        stats = AlignmentStats({"chr": 10})
        stats.add(sam("a", 0, 10, "1M", "A"))
        self.assertEqual(self.finish(stats)["coverage"]["covered_positions"], 1)

    def test_fastq_validation(self):
        for contents, expected_reads, expected_errors in [
            (b"@a\nACGT\n+\nIIII\n", 1, 0),
            (b"@a\nACGT\n+\nIII\n", 0, 1),
            (b"@a\nACGT\n+\n", 0, 1),
        ]:
            with self.subTest(contents=contents), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "reads.fastq.gz"
                path.write_bytes(gzip.compress(contents))
                errors, counts = [], Counter()
                feed_reads(path, io.BytesIO(), errors, counts)
                self.assertEqual(counts["reads"], expected_reads)
                self.assertEqual(len(errors), expected_errors)


if __name__ == "__main__":
    unittest.main()
