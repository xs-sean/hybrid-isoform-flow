import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from idp_beta.adapters import cigar_to_idp, convert_sam_to_idp_inputs, parse_cigar, query_length


class AdapterTests(unittest.TestCase):
    def test_parse_cigar_normalizes_equals_and_mismatch(self):
        parts = parse_cigar("5=1X4M10N10M")
        self.assertEqual(parts, [(5, "M"), (1, "M"), (4, "M"), (10, "N"), (10, "M")])
        self.assertEqual(query_length(parts), 20)

    def test_cigar_to_idp_rejects_soft_clipped_reads(self):
        self.assertIsNone(cigar_to_idp("5S95M"))
        self.assertEqual(cigar_to_idp("10M20N10M"), "10M20N10M")

    def test_convert_sam_to_idp_inputs_writes_sam_and_junction_bed(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            sam = root / "input.sam"
            sam.write_text(
                "@SQ\tSN:chrToy\tLN:100\n"
                "read1\t0\tchrToy\t1\t60\t10M20N10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\n"
                "read2\t256\tchrToy\t1\t60\t10M20N10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\n"
                "read3\t0\tchrOther\t1\t60\t20M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\n"
            )
            out_sam = root / "out.sam"
            out_bed = root / "out.bed"
            stats = convert_sam_to_idp_inputs(sam, out_sam, out_bed, read_length=20, chromosome="chrToy")
            self.assertEqual(stats.total_records, 3)
            self.assertEqual(stats.kept_records, 1)
            self.assertEqual(stats.junctions, 1)
            self.assertIn("10M20N10M", out_sam.read_text())
            bed_fields = out_bed.read_text().strip().split("\t")
            self.assertEqual(bed_fields[0], "chrToy")
            self.assertEqual(bed_fields[3], "junction_[1](1/0)")
            self.assertEqual(bed_fields[9], "2")


if __name__ == "__main__":
    unittest.main()
