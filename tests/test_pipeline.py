import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from idp_beta.config import Sample
from idp_beta.pipeline import gpd_to_gtf, make_idp_configs, merge_expression_tables, merge_gpd_outputs, summarize_outputs


class PipelineTests(unittest.TestCase):
    def test_make_idp_configs_selects_long_read_key(self):
        with TemporaryDirectory() as tmp:
            sample = Sample("S1", "control", long_read=Path("lr.gpd"), sr_sam=Path("sr.sam"), sr_junction_bed=Path("sr.bed"))
            refs = {
                "genome_fasta": "genome.fa",
                "gmap_database": "gmap_index",
                "reference_gpd": "ref.gpd",
                "allref_gpd": "allref.gpd",
            }
            idp = {"Nthread": 1, "read_length": 150}
            written = make_idp_configs([sample], refs, idp, Path(tmp))
            text = written[0].read_text()
            self.assertIn("LR_gpd_pathfilename = lr.gpd", text)
            self.assertIn("SR_sam_pathfilename = sr.sam", text)

    def test_make_idp_configs_uses_workflow_outputs_for_raw_samples(self):
        with TemporaryDirectory() as tmp:
            sample = Sample("S1", "control", isoseq_reads=Path("flnc.fa"), short_read_1=Path("r1.fq"), short_read_2=Path("r2.fq"))
            refs = {
                "genome_fasta": "genome.fa",
                "gmap_database": "gmap_index",
                "reference_gpd": "ref.gpd",
                "allref_gpd": "allref.gpd",
            }
            written = make_idp_configs([sample], refs, {"Nthread": 1}, Path(tmp), chromosome="chr1")
            text = written[0].read_text()
            self.assertIn("LR_psl_pathfilename = results/long_read/gmap/S1.chr1.psl", text)
            self.assertIn("SR_jun_pathfilename = results/idp_inputs/S1.chr1.splicemap-like.junctions.bed", text)

    def test_merge_and_summarize_outputs(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            output_dir = root / "idp" / "S1"
            output_dir.mkdir(parents=True)
            for name in ["isoform_detection.gpd", "isoform_prediction.gpd", "isoform.gpd"]:
                (output_dir / name).write_text("gene\ttx\tchr1\t+\t0\t10\t0\t10\t1\t0,\t10,\n")
            merged = merge_gpd_outputs(root / "idp", root / "gpd")
            self.assertEqual(len(merged), 3)
            summary = root / "summary.tsv"
            summarize_outputs(root, summary)
            self.assertIn("isoform_detection.gpd", summary.read_text())

    def test_merge_expression_tables_fills_missing_transcripts(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "idp" / "S1.chr1").mkdir(parents=True)
            (root / "idp" / "S2.chr1").mkdir(parents=True)
            (root / "idp" / "S1.chr1" / "isoform.exp").write_text("tx1\t1\n")
            (root / "idp" / "S2.chr1" / "isoform.exp").write_text("tx2\t2\n")
            samples = [Sample("S1", "control"), Sample("S2", "ipf")]
            out = root / "expression.tsv"
            merge_expression_tables(root / "idp", out, samples)
            self.assertEqual(out.read_text().splitlines(), ["transcript_id\tS1\tS2", "tx1\t1\t0", "tx2\t0\t2"])

    def test_gpd_to_gtf_normalizes_idp_gpd_before_conversion(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            gpd_dir = root / "gpd"
            gpd_dir.mkdir()
            (gpd_dir / "sample.isoform.gpd").write_text("geneA\tisoA\tchr1\t+\t0\t10\t0\t10\t1\t0,\t10,\textra\n")
            commands = gpd_to_gtf(gpd_dir, root / "gtf", Path("genePredToGtf"), dry_run=True)
            self.assertEqual(commands, [["genePredToGtf", "file", str(root / "gtf" / "sample.isoform.genePred"), str(root / "gtf" / "sample.isoform.gtf")]])
            self.assertEqual((root / "gtf" / "sample.isoform.genePred").read_text(), "geneA\tchr1\t+\t0\t10\t0\t10\t1\t0,\t10,\n")


if __name__ == "__main__":
    unittest.main()
