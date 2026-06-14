import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from idp_beta.config import ConfigError, read_samples, read_simple_yaml, validate_configs


class ConfigTests(unittest.TestCase):
    def test_read_samples_rejects_duplicate_ids(self):
        with TemporaryDirectory() as tmp:
            sample_sheet = Path(tmp) / "samples.tsv"
            sample_sheet.write_text(
                "sample_id\tcondition\tlong_read\tsr_sam\tsr_junction_bed\n"
                "S1\tcontrol\tlr.fa\tsr.sam\tsr.bed\n"
                "S1\tipf\tlr2.fa\tsr2.sam\tsr2.bed\n"
            )
            with self.assertRaisesRegex(ConfigError, "duplicates sample_id"):
                read_samples(sample_sheet)

    def test_read_samples_allows_optional_tpm(self):
        with TemporaryDirectory() as tmp:
            sample_sheet = Path(tmp) / "samples.tsv"
            sample_sheet.write_text(
                "sample_id\tcondition\tisoseq_reads\tshort_read_1\tshort_read_2\texternal_tpm\n"
                "S1\tcontrol\tflnc.fa\tr1.fq\tr2.fq\tmatrix.tsv\n"
            )
            samples = read_samples(sample_sheet)
            self.assertEqual(samples[0].sample_id, "S1")
            self.assertEqual(samples[0].external_tpm, Path("matrix.tsv"))
            self.assertTrue(samples[0].has_raw_inputs)

    def test_read_samples_keeps_processed_input_compatibility(self):
        with TemporaryDirectory() as tmp:
            sample_sheet = Path(tmp) / "samples.tsv"
            sample_sheet.write_text(
                "sample_id\tcondition\tlong_read\tsr_sam\tsr_junction_bed\ttpm\n"
                "S1\tcontrol\tlr.gpd\tsr.sam\tsr.bed\tmatrix.tsv\n"
            )
            samples = read_samples(sample_sheet)
            self.assertEqual(samples[0].long_read, Path("lr.gpd"))
            self.assertEqual(samples[0].sr_sam, Path("sr.sam"))
            self.assertEqual(samples[0].sr_junction_bed, Path("sr.bed"))
            self.assertTrue(samples[0].has_processed_inputs)

    def test_simple_yaml_nested_values(self):
        with TemporaryDirectory() as tmp:
            config = Path(tmp) / "config.yml"
            config.write_text("root:\n  answer: 42\nflag: true\nname: test\nitems: [chr1, chr2]\nempty:\n")
            self.assertEqual(read_simple_yaml(config), {"root": {"answer": 42}, "flag": True, "name": "test", "items": ["chr1", "chr2"], "empty": ""})

    def test_validate_configs_rejects_non_grch38_build(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            samples = root / "samples.tsv"
            refs = root / "references.yml"
            idp = root / "idp.yml"
            samples.write_text(
                "sample_id\tcondition\tisoseq_reads\tshort_read_1\tshort_read_2\n"
                "S1\tcontrol\tflnc.fa\tr1.fq\tr2.fq\n"
            )
            refs.write_text(
                "build: hg19\n"
                "genome_fasta: genome.fa\n"
                "chromosome_sizes: chrom.sizes\n"
                "hisat2_index: hisat2/index\n"
                "gmap_database: gmap/db\n"
                "reference_gpd: ref.gpd\n"
                "allref_gpd: allref.gpd\n"
                "genePredToGtf: genePredToGtf\n"
            )
            idp.write_text(
                "Nthread: 10\n"
                "Njun_limit: 10\n"
                "Niso_limit: 20\n"
                "L_exon_limit: 1700\n"
                "L_min_intron: 68\n"
                "read_length: 150\n"
                "min_junction_overlap_len: 10\n"
                "estimator_choice: MLE\n"
            )
            with self.assertRaisesRegex(ConfigError, "GRCh38/hg38"):
                validate_configs(samples, refs, idp, check_files=False)


if __name__ == "__main__":
    unittest.main()
