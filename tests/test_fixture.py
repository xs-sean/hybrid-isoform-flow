import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from idp_beta.fixture import write_synthetic_fixture


class FixtureTests(unittest.TestCase):
    def test_write_synthetic_fixture_creates_public_toy_inputs(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            written = write_synthetic_fixture(root)
            self.assertGreater(len(written), 5)
            sample_sheet = root / "config" / "samples.tsv"
            self.assertTrue(sample_sheet.exists())
            self.assertIn("ToyCon", sample_sheet.read_text())
            all_text = "\n".join(path.read_text(errors="replace") for path in written if path.is_file())
            self.assertNotIn("PRIVATE_MANUSCRIPT_SENTINEL", all_text)
            self.assertNotIn("PRIVATE_ABSOLUTE_PATH_SENTINEL", all_text)


if __name__ == "__main__":
    unittest.main()
