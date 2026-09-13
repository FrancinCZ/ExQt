import json
from pathlib import Path
import tempfile
import unittest

from scale_reporting import load_scale_metadata, scale_correction_label


class ScaleReportingTests(unittest.TestCase):
    def test_factors_and_dimensions(self):
        for factor, volume, area in [(1, 1, 1), (8, 512, 64), (10, 1000, 100), (2.5, 15.625, 6.25)]:
            metadata = {"parameters": {"expansion_factor": factor}}
            for column, power, divisor in [("volume_bio_um3", 3, volume), ("area_bio_um2", 2, area)]:
                label = scale_correction_label(metadata, [column])
                self.assertIn(f"{factor:g}^{power} (= / {divisor:g})", label)

    def test_invalid_never_defaults_to_ten(self):
        for value in [None, 0, -2, True, "bad", float("nan"), float("inf"), 1e200]:
            self.assertIn("not verified", scale_correction_label({"parameters": {"expansion_factor": value}}, ["volume_bio_um3"]))
        self.assertIn("unavailable", scale_correction_label({}, []))

    def test_saved_sidecar_and_missing_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            csv = Path(directory) / "batch.csv"
            self.assertEqual(load_scale_metadata(csv), {})
            sidecar = csv.with_name("batch_metadata.json")
            sidecar.write_text(json.dumps({"parameters": {"expansion_factor": 8}}))
            self.assertIn("512", scale_correction_label(load_scale_metadata(csv), ["volume_bio_um3"]))
            sidecar.write_text("invalid")
            self.assertEqual(load_scale_metadata(csv), {})


if __name__ == "__main__":
    unittest.main()
