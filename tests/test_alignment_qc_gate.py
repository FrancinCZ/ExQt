from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import tifffile

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))

from Batch import validate_alignment_qc


class AlignmentQCGateTests(unittest.TestCase):
    def test_non_pass_sidecar_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            image = folder / "sample.tif"
            tifffile.imwrite(image, np.ones((3, 16, 16), dtype=np.uint16), ome=True, metadata={"axes": "ZYX"})
            (folder / "sample_drift.csv").write_text(
                "status,step_shift_y_px,step_shift_x_px,step_magnitude_px\n"
                "FAIL,0,0,0\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "excluded from statistics"):
                validate_alignment_qc(image)

    def test_missing_sidecar_is_allowed_for_original_unaligned_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            image = Path(directory) / "original.tif"
            self.assertIsNone(validate_alignment_qc(image))


if __name__ == "__main__":
    unittest.main(verbosity=2)
