from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from scipy import ndimage as ndi
import tifffile

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))

from stack_aligner import AlignmentConfig, _read_tiff, align_tiff_folder, align_tiff_stack


def _plane(shape=(80, 92)):
    y, x = np.indices(shape, dtype=float)
    image = (
        3200 * np.exp(-((y - 20) ** 2 + (x - 28) ** 2) / 85)
        + 2700 * np.exp(-((y - 56) ** 2 + (x - 64) ** 2) / 115)
        + 1900 * np.exp(-((y - 40) ** 2 + (x - 48) ** 2) / 50)
    )
    return image.astype(np.uint16)


def _stack():
    base = _plane()
    positions = np.array([[0, 0], [1, -1], [2, -2], [2, -3]], dtype=float)
    reference = np.stack([ndi.shift(base * .7 + 400, p, order=1, mode="constant", prefilter=False) for p in positions]).astype(np.uint16)
    signal = np.stack([ndi.shift(base, p, order=1, mode="constant", prefilter=False) for p in positions]).astype(np.uint16)
    return np.stack([reference, signal]), positions


def _write(path, data, axes):
    tifffile.imwrite(path, data, ome=True, metadata={"axes": axes, "PhysicalSizeX": .058, "PhysicalSizeZ": .250})


class LightAlignerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.input = self.root / "input"
        self.output = self.root / "aligned"
        self.input.mkdir()

    def tearDown(self):
        self.tmp.cleanup()

    def test_one_stack_writes_tiff_mask_csv_and_plot(self):
        data, positions = _stack()
        source = self.input / "sample.tif"
        mask = self.input / "sample_Mask.tif"
        _write(source, data, "CZYX")
        _write(mask, (data[1] > 600).astype(np.uint8), "ZYX")

        result = align_tiff_stack(source, self.output, reference_channel=0, mask_path=mask)
        self.assertEqual(result.alignment.status, "PASS")
        self.assertTrue(result.output.is_file())
        self.assertTrue(result.mask_output.is_file())
        self.assertTrue(result.drift_csv.is_file())
        self.assertTrue(result.drift_plot.is_file())
        self.assertTrue(np.allclose(result.alignment.cumulative_shifts_yx, -positions, atol=.3))

        with tifffile.TiffFile(result.output) as tif:
            self.assertEqual(tif.series[0].axes, "CZYX")
            self.assertEqual(tif.asarray().shape, (2, 4, 82, 95))
            self.assertIn("PhysicalSizeX=\"0.058\"", tif.ome_metadata)
        with tifffile.TiffFile(result.mask_output) as tif:
            self.assertEqual(tif.series[0].axes, "ZYX")
            self.assertEqual(tif.asarray().shape, (4, 82, 95))
            self.assertTrue(set(np.unique(tif.asarray())).issubset({0, 1}))


    def test_folder_alignment_creates_one_summary(self):
        data, _ = _stack()
        _write(self.input / "a.tif", data, "CZYX")
        _write(self.input / "b.tiff", data, "CZYX")
        messages = []
        batch = align_tiff_folder(self.input, self.output, progress_callback=messages.append)
        self.assertEqual(len(batch.records), 2)
        self.assertEqual(set(batch.records["status"]), {"PASS"})
        self.assertTrue(batch.summary_csv.is_file())
        self.assertTrue(len(messages) >= 2)

    def test_review_is_exported_for_inspection(self):
        data, _ = _stack()
        source = self.input / "review.tif"
        _write(source, data, "CZYX")
        result = align_tiff_stack(
            source, self.output, reference_channel=0,
            config=AlignmentConfig(min_post_correlation=0.999999999),
        )
        self.assertEqual(result.alignment.status, "REVIEW")
        self.assertTrue(result.output.is_file())
        self.assertTrue(result.drift_csv.is_file())
        self.assertTrue(result.drift_plot.is_file())

    def test_bad_mask_is_rejected_before_any_image_is_written(self):
        data, _ = _stack()
        source = self.input / "bad_mask.tif"
        mask = self.input / "bad_mask_Mask.tif"
        _write(source, data, "CZYX")
        _write(mask, np.zeros((2, 3, 80, 92), dtype=np.uint8), "CZYX")
        with self.assertRaises(ValueError):
            align_tiff_stack(source, self.output, mask_path=mask)
        self.assertFalse((self.output / source.name).exists())

    def test_existing_output_is_never_overwritten(self):
        data, _ = _stack()
        source = self.input / "repeat.tif"
        _write(source, data, "CZYX")
        align_tiff_stack(source, self.output)
        with self.assertRaises(FileExistsError):
            align_tiff_stack(source, self.output)

    def test_unknown_axes_are_rejected(self):
        source = self.input / "unknown_axes.tif"
        # tifffile cannot write an invalid OME axes string, so this is covered by
        # the public reader contract through the explicit axis validation tests.
        data, _ = _stack()
        with self.assertRaises(ValueError):
            from stack_aligner import apply_xy_shifts
            apply_xy_shifts(data, "QZYX", np.zeros((4, 2)))

    def test_scifio_image_index_stack_is_treated_as_z(self):
        data = np.zeros((4, 7, 8), dtype=np.uint8)

        class FakeTiff:
            series = [SimpleNamespace(axes="IYX", asarray=lambda: data)]
            ome_metadata = None

            def __enter__(self):
                return self

            def __exit__(self, *_):
                return False

        with patch("stack_aligner.tifffile.TiffFile", return_value=FakeTiff()):
            mask, axes, _ = _read_tiff(self.input / "scifio_mask.tif")
        self.assertEqual(axes, "ZYX")
        self.assertIs(mask, data)

    def test_same_folder_is_rejected(self):
        data, _ = _stack()
        source = self.input / "sample.tif"
        _write(source, data, "CZYX")
        with self.assertRaises(ValueError):
            align_tiff_stack(source, self.input)

    def test_flat_stack_is_blocked_and_has_zero_drift(self):
        source = self.input / "flat.tif"
        _write(source, np.ones((3, 32, 32), dtype=np.uint16), "ZYX")
        result = align_tiff_stack(source, self.output)
        self.assertEqual(result.alignment.status, "FAIL")
        self.assertIsNone(result.output)
        self.assertFalse((self.output / source.name).exists())
        self.assertTrue(result.drift_csv.is_file())
        self.assertTrue(result.drift_plot.is_file())
        self.assertTrue(np.array_equal(result.alignment.cumulative_shifts_yx, np.zeros((3, 2))))


if __name__ == "__main__":
    unittest.main(verbosity=2)
