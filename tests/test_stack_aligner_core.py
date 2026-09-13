from __future__ import annotations

import itertools
from pathlib import Path
import sys
import unittest

import numpy as np

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))

from stack_aligner import AlignmentConfig, apply_xy_shifts, estimate_xy_drift, extract_reference_stack


class StackAlignerCoreTests(unittest.TestCase):
    def test_channel_extraction_and_shift_application_for_all_czyx_axis_orders(self):
        rng = np.random.default_rng(42)
        zyx = rng.normal(size=(3, 15, 17)).astype(np.float32)
        shifts = np.array([[0.0, 0.0], [1.0, -2.0], [-1.0, 1.0]])

        for axis_tuple in itertools.permutations("CZYX"):
            axes = "".join(axis_tuple)
            shape = tuple({"C": 2, "Z": 3, "Y": 15, "X": 17}[axis] for axis in axes)
            data = np.empty(shape, dtype=np.float32)
            for channel in range(2):
                indexer = [slice(None)] * 4
                indexer[axes.index("C")] = channel
                data[tuple(indexer)] = np.transpose(
                    zyx + channel * 100.0,
                    ["ZYX".index(axis) for axis in axes.replace("C", "")],
                )

            extracted = extract_reference_stack(data, axes, channel_index=1)
            self.assertTrue(np.array_equal(extracted, zyx + 100.0), axes)
            shifted = apply_xy_shifts(data, axes, shifts, interpolation_order=0)
            roundtrip = apply_xy_shifts(shifted, axes, -shifts, interpolation_order=0)
            central = [slice(None)] * 4
            central[axes.index("Y")] = slice(3, -3)
            central[axes.index("X")] = slice(3, -3)
            self.assertTrue(np.array_equal(roundtrip[tuple(central)], data[tuple(central)]), axes)

    def test_textureless_stack_fails_without_false_shifts(self):
        result = estimate_xy_drift(np.zeros((3, 48, 48), dtype=np.uint16))
        self.assertEqual(result.status, "FAIL")
        self.assertEqual(set(result.metrics["reason"]), {"insufficient_texture"})
        self.assertTrue(np.array_equal(result.cumulative_shifts_yx, np.zeros((3, 2))))

    def test_time_series_and_non_finite_shifts_are_rejected(self):
        array = np.zeros((2, 1, 3, 20, 20), dtype=np.uint16)
        with self.assertRaises(ValueError):
            extract_reference_stack(array, "TCZYX", channel_index=0)
        with self.assertRaises(ValueError):
            apply_xy_shifts(
                np.zeros((3, 20, 20), dtype=np.uint16),
                "ZYX",
                np.array([[0.0, 0.0], [np.nan, 0.0], [0.0, 0.0]]),
            )

    def test_expand_canvas_preserves_total_intensity_and_dimensions(self):
        # 3 planes, 15x15 image with foreground block
        img = np.zeros((3, 15, 15), dtype=np.float32)
        img[:, 4:10, 4:10] = 50.0
        shifts = np.array([[0.0, 0.0], [2.0, -3.0], [-1.0, 4.0]])

        # Default unexpanded: same shape
        unexpanded = apply_xy_shifts(img, "ZYX", shifts, expand_canvas=False)
        self.assertEqual(unexpanded.shape, (3, 15, 15))

        # Expanded: canvas expands by max shifts in Y and X
        expanded = apply_xy_shifts(img, "ZYX", shifts, expand_canvas=True)
        # Y shifts: min=-1 (pad_top=1), max=2 (pad_bottom=2) -> H = 15 + 1 + 2 = 18
        # X shifts: min=-3 (pad_left=3), max=4 (pad_right=4) -> W = 15 + 3 + 4 = 22
        self.assertEqual(expanded.shape, (3, 18, 22))
        # Total sum of intensity is fully preserved
        for z in range(3):
            self.assertAlmostEqual(float(img[z].sum()), float(expanded[z].sum()), places=2)

    def test_residual_misalignment_after_translation_is_flagged(self):
        # Two planes related by a 1px shift but with added structured noise so that
        # correlation_after stays well below min_residual_correlation=0.99.
        # This tests that the detection branch fires when translation succeeds
        # (shift within max_step_px, finite, bidirectional OK) but residual
        # correlation is still too low — the expected real-world trigger for
        # gel tilt, rotation, or local warp that pure XY translation cannot fix.
        rng = np.random.default_rng(7)
        base = rng.standard_normal((64, 64)).astype(np.float32)
        from scipy.ndimage import shift as ndi_shift
        shifted = ndi_shift(base, [1.0, 0.0], mode="wrap").astype(np.float32)
        # 30% additive noise lowers correlation_after to ~0.96, below the 0.99 threshold
        shifted += rng.standard_normal((64, 64)).astype(np.float32) * 0.3

        stack = np.stack([base, shifted, base])
        config = AlignmentConfig(
            min_residual_correlation=0.99,
            min_post_correlation=0.05,
            max_step_px=8.0,
        )
        result = estimate_xy_drift(stack, config)

        reasons_z0_z1 = [r["reason"] for r in result.metrics.to_dict("records") if r["reference_z"] == 0]
        self.assertTrue(
            any("residual_misalignment_after_translation" in r for r in reasons_z0_z1),
            "Expected residual_misalignment flag, got: %s" % reasons_z0_z1,
        )
        # Translation was accepted so overall is REVIEW, never FAIL
        self.assertEqual(result.status, "REVIEW")



if __name__ == "__main__":
    unittest.main(verbosity=2)

