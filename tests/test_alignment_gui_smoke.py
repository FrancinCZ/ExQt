from __future__ import annotations

import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))

from PySide6.QtWidgets import QApplication, QDialog, QWidget
from App import AlignmentOptionsDialog, ExQt


class AlignmentGuiSmokeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_alignment_dialog_and_window_actions(self):
        with tempfile.TemporaryDirectory() as input_folder:
            dialog = AlignmentOptionsDialog(input_folder)
            output_folder, reference_channel, align_masks, max_step_px, max_fail_pct, expand_canvas = dialog.values()
            self.assertTrue(output_folder)
            self.assertGreaterEqual(reference_channel, 0)
            self.assertTrue(align_masks)
            self.assertGreater(max_step_px, 0)
            self.assertGreater(max_fail_pct, 0)
            self.assertTrue(expand_canvas)

            # Auto-detect button and result label must exist
            self.assertTrue(hasattr(dialog, "auto_detect_btn"))
            self.assertTrue(dialog.auto_detect_btn.isEnabled())
            self.assertTrue(hasattr(dialog, "channel_info_label"))
            dialog.close()

        with patch("napari.Viewer") as mock_viewer:
            mock_inst = MagicMock()
            mock_inst.window._qt_window = QWidget()
            mock_viewer.return_value = mock_inst
            window = ExQt()
            self.assertEqual(window.align_stacks_action.text(), "Align Z-stacks...")
            self.assertTrue(window.align_stacks_action.isEnabled())
            self.assertTrue(hasattr(window, "align_progress_bar"))
            window.close()
            self.app.processEvents()

    def test_napari_layer_scale_for_3d_mode(self):
        import numpy as np
        import tifffile
        from Batch import process_condensates

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            raw_path = tmp / "sample.tif"
            mask_path = tmp / "sample_Mask.tif"

            # 4 Z-slices, 2 channels (C=0 DAPI, C=1 Signal), 20x20 YX
            stack = np.ones((4, 2, 20, 20), dtype=np.uint16) * 100
            mask = np.zeros((4, 20, 20), dtype=np.uint16)
            mask[1:3, 8:12, 8:12] = 1

            tifffile.imwrite(raw_path, stack, metadata={"axes": "ZCYX"})
            tifffile.imwrite(mask_path, mask, metadata={"axes": "ZYX"})

            received_layers = []
            process_condensates(
                tif_path=raw_path,
                mask_path=mask_path,
                mode="3d",
                pixel_size_nm=50.0,
                z_step_nm=200.0,
                signal_channel=1,
                dapi_channel=0,
                auto_roi=True,
                send_layer_func=received_layers.append,
            )

            # Check that 3D preview layers received scale=(200/50, 1, 1) = (4.0, 1, 1)
            image_layers = [l for l in received_layers if l.get("type") == "image"]
            self.assertTrue(len(image_layers) >= 1)
            for layer in image_layers:
                self.assertIn("kwargs", layer)
                self.assertEqual(layer["kwargs"].get("scale"), (4.0, 1, 1))

            label_layers = [l for l in received_layers if l.get("type") == "labels"]
            self.assertTrue(len(label_layers) >= 1)
            for layer in label_layers:
                self.assertIn("kwargs", layer)
                self.assertEqual(layer["kwargs"].get("scale"), (4.0, 1, 1))

    def test_mode_a_layer_mean_intensities(self):
        import numpy as np
        import tifffile
        from Batch import process_condensates

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            raw_path = tmp / "sample.tif"
            mask_path = tmp / "sample_Mask.tif"

            # 10 Z-slices, 2 channels, 30x30 YX
            stack = np.zeros((10, 2, 30, 30), dtype=np.uint16)
            # Create a 3D sphere/box object with graded intensity (core brighter than shell)
            mask = np.zeros((10, 30, 30), dtype=np.uint16)
            mask[2:8, 5:25, 5:25] = 1
            # Signal channel (C=1): realistic intensity structure:
            #   100 = camera dark level (baseline across entire FOV)
            #   400 = nucleoplasm signal (region around condensate within active Z)
            #   500 = condensate body (within mask)
            #   2000 = bright core of the condensate
            stack[:, 1, :, :] = 100         # camera dark level
            stack[2:8, 1, :, :] = 400       # nucleoplasm in active Z-slices
            stack[2:8, 1, 5:25, 5:25] = 500 # condensate body (overlaps mask)
            stack[4:6, 1, 12:18, 12:18] = 2000  # bright core of the condensate

            tifffile.imwrite(raw_path, stack, metadata={"axes": "ZCYX"})
            tifffile.imwrite(mask_path, mask, metadata={"axes": "ZYX"})

            df = process_condensates(
                tif_path=raw_path,
                mask_path=mask_path,
                mode="3d",
                pixel_size_nm=50.0,
                z_step_nm=50.0,
                signal_channel=1,
                dapi_channel=0,
                auto_roi=True,
                mode_a_enabled=True,
                mode_a_min_core_voxels=5,
            )

            self.assertIsNotNone(df)
            self.assertFalse(df.empty)
            self.assertIn("mean_intensity_core", df.columns)
            self.assertIn("mean_intensity_middle", df.columns)
            self.assertIn("mean_intensity_shell", df.columns)
            self.assertIn("Delta_intensity_core_shell", df.columns)

            # Core has the 2000 intensity voxels, so core mean > shell mean
            core_val = df["mean_intensity_core"].iloc[0]
            shell_val = df["mean_intensity_shell"].iloc[0]
            self.assertGreater(core_val, shell_val)

            # Test partition coefficient and condensate classification
            self.assertIn("partition_coefficient", df.columns)
            self.assertIn("condensate_class", df.columns)
            self.assertIn("nucleoplasm_mean_intensity", df.columns)
            part_val = df["partition_coefficient"].iloc[0]
            self.assertTrue(np.isfinite(part_val))
            self.assertGreater(part_val, 1.0)
            self.assertIn(df["condensate_class"].iloc[0], ["Globular (Low FA) / Core-Enriched", "Unconstrained Droplet (Spherical LLPS)"])

            # Test partitioning plots export
            from partitioning_plots import export_partitioning_analysis
            plots_dir = tmp / "partitioning_output"
            summary = export_partitioning_analysis(df, plots_dir, "test_sample")
            self.assertTrue((plots_dir / "test_sample_partitioning_analysis.png").is_file())
            self.assertTrue((plots_dir / "test_sample_partitioning_summary.json").is_file())
            self.assertGreater(summary["total_objects"], 0)



if __name__ == "__main__":
    unittest.main()



