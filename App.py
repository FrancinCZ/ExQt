import os
import json
import threading
from datetime import datetime
from Batch import (
    _interpolate_or_extrude_roi,
    get_metadata_from_tif,
    matching_mask_path as _matching_mask_path,
    process_condensates,
    roi_output_path,
    source_tiff_files as _source_tiff_files,
)
from rezim_a_metrics import GRADIENT_ALPHA, GRADIENT_TEST, MODE_A_LAYER_SCHEME
from provenance import code_provenance, file_sha256, input_file_record
from reference_values import REFERENCE_VALUES
from defaults import DEFAULT_SETTINGS
from lif_metadata import offset_from_lif
from stack_aligner import AlignmentConfig, align_tiff_folder
import sys
import napari
from PySide6.QtWidgets import (QApplication, QLabel, QMainWindow, QPushButton, 
                                QVBoxLayout, QHBoxLayout, QWidget, QFormLayout, 
                                QSpinBox, QDoubleSpinBox, QComboBox, QCheckBox, QLineEdit, QFileDialog, QDialog, QDialogButtonBox,
                                QMessageBox, QProgressBar, QSplitter)
from PySide6.QtGui import QAction
import qdarktheme
from PySide6.QtCore import QSettings, QThread, Signal, Qt
import numpy as np
import tifffile
import pandas as pd
from pathlib import Path
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from postprocessing import (
    QCPolicyMismatchError,
    _prepare_reporting_frames,
    detect_red_flags,
    generate_excel_stats,
    generate_plots,
    generate_rezim_a_plots,
    merge_statistics_folder,
)
from calibration_policy import (
    CALIBRATION_SOURCE_GUI,
    CALIBRATION_SOURCE_TIFF,
    looks_pre_expansion,
    resolve_file_calibration,
    summarize_calibrations,
    validate_calibration,
)
from size_preview import collect_size_preview


REPORT_DEFAULTS = {
    "report_excel": True,
    "report_primary_csv": True,
    "report_excluded_csv": True,
    "report_raw_audit_csv": False,
    "report_standard_plots": True,
    "report_mode_a_plots": True,
    "report_partitioning_plots": True,
}
#Radial FA Profiling settings passed unchanged from DEFAULT_SETTINGS/params to Batch.
MODE_A_KEYS = (
    "mode_a_min_core_voxels", "mode_a_exclude_split_slices",
    "mode_a_z_split_min_component_voxels", "mode_a_z_split_min_component_fraction",
    "mode_a_z_split_pass_fraction", "mode_a_z_split_review_fraction",
)


def _safe_float(value, default):
    if value is None:
        return default
    try:
        if isinstance(value, str):
            value = value.replace(',', '.')
        return float(value)
    except (ValueError, TypeError):
        return default


def apply_theme(dark):
    app = QApplication.instance()
    if app is not None:
        app.setStyleSheet(qdarktheme.load_stylesheet("dark" if dark else "light"))


def _safe_bool(value, default=False):
    if value is None:
        return default
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


#Offset source stored by AdvancedSettingsDialog; anything unreadable counts as a manual value.
def _load_offset_source(settings):
    try:
        source = json.loads(settings.value("adv_detector_offset_source") or "")
    except (TypeError, ValueError):
        return {"method": "manual"}
    return source if isinstance(source, dict) and "method" in source else {"method": "manual"}


class AnalysisAborted(Exception):
    """Raised inside the worker when the user stops the batch during ROI drawing."""


class AnalysisWorker(QThread):
    #Run batch analysis off the GUI thread and relay UI-safe signals.
    layer_ready = Signal(dict)
    progress = Signal(str)
    request_roi_signal = Signal(dict)
    request_review_signal = Signal()

    def __init__(self, params):
        super().__init__()
        self.params = params
        self.roi_event = threading.Event()
        self.review_event = threading.Event()
        self.user_roi_data = None
        self.user_roi_source = "manual"
        self.abort_requested = False

    #Pause the worker until the GUI returns a painted ROI mask.
    def request_roi_callback(self, img_shape, is_3d):
        # Reset first so an ROI from the previous image can never be returned for this one.
        self.user_roi_data = None
        self.user_roi_source = "manual"
        self.request_roi_signal.emit({"shape": img_shape, "is_3d": is_3d})
        self.roi_event.wait()
        self.roi_event.clear()
        if self.abort_requested:
            raise AnalysisAborted()
        if self.user_roi_data is None:
            raise RuntimeError("No ROI was received from the GUI for this image.")
        return self.user_roi_data

    #Wake any pending ROI/review wait and mark the batch for shutdown.
    def request_abort(self):
        self.abort_requested = True
        self.review_event.set()
        self.roi_event.set() 

    def run(self):
        p = self.params
        try:
            folder_path = Path(p["input_folder"])
            if not folder_path.is_dir():
                self.progress.emit("Error: Input folder does not exist or is not valid.")
                self.progress.emit("Done")
                return
            output_folder = Path(p["output_folder"]) if p["output_folder"] else folder_path
            frames, calibrations, input_records, excluded_files = [], {}, {}, []

            def discard(name):
                self.progress.emit(f"Discarded and stopped at: {name}")
                excluded_files.append({"name": name, "reason": "discarded by user"})

            for raw_tif in _source_tiff_files(folder_path):
                if self.abort_requested:
                    self.progress.emit("Analysis stopped by user.")
                    break
                mask_file = _matching_mask_path(raw_tif)
                if not mask_file.exists():
                    self.progress.emit(f"Skipping {raw_tif.name}: Mask not found.")
                    excluded_files.append({"name": raw_tif.name, "reason": "mask not found"})
                    continue
                try:
                    #The GUI's calibration choice decides; the worker never re-reads TIFF tags.
                    calibration = resolve_file_calibration(
                        raw_tif.name, p.get("calibration_source", CALIBRATION_SOURCE_GUI),
                        p["pixel_size_nm"], p["z_step_nm"], p.get("detected_metadata_by_file", {}),
                    )
                    self.progress.emit(
                        f"Processing: {raw_tif.name} (XY={calibration['pixel_size_nm']:g} nm, "
                        f"Z={calibration['z_step_nm']:g} nm, source={calibration['source']})"
                    )
                    df_file = process_condensates(
                        tif_path=raw_tif, mask_path=mask_file, mode=p["mode"],
                        expansion_factor=p["expansion_factor"],
                        min_voxels=p.get("min_voxels", DEFAULT_SETTINGS["raw_min_voxels"]),
                        auto_roi=p["auto_roi"],
                        send_layer_func=self.layer_ready.emit if p.get("show_napari", True) else None,
                        request_roi_func=None if p.get("auto_roi", False) else self.request_roi_callback,
                        pixel_size_nm=calibration["pixel_size_nm"], z_step_nm=calibration["z_step_nm"],
                        signal_channel=p["signal_channel"], dapi_channel=p["dapi_channel"],
                        mode_a_enabled=p.get("mode_a_enabled", False),
                        **{key: p.get(key, DEFAULT_SETTINGS[key]) for key in MODE_A_KEYS},
                        detector_offset_adu=p.get("detector_offset_adu", DEFAULT_SETTINGS["adv_detector_offset"]),
                        detector_offset_source=p.get("detector_offset_source", {}).get("method", "manual"),
                        roi_output_dir=output_folder,
                    )
                    #Abort is honoured regardless of review: the current image is never kept.
                    if self.abort_requested:
                        discard(raw_tif.name)
                        break
                    measured = df_file is not None and not df_file.empty
                    if measured:
                        if not p.get("auto_roi", False):
                            df_file["roi_source"] = self.user_roi_source
                        record = input_file_record(raw_tif, mask_file)
                        roi_path = roi_output_path(raw_tif, mask_file, output_folder)
                        if roi_path.is_file():
                            record.update(roi_file=roi_path.name, roi_sha256=file_sha256(roi_path))
                    if p.get("review_each_image", False):
                        self.progress.emit(f"Review: {raw_tif.name}")
                        self.request_review_signal.emit()
                        self.review_event.wait()
                        self.review_event.clear()
                        if self.abort_requested:
                            discard(raw_tif.name)
                            break
                    if measured:
                        frames.append(df_file)
                        calibrations[raw_tif.name] = calibration
                        input_records[raw_tif.name] = record
                except AnalysisAborted:
                    discard(raw_tif.name)
                    break
                except Exception as e:
                    self.progress.emit(f"Error: {raw_tif.name}: {e}")
                    excluded_files.append({"name": raw_tif.name, "reason": str(e)})

            if frames:
                self._write_outputs(folder_path, output_folder, frames, calibrations, input_records, excluded_files)
            else:
                self.progress.emit("No data available for processing.")
            self.progress.emit("Done")
        except Exception as e:
            self.progress.emit(f"Error when loading: {e}")
            self.progress.emit("Done")

    #Save the batch CSV and metadata.json, then the reports the user selected.
    def _write_outputs(self, folder_path, output_folder, frames, calibrations, input_records, excluded_files):
        p = self.params
        final_df = pd.concat(frames, ignore_index=True)
        output_folder.mkdir(parents=True, exist_ok=True)
        output_csv = output_folder / f"{folder_path.name}_Output_Batch_{p['mode']}.csv"
        final_df.to_csv(output_csv, index=False)
        self.progress.emit(f"CSV saved: {output_csv.name}")

        # Red flags are reported loudly and never used to filter.
        red_flags = detect_red_flags(final_df)
        for flag in red_flags:
            self.progress.emit(f"RED FLAG: {flag['message']}")

        metadata = self._run_metadata(final_df, red_flags, calibrations, input_records, excluded_files)
        try:
            with open(output_folder / f"{output_csv.stem}_metadata.json", "w", encoding="utf-8") as f:
                json.dump(metadata, f, indent=4)
            self.progress.emit("Metadata saved.")
        except Exception as e:
            self.progress.emit(f"Error saving metadata: {e}")

        min_size = p.get("plot_min_size", DEFAULT_SETTINGS["plot_min_size"])
        max_size = p.get("plot_max_size", DEFAULT_SETTINGS["plot_max_size"])

        def wanted(key):
            return p.get("generate_reports", False) and p.get(key, REPORT_DEFAULTS[key])

        def report(error_label, action):
            try:
                self.progress.emit(action())
            except Exception as e:
                self.progress.emit(f"Error generating {error_label}: {e}")

        table_keys = ("report_excel", "report_primary_csv", "report_excluded_csv", "report_raw_audit_csv")
        if any(wanted(key) for key in table_keys):
            def tables():
                result = generate_excel_stats(
                    str(output_csv), min_size=min_size, max_size=max_size,
                    generate_excel=p.get("report_excel", True),
                    generate_primary_csv=p.get("report_primary_csv", True),
                    generate_excluded_csv=p.get("report_excluded_csv", True),
                    generate_raw_audit_csv=p.get("report_raw_audit_csv", False),
                )
                names = [n for n in ("excel", "primary_csv", "excluded_csv", "all_objects_csv") if result.get(n)]
                return f"Selected tables generated: {', '.join(names)}"
            report("report tables", tables)
        if wanted("report_standard_plots"):
            def standard_plots():
                generate_plots(str(output_csv), min_size=min_size, max_size=max_size)
                return "Graphs were generated."
            report("graphs", standard_plots)
        if wanted("report_mode_a_plots") and p.get("mode_a_enabled") and p["mode"] == "3d":
            def mode_a_plots():
                generate_rezim_a_plots(str(output_csv), min_size=min_size, max_size=max_size)
                return "Radial FA Profiling plots were generated."
            report("Radial FA Profiling plots", mode_a_plots)
        if wanted("report_partitioning_plots"):
            def partitioning():
                from partitioning_plots import export_partitioning_analysis
                # The report must see primary_qc_valid, which only exists after QC preparation.
                qc_frame = _prepare_reporting_frames(final_df, min_size, max_size)[3]
                export_partitioning_analysis(qc_frame, output_folder, file_stem=output_csv.stem,
                                            min_size=min_size, max_size=max_size)
                return "Partitioning & Classification plots (K_part) were generated."
            report("Partitioning plots", partitioning)

    #Everything needed to reproduce and interpret the run.
    def _run_metadata(self, final_df, red_flags, calibrations, input_records, excluded_files):
        p = self.params
        #Batch-level calibration only when every file used the same one; per-file truth is below.
        pairs = {(c["pixel_size_nm"], c["z_step_nm"]) for c in calibrations.values()}
        uniform_xy, uniform_z = pairs.pop() if len(pairs) == 1 else (None, None)
        setting = lambda key: p.get(key, DEFAULT_SETTINGS[key])
        review_files = (
            sorted(final_df.loc[final_df["alignment_status"].astype(str).eq("REVIEW"), "filename"].unique().tolist())
            if "alignment_status" in final_df.columns else []
        )
        return {
            "timestamp": datetime.now().isoformat(),
            "software": "ExQt",
            "provenance": code_provenance(),
            "red_flags": red_flags,
            "reference_values": REFERENCE_VALUES,
            "parameters": {
                "mode": p["mode"],
                "auto_roi": p.get("auto_roi"),
                "expansion_factor": p["expansion_factor"],
                "min_voxels": p.get("min_voxels", DEFAULT_SETTINGS["raw_min_voxels"]),
                "pixel_size_nm": uniform_xy,
                "z_step_nm": uniform_z,
                "gui_pixel_size_nm": p["pixel_size_nm"],
                "gui_z_step_nm": p["z_step_nm"],
                "calibration_source": p.get("calibration_source", CALIBRATION_SOURCE_GUI),
                "applied_calibration_by_file": calibrations,
                "calibration_confirmation": p.get("calibration_confirmation", "unknown"),
                "detected_metadata_by_file": p.get("detected_metadata_by_file", {}),
                "plot_min_size": setting("plot_min_size"),
                "plot_max_size": setting("plot_max_size"),
            },
            "files": {
                "analysed_files": list(calibrations),
                "input_files": list(input_records.values()),
                "alignment_review_files": review_files,
                "excluded_files": excluded_files,
            },
            "classification": {"gradient_test": GRADIENT_TEST, "gradient_alpha": GRADIENT_ALPHA},
            "partitioning": {
                "K_offset_method": "explicit_setting",
                "detector_offset_adu": p.get("detector_offset_adu", DEFAULT_SETTINGS["adv_detector_offset"]),
                "offset_source": p.get("detector_offset_source", {"method": "manual"}),
                "auto_roi": p.get("auto_roi"),
                "auto_roi_K_policy": "K not computed (NaN, K_invalid_reason=auto_roi_no_nucleoplasm)",
                "padding_policy": "voxels flagged 0 in <stem>_alignment_valid.tif are excluded from nucleoplasm statistics",
            },
            "channels": {"signal_channel": p["signal_channel"], "dapi_channel": p["dapi_channel"]},
            "mode_a": {
                "enabled": p.get("mode_a_enabled", False),
                "available_in_mode": "3d",
                "min_core_voxels": setting("mode_a_min_core_voxels"),
                "minimum_voxel_policy": "same minimum applied to object, shell, middle, and core FA validity; also used as minimum core size",
                "exclude_split_slices_from_primary": setting("mode_a_exclude_split_slices"),
                "require_z_topology_pass_for_primary": setting("mode_a_exclude_split_slices"),
                "z_split_policy": "substantial_component_fraction_v1",
                **{key.replace("mode_a_", ""): setting(key) for key in MODE_A_KEYS if key.startswith("mode_a_z_split_")},
                "layer_scheme": MODE_A_LAYER_SCHEME,
                "sampling_order": "Z,Y,X",
                "metrics": [
                    "A_object", "A_shell", "A_middle", "A_core",
                    "Delta_A_middle_shell", "Delta_A_core_middle", "Delta_A_core_shell",
                ],
                "interpretation": "Geometric structural-response metrics; not direct stiffness, liquidity, or viscosity measurements.",
            },
            "reporting": {
                "enabled": p.get("generate_reports", False),
                **{key.replace("report_", ""): p.get(key, default) for key, default in REPORT_DEFAULTS.items()},
            },
        }

class AlignmentWorker(QThread):
    """Run the small TIFF aligner without freezing the GUI."""
    progress = Signal(int, str)
    completed = Signal(object)
    failed = Signal(str)

    def __init__(self, input_folder, output_folder, reference_channel, align_masks, config=None):
        super().__init__()
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.reference_channel = reference_channel
        self.align_masks = align_masks
        self.config = config

    def run(self):
        try:
            result = align_tiff_folder(
                self.input_folder,
                self.output_folder,
                reference_channel=self.reference_channel,
                align_masks=self.align_masks,
                config=self.config,
                progress_callback=self.progress.emit,
            )
            self.completed.emit(result)
        except Exception as error:
            self.failed.emit(str(error))


class ChannelDetectionWorker(QThread):
    result = Signal(list)   #[(channel_index, mean_corr_pct), ...]
    error  = Signal(str)

    def __init__(self, input_folder: str, n_sample: int = 20):
        super().__init__()
        self.input_folder = input_folder
        self.n_sample      = n_sample

    def run(self):
        try:
            import tifffile as _tifffile
            import numpy as _np
            from scipy.ndimage import gaussian_filter as _gf
            from pathlib import Path as _Path

            folder  = _Path(self.input_folder)
            sources = sorted(
                p for p in folder.iterdir()
                if p.is_file() and p.suffix.lower() in {".tif", ".tiff"}
                and not p.stem.lower().endswith("_mask")
            )
            if not sources:
                self.error.emit("No TIFF files found in the input folder.")
                return

            with _tifffile.TiffFile(sources[0]) as tif:
                series = tif.series[0]
                data   = series.asarray()
                axes   = series.axes

            if "C" not in axes or "Z" not in axes:
                self.error.emit(
                    f"Cannot auto-detect: axes '{axes}' must contain both C and Z."
                )
                return

            c_ax = axes.index("C")
            z_ax = axes.index("Z")
            n_c  = data.shape[c_ax]
            n_z  = data.shape[z_ax]

            #Sample planes from the middle 80 % to skip empty leading/trailing edges.
            start = int(n_z * 0.10)
            end   = int(n_z * 0.90)
            step  = max(1, (end - start) // self.n_sample)
            z_idx = list(range(start, end, step))[: self.n_sample]
            if len(z_idx) < 2:
                self.error.emit("Stack too short for channel detection.")
                return

            def _prep(plane, hp_sigma=8.0):
                img = plane.astype(_np.float64)
                p1, p99 = _np.percentile(img, [1, 99])
                if p99 <= p1:
                    return None
                img = _np.clip(img, p1, p99)
                img = img - _gf(img, hp_sigma)
                std = img.std()
                if std < 1e-6:
                    return None
                return (img - img.mean()) / std

            def _pearson(a, b):
                a = a - a.mean(); b = b - b.mean()
                d = _np.sqrt((a * a).sum() * (b * b).sum())
                return float((a * b).sum() / d) if d > 0 else _np.nan

            results = []
            for ch in range(n_c):
                corrs = []
                for zi in z_idx:
                    #Compare truly adjacent planes (zi and zi+1), same as the aligner does.
                    if zi + 1 >= n_z:
                        continue
                    idx_a = [slice(None)] * data.ndim
                    idx_b = [slice(None)] * data.ndim
                    idx_a[c_ax] = ch; idx_a[z_ax] = zi
                    idx_b[c_ax] = ch; idx_b[z_ax] = zi + 1
                    prep_a = _prep(data[tuple(idx_a)])
                    prep_b = _prep(data[tuple(idx_b)])
                    if prep_a is not None and prep_b is not None:
                        v = _pearson(prep_a, prep_b)
                        if _np.isfinite(v):
                            corrs.append(v)
                mean_pct = float(_np.mean(corrs)) * 100.0 if corrs else 0.0
                results.append((ch, mean_pct))


            self.result.emit(results)
        except Exception as exc:
            self.error.emit(str(exc))


class AlignmentOptionsDialog(QDialog):
    """Keep alignment settings deliberately small: output, reference and masks."""
    def __init__(self, input_folder, parent=None):
        super().__init__(parent)
        self.input_folder = input_folder
        self.setWindowTitle("Align Z-stacks")
        self.setMinimumWidth(520)
        layout = QVBoxLayout(self)
        explanation = QLabel(
            "The aligner estimates XY translation from one reference channel and applies the same "
            "correction to all channels. Choose the channel with the highest inter-plane correlation "
            "— in ExM datasets this is often the signal channel, not DAPI, because expansion makes "
            "DAPI sparse and unreliable for phase correlation. "
            "Use Auto-detect to measure correlations automatically. "
            "Original TIFF files remain unchanged."
        )
        explanation.setWordWrap(True)
        layout.addWidget(explanation)

        form = QFormLayout()
        self.reference_channel_spin = QSpinBox()
        self.reference_channel_spin.setRange(0, 99)
        self.reference_channel_spin.setValue(0)
        self.reference_channel_spin.setToolTip(
            "Use the channel with the most continuous, uniform signal throughout the Z stack. "
            "For ExM data, DAPI is often a poor choice (sparse after expansion). "
            "Click 'Auto-detect' to measure inter-plane correlation for each channel automatically."
        )

        ref_row = QHBoxLayout()
        ref_row.addWidget(self.reference_channel_spin)
        self.auto_detect_btn = QPushButton("Auto-detect…")
        self.auto_detect_btn.setToolTip(
            "Measure average inter-plane correlation for each channel using the same "
            "preprocessing as the aligner, then automatically select the best channel."
        )
        self.auto_detect_btn.clicked.connect(self.run_auto_detect)
        ref_row.addWidget(self.auto_detect_btn)
        form.addRow("Reference channel (0-based):", ref_row)

        self.channel_info_label = QLabel("")
        self.channel_info_label.setWordWrap(True)
        self.channel_info_label.setStyleSheet("color: gray; font-size: 10px;")
        form.addRow("", self.channel_info_label)

        self.max_step_spin = QDoubleSpinBox()
        self.max_step_spin.setRange(1.0, 50.0)
        self.max_step_spin.setSingleStep(0.5)
        self.max_step_spin.setDecimals(1)
        self.max_step_spin.setValue(8.0)
        self.max_step_spin.setToolTip(
            "Maximum allowed XY shift between adjacent Z planes (pixels). "
            "Steps larger than this are rejected (shift held at zero for that slice). "
            "Increase if the aligner reports many 'step_above_limit' failures, "
            "but inspect the drift CSV first to confirm the large steps are real drift."
        )
        form.addRow("Max step per Z plane (px):", self.max_step_spin)

        self.max_fail_spin = QSpinBox()
        self.max_fail_spin.setRange(5, 100)
        self.max_fail_spin.setSingleStep(5)
        self.max_fail_spin.setValue(20)
        self.max_fail_spin.setSuffix(" %")
        self.max_fail_spin.setToolTip(
            "Maximum fraction of mid-stack steps that may fail before the file is blocked from export. "
            "Empty slices at the top/bottom of the stack (before/after the sample) are excluded from this count. "
            "Raise only if you know many mid-stack steps are unreliable but drift is still worth correcting."
        )
        form.addRow("Max mid-stack fail tolerance:", self.max_fail_spin)
        layout.addLayout(form)

        self.align_masks_check = QCheckBox("Also align matching *_Mask.tif / *_Mask.tiff files")
        self.align_masks_check.setChecked(True)
        layout.addWidget(self.align_masks_check)

        self.expand_canvas_check = QCheckBox("Expand canvas (zero voxel loss at boundaries)")
        self.expand_canvas_check.setChecked(True)
        self.expand_canvas_check.setToolTip(
            "Enlarges the output field of view by the exact maximum cumulative drift so no boundary voxels are clipped."
        )
        layout.addWidget(self.expand_canvas_check)

        output_layout = QHBoxLayout()
        self.output_edit = QLineEdit(str(Path(input_folder) / "Aligned"))
        output_layout.addWidget(self.output_edit)
        browse = QPushButton("Browse...")
        browse.clicked.connect(self.choose_output)
        output_layout.addWidget(browse)
        layout.addWidget(QLabel("Output folder:"))
        layout.addLayout(output_layout)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def choose_output(self):
        folder = QFileDialog.getExistingDirectory(self, "Choose alignment output folder", self.output_edit.text())
        if folder:
            self.output_edit.setText(folder)

    def run_auto_detect(self):
        """Start ChannelDetectionWorker; disable the button while the thread runs."""
        self.auto_detect_btn.setEnabled(False)
        self.auto_detect_btn.setText("Detecting…")
        self.channel_info_label.setText("Scanning channels — this may take a few seconds…")
        self._detection_worker = ChannelDetectionWorker(self.input_folder)
        self._detection_worker.result.connect(self.show_channel_results)
        self._detection_worker.error.connect(self.show_detection_error)
        self._detection_worker.start()

    def show_channel_results(self, results):
        """Fill the spinbox with the best channel and display a per-channel summary."""
        self.auto_detect_btn.setEnabled(True)
        self.auto_detect_btn.setText("Auto-detect…")
        if not results:
            self.channel_info_label.setText("No channels found.")
            return
        best_ch, _ = max(results, key=lambda x: x[1])
        self.reference_channel_spin.setValue(best_ch)
        parts = [
            f"Ch {ch}: {pct:.0f} %{'  ← recommended' if ch == best_ch else ''}"
            for ch, pct in results
        ]
        self.channel_info_label.setText("  |  ".join(parts))

    def show_detection_error(self, text):
        self.auto_detect_btn.setEnabled(True)
        self.auto_detect_btn.setText("Auto-detect…")
        self.channel_info_label.setText(f"Detection failed: {text}")

    def values(self):
        return (
            self.output_edit.text().strip(),
            self.reference_channel_spin.value(),
            self.align_masks_check.isChecked(),
            self.max_step_spin.value(),
            self.max_fail_spin.value(),
            self.expand_canvas_check.isChecked(),
        )



class ReportOptionsDialog(QDialog):

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Generated Output Options")
        self.setMinimumWidth(470)
        self.settings = QSettings("MyLab", "ExQt")

        layout = QVBoxLayout(self)
        note = QLabel(
            "The original *_Output_Batch_*.csv is always saved as the complete machine/audit table.\n"
            "Choose only the additional human-facing reports you need."
        )
        note.setWordWrap(True)
        layout.addWidget(note)

        self.controls = {}
        choices = [
            ("report_excel", "Clean Excel report", "Summary plus primary, excluded, raw and QC-policy sheets."),
            ("report_primary_csv", "Primary QC-valid CSV", "Compact table used for primary analysis."),
            ("report_excluded_csv", "QC-excluded CSV", "Reason-focused table for troubleshooting exclusions."),
            ("report_standard_plots", "Standard descriptive plots", "Volume, intensity and density overview."),
            ("report_mode_a_plots", "Radial FA Profiling plots", "Generated only when Radial FA Profiling is active in 3D."),
            ("report_partitioning_plots", "Partitioning & Classification plots (K_part)", "Size vs K_part and FA vs K_part morphology profiling."),
            ("report_raw_audit_csv", "Extra full raw audit CSV", "Usually unnecessary: duplicates the source table with reporting flags and diameters."),
        ]
        for key, label, tooltip in choices:
            checkbox = QCheckBox(label)
            checkbox.setToolTip(tooltip)
            checkbox.setChecked(_safe_bool(
                self.settings.value(key, REPORT_DEFAULTS[key]), REPORT_DEFAULTS[key]
            ))
            self.controls[key] = checkbox
            layout.addWidget(checkbox)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def accept(self):
        for key, checkbox in self.controls.items():
            self.settings.setValue(key, checkbox.isChecked())
        super().accept()


class SizePreviewDialog(QDialog):
    #Explore the pre-ROI mask-size distribution and publish chosen report bounds.
    range_set = Signal(float, float)

    def __init__(self, size_table, minimum, maximum, calibration_summary="", parent=None):
        super().__init__(parent)
        self.setWindowTitle("Choose analyzed biological-size range")
        self.resize(1050, 650)
        self.values = size_table["biological_size"].to_numpy(dtype=float)
        self.unit = str(size_table["unit"].iloc[0])
        self.file_count = int(size_table["filename"].nunique())
        self._dragged_boundary = None
        self._histogram_axis = None
        self._count_axis = None
        #None shows the complete measured distribution. Pressing Set range stores a focused display interval without changing any source data.
        self._view_bounds = None

        layout = QVBoxLayout(self)
        note = QLabel(
            "Preview of all supplied masks after the raw pixel/voxel noise filter, but before manual ROI. "
            "Drag the blue/red handles or enter exact values. Set range updates the main window and fits both graphs "
        )
        note.setWordWrap(True)
        layout.addWidget(note)
        if calibration_summary:
            calibration_label = QLabel(calibration_summary)
            calibration_label.setWordWrap(True)
            calibration_label.setStyleSheet("font-weight: bold; margin: 4px 0;")
            layout.addWidget(calibration_label)

        controls = QHBoxLayout()
        self.minimum_spin = QDoubleSpinBox()
        self.minimum_spin.setRange(0.0, 100000.0)
        self.minimum_spin.setDecimals(5)
        self.minimum_spin.setValue(float(minimum))
        self.maximum_spin = QDoubleSpinBox()
        self.maximum_spin.setRange(0.00001, 100000.0)
        self.maximum_spin.setDecimals(5)
        self.maximum_spin.setValue(float(maximum))
        controls.addWidget(QLabel(f"Minimum ({self.unit}):"))
        controls.addWidget(self.minimum_spin)
        controls.addWidget(QLabel(f"Maximum ({self.unit}):"))
        controls.addWidget(self.maximum_spin)
        controls.addStretch()
        self.selection_label = QLabel()
        controls.addWidget(self.selection_label)
        layout.addLayout(controls)
        self.set_status_label = QLabel(
            f"Main-window range: {float(minimum):g}-{float(maximum):g} {self.unit}"
        )
        self.set_status_label.setStyleSheet("color: #6c757d;")
        layout.addWidget(self.set_status_label)

        self.figure = Figure(figsize=(10, 5), constrained_layout=True)
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas, stretch=1)

        self.buttons = QDialogButtonBox(QDialogButtonBox.Close)
        self.set_button = self.buttons.addButton("Set range", QDialogButtonBox.ApplyRole)
        self.set_button.setToolTip(
            "Apply the displayed range and fit both graph axes to it while keeping this preview open."
        )
        self.set_button.clicked.connect(self._publish_range)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)

        self.minimum_spin.valueChanged.connect(self._redraw)
        self.maximum_spin.valueChanged.connect(self._redraw)
        self.canvas.mpl_connect("button_press_event", self._graph_pressed)
        self.canvas.mpl_connect("motion_notify_event", self._graph_dragged)
        self.canvas.mpl_connect("button_release_event", self._graph_released)
        self._redraw()

    def selected_bounds(self):
        return self.minimum_spin.value(), self.maximum_spin.value()

    def _publish_range(self):
        minimum, maximum = self.selected_bounds()
        if minimum < maximum:
            self._view_bounds = (minimum, maximum)
            self.range_set.emit(minimum, maximum)
            self.set_status_label.setText(
                f"Range set and graph fitted to {minimum:g}-{maximum:g} {self.unit}."
            )
            self.set_status_label.setStyleSheet("color: #198754; font-weight: bold;")
            self._redraw()

    def _nearest_boundary(self, event, require_handle=False):
        if event.inaxes is None or event.xdata is None:
            return None
        minimum, maximum = self.selected_bounds()
        candidates = {"minimum": minimum}
        if event.inaxes is self._histogram_axis:
            candidates["maximum"] = maximum
        elif event.inaxes is not self._count_axis:
            return None

        distances = {
            name: abs(event.x - event.inaxes.transData.transform((value, 0))[0])
            for name, value in candidates.items()
        }
        nearest = min(distances, key=distances.get)
        if require_handle and distances[nearest] > 18:
            return None
        return nearest

    def _set_boundary(self, boundary, x_value):
        x_value = max(0.0, float(x_value))
        step = 10 ** (-self.minimum_spin.decimals())
        if boundary == "minimum":
            self.minimum_spin.setValue(min(x_value, self.maximum_spin.value() - step))
        elif boundary == "maximum":
            self.maximum_spin.setValue(max(x_value, self.minimum_spin.value() + step))

    def _graph_pressed(self, event):
        if event.button != 1 or event.inaxes is None or event.xdata is None:
            return
        self._dragged_boundary = self._nearest_boundary(event, require_handle=True)
        if self._dragged_boundary is None:
            #Preserve the original quick-click behaviour away from a handle.
            self._set_boundary(self._nearest_boundary(event), event.xdata)

    def _graph_dragged(self, event):
        if self._dragged_boundary is not None and event.inaxes is not None and event.xdata is not None:
            self._set_boundary(self._dragged_boundary, event.xdata)

    def _graph_released(self, _event):
        self._dragged_boundary = None

    def _redraw(self):
        minimum, maximum = self.selected_bounds()
        valid_range = minimum < maximum
        selected = (self.values >= minimum) & (self.values <= maximum) if valid_range else np.zeros_like(self.values, dtype=bool)
        self.set_button.setEnabled(valid_range)
        if valid_range:
            self.selection_label.setText(
                f"Selected: {int(selected.sum())}/{len(self.values)} objects from {self.file_count} file(s)"
            )
            self.selection_label.setStyleSheet("")
        else:
            self.selection_label.setText("Minimum must be smaller than maximum.")
            self.selection_label.setStyleSheet("color: #d9534f; font-weight: bold;")

        self.figure.clear()
        histogram_axis = self.figure.add_subplot(1, 2, 1)
        count_axis = self.figure.add_subplot(1, 2, 2)
        self._histogram_axis = histogram_axis
        self._count_axis = count_axis

        if self._view_bounds is not None:
            view_minimum, view_maximum = self._view_bounds
            view_span = max(view_maximum - view_minimum, 10 ** (-self.minimum_spin.decimals()))
            view_padding = view_span * 0.06
            view_left = max(0.0, view_minimum - view_padding)
            view_right = view_maximum + view_padding
            displayed = selected & (self.values >= view_minimum) & (self.values <= view_maximum)
            displayed_values = self.values[displayed]
            bin_count = int(np.clip(np.sqrt(max(len(displayed_values), 1)) * 2, 6, 30))
            bin_edges = np.linspace(view_minimum, view_maximum, bin_count + 1)
            histogram_axis.hist(
                displayed_values,
                bins=bin_edges,
                color="#2a82da",
                alpha=0.85,
                label="Objects in selected range",
            )
        else:
            view_left = None
            view_right = None
            bin_count = int(np.clip(np.sqrt(max(len(self.values), 1)) * 2, 8, 40))
            bin_edges = np.histogram_bin_edges(self.values, bins=bin_count)
            histogram_axis.hist(self.values, bins=bin_edges, color="#9aa9bc", alpha=0.65, label="All objects")
            if selected.any():
                histogram_axis.hist(
                    self.values[selected], bins=bin_edges, color="#2a82da", alpha=0.85, label="Selected range"
                )
        histogram_axis.axvline(minimum, color="#1f77b4", linestyle="--", linewidth=2, label="Minimum")
        histogram_axis.axvline(maximum, color="#d9534f", linestyle="--", linewidth=2, label="Maximum")
        handle_y = histogram_axis.get_ylim()[1] * 0.96
        histogram_axis.scatter(
            [minimum, maximum], [handle_y, handle_y],
            color=["#1f77b4", "#d9534f"], s=90, zorder=5, edgecolor="white", linewidth=1.2,
        )
        histogram_axis.set_title("A) Mask-size distribution")
        histogram_axis.set_xlabel(f"Biological size ({self.unit})")
        histogram_axis.set_ylabel("Object count")
        histogram_axis.legend(fontsize=8)

        if view_left is not None and view_right is not None:
            histogram_axis.set_xlim(view_left, view_right)

        if self._view_bounds is not None:
            lower, upper = self._view_bounds
        else:
            lower = max(0.0, min(float(self.values.min()), minimum) * 0.9)
            #The candidate-minimum curve only needs the measured data domain. A deliberately very large maximum must not flatten the useful part of the graph.
            upper = max(float(self.values.max()), minimum)
        if upper <= lower:
            upper = lower + 1.0
        thresholds = np.linspace(lower, upper, 250)
        retained = np.array([np.count_nonzero((self.values >= x) & (self.values <= maximum)) for x in thresholds])
        count_axis.plot(thresholds, retained, color="#6f42c1", linewidth=2)
        count_axis.scatter(
            [minimum], [int(selected.sum())], color="#2a82da", s=90,
            zorder=5, edgecolor="white", linewidth=1.2,
        )
        count_axis.axvline(minimum, color="#1f77b4", linestyle="--", linewidth=1.5)
        count_axis.set_title("B) Retained count vs. minimum size")
        count_axis.set_xlabel(f"Candidate minimum ({self.unit}); maximum fixed at {maximum:g}")
        count_axis.set_ylabel("Objects retained")
        count_axis.set_ylim(bottom=0)
        if view_left is not None and view_right is not None:
            count_axis.set_xlim(view_left, view_right)
        self.canvas.draw_idle()


class AdvancedSettingsDialog(QDialog):
    def __init__(self, parent=None, mode="3d"):
        super().__init__(parent)
        self.setWindowTitle("Settings")
        self.setMinimumWidth(460)

        self.settings = QSettings("MyLab", "ExQt")

        layout = QVBoxLayout()
        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)

        self.pixel_size_spin = QDoubleSpinBox()
        self.pixel_size_spin.setRange(1.0, 2000.0)
        self.pixel_size_spin.setSingleStep(0.1)
        self.pixel_size_spin.setDecimals(1)
        form.addRow("Pixel Size XY (nm):", self.pixel_size_spin)

        self.z_step_spin = QDoubleSpinBox()
        self.z_step_spin.setRange(1.0, 5000.0)
        self.z_step_spin.setSingleStep(0.1)
        self.z_step_spin.setDecimals(1)
        form.addRow("Z-step (nm):", self.z_step_spin)
        self.z_step_label = form.labelForField(self.z_step_spin)

        self.signal_spin = QSpinBox()
        form.addRow("Signal Channel:", self.signal_spin)

        self.dapi_spin = QSpinBox()
        form.addRow("DAPI Channel:", self.dapi_spin)

        self.detector_offset_spin = QDoubleSpinBox()
        self.detector_offset_spin.setRange(0.0, 65535.0)
        self.detector_offset_spin.setDecimals(1)
        self.detector_offset_spin.setToolTip(
            "Detector/camera offset in ADU subtracted from condensate and nucleoplasm intensity for K_part. "
            "Take it from the acquisition settings or a dark frame; 0 = no subtraction "
            "(e.g. photon-counting HyD). It is never estimated from the image."
        )
        self.detector_offset_lif_button = QPushButton("From .lif…")
        self.detector_offset_lif_button.setToolTip(
            "Read the detector settings of a Leica .lif file. If every active detector is in "
            "photon-counting mode, the offset is 0 ADU and the file is recorded as its source."
        )
        self.detector_offset_lif_button.clicked.connect(self.read_offset_from_lif)
        offset_row = QWidget()
        offset_layout = QHBoxLayout(offset_row)
        offset_layout.setContentsMargins(0, 0, 0, 0)
        offset_layout.addWidget(self.detector_offset_spin, stretch=1)
        offset_layout.addWidget(self.detector_offset_lif_button)
        form.addRow("Detector offset (ADU):", offset_row)
        self._offset_source = {"method": "manual"}

        self.raw_min_voxels_spin = QSpinBox()
        self.raw_min_voxels_spin.setRange(1, 10000)
        self.raw_min_voxels_spin.setToolTip(
            "Early connected-component noise floor before calibrated biological-size filtering. "
            "This is separate from the analyzed µm²/µm³ range and from Radial FA Profiling layer validity."
        )
        form.addRow("Noise filter (voxels):", self.raw_min_voxels_spin)

        self.mode_a_enabled_check = QCheckBox("Enable Radial FA Profiling (3D)")
        self.mode_a_enabled_check.setToolTip(
            "Adds geometric core-shell Fractional Anisotropy metrics. "
            "It does not directly measure stiffness, liquidity, or viscosity."
        )
        form.addRow("", self.mode_a_enabled_check)

        self.mode_a_min_core_spin = QSpinBox()
        self.mode_a_min_core_spin.setRange(1, 100000)
        self.mode_a_min_core_spin.setToolTip(
            "Minimum voxels required for a layer FA to be valid; the same value also sets the minimum core size."
        )
        form.addRow("Min. voxels per layer:", self.mode_a_min_core_spin)
        self.mode_a_label = form.labelForField(self.mode_a_min_core_spin)

        self.mode_a_exclude_split_check = QCheckBox(
            "Require Z-topology PASS for primary comparison"
        )
        self.mode_a_exclude_split_check.setToolTip(
            "Ignores tiny detached components, grades persistent substantial splitting as PASS/REVIEW/FAIL, "
            "and keeps REVIEW/FAIL objects in the CSV outside the primary comparison."
        )
        form.addRow("", self.mode_a_exclude_split_check)

        self.dark_mode_check = QCheckBox("Dark mode")
        form.addRow("Appearance:", self.dark_mode_check)

        layout.addLayout(form)

        self.load_adv_settings()
        self._apply_mode_state(mode)

        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)

        self.setLayout(layout)

    def _apply_mode_state(self, mode):
        is_3d = mode == "3d"
        self.z_step_spin.setEnabled(is_3d)
        if self.z_step_label is not None:
            self.z_step_label.setEnabled(is_3d)
        if self.mode_a_label is not None:
            self.mode_a_label.setEnabled(is_3d)
        self.mode_a_enabled_check.setEnabled(is_3d)
        self.mode_a_min_core_spin.setEnabled(is_3d)
        self.mode_a_exclude_split_check.setEnabled(is_3d)
        tooltip = "" if is_3d else "Z-step is only used in 3D mode and is ignored for the current Process mode."
        self.z_step_spin.setToolTip(tooltip)
        if not is_3d:
            self.mode_a_enabled_check.setToolTip("Radial FA Profiling is available only in 3D mode.")

    def load_adv_settings(self):
        self.pixel_size_spin.setValue(_safe_float(self.settings.value("adv_pixel_size", DEFAULT_SETTINGS["adv_pixel_size"]), DEFAULT_SETTINGS["adv_pixel_size"]))
        self.z_step_spin.setValue(_safe_float(self.settings.value("adv_z_step", DEFAULT_SETTINGS["adv_z_step"]), DEFAULT_SETTINGS["adv_z_step"]))
        self.signal_spin.setValue(int(self.settings.value("adv_signal_ch", DEFAULT_SETTINGS["adv_signal_ch"])))
        self.dapi_spin.setValue(int(self.settings.value("adv_dapi_ch", DEFAULT_SETTINGS["adv_dapi_ch"])))
        self.raw_min_voxels_spin.setValue(int(self.settings.value("raw_min_voxels", DEFAULT_SETTINGS["raw_min_voxels"])))
        self.detector_offset_spin.setValue(_safe_float(self.settings.value("adv_detector_offset", DEFAULT_SETTINGS["adv_detector_offset"]), DEFAULT_SETTINGS["adv_detector_offset"]))
        self._offset_source = _load_offset_source(self.settings)
        self.dark_mode_check.setChecked(_safe_bool(self.settings.value("dark_mode", False)))
        self.mode_a_enabled_check.setChecked(_safe_bool(self.settings.value("mode_a_enabled", DEFAULT_SETTINGS["mode_a_enabled"])))
        self.mode_a_min_core_spin.setValue(int(self.settings.value("mode_a_min_core_voxels", DEFAULT_SETTINGS["mode_a_min_core_voxels"])))
        self.mode_a_exclude_split_check.setChecked(_safe_bool(self.settings.value("mode_a_exclude_split_slices", DEFAULT_SETTINGS["mode_a_exclude_split_slices"])))

    def accept(self):
        self.settings.setValue("adv_pixel_size", self.pixel_size_spin.value())
        self.settings.setValue("adv_z_step", self.z_step_spin.value())
        self.settings.setValue("adv_signal_ch", self.signal_spin.value())
        self.settings.setValue("adv_dapi_ch", self.dapi_spin.value())
        self.settings.setValue("raw_min_voxels", self.raw_min_voxels_spin.value())
        self.settings.setValue("adv_detector_offset", self.detector_offset_spin.value())
        source = self._offset_source
        if source.get("method") == "lif_photon_counting" and self.detector_offset_spin.value() != 0.0:
            source = {"method": "manual"}  # the user changed the value after reading the LIF file
        self.settings.setValue("adv_detector_offset_source", json.dumps(source))
        self.settings.setValue("dark_mode", self.dark_mode_check.isChecked())
        apply_theme(self.dark_mode_check.isChecked())
        self.settings.setValue("mode_a_enabled", self.mode_a_enabled_check.isChecked())
        self.settings.setValue("mode_a_min_core_voxels", self.mode_a_min_core_spin.value())
        self.settings.setValue("mode_a_exclude_split_slices", self.mode_a_exclude_split_check.isChecked())
        super().accept()

    #Read the detector offset from a Leica .lif header; only photon counting sets a value.
    def read_offset_from_lif(self):
        path, _ = QFileDialog.getOpenFileName(self, "Leica LIF file", "", "Leica LIF (*.lif)")
        if not path:
            return
        try:
            value, source = offset_from_lif(path)
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, "Detector offset from LIF", str(error))
            return
        lines = "\n".join(
            f"  • {d['name']} (channel {d['channel']}): {d['mode'] or 'unknown mode'}, "
            f"Offset {d['offset']}, Gain {d['gain']}"
            for d in source["detectors"]
        )
        details = f"{source['file']} – active detectors:\n{lines or '  (none)'}"
        if value is None:
            QMessageBox.warning(
                self, "Detector offset from LIF",
                f"{details}\n\nNo offset was set: {source['message']}",
            )
            return
        self.detector_offset_spin.setValue(value)
        self._offset_source = source
        QMessageBox.information(
            self, "Detector offset from LIF",
            f"{details}\n\nAll active detectors count photons, so the offset is {value:g} ADU. "
            "The file is recorded as the source of this value when you press OK.",
        )

class ExQt(QMainWindow): 
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ExQt: Analysis of nuclear condensates")
        self.resize(1200, 800) 
        self.create_menu()

        self.left_panel_layout = QVBoxLayout() 

        self.input_label = QLabel("Input folder:")
        self.left_panel_layout.addWidget(self.input_label)

        folder_layout = QHBoxLayout()
        self.folder_input = QLineEdit()
        self.folder_input.setPlaceholderText("Choose folder with TIF images and TIF masks...")
        self.folder_input.setReadOnly(True)
        
        self.btn_browse = QPushButton("Browse")
        self.btn_browse.clicked.connect(self.choose_folder)
        
        folder_layout.addWidget(self.folder_input)
        folder_layout.addWidget(self.btn_browse)
        
        self.left_panel_layout.addLayout(folder_layout)

        self.output_label = QLabel("Output folder:")
        self.left_panel_layout.addWidget(self.output_label)

        self.output_layout = QHBoxLayout()

        self.output_path_edit = QLineEdit()
        self.output_path_edit.setPlaceholderText("Choose folder for .csv output...")
        self.output_layout.addWidget(self.output_path_edit)

        self.btn_output_browse = QPushButton("Browse...")
        self.btn_output_browse.clicked.connect(self.choose_output_folder)
        self.output_layout.addWidget(self.btn_output_browse)

        self.left_panel_layout.addLayout(self.output_layout)

        self.same_folder_checkbox = QCheckBox("Save next to input")
        self.same_folder_checkbox.setToolTip("Write the results into the input folder instead of a separate output folder.")
        self.same_folder_checkbox.toggled.connect(self.switch_same_folder)
        self.left_panel_layout.addWidget(self.same_folder_checkbox)
        
        form_layout = QFormLayout()

        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["3d", "2d", "single_slice"])
        form_layout.addRow("Process mode:", self.mode_combo)

        self.exp_factor_spin = QDoubleSpinBox()
        self.exp_factor_spin.setMinimum(0.1)
        # 3 decimals: a rounded ExF (4.25 -> 4.3) changes volumes by ExF^3 (~3.6 %).
        self.exp_factor_spin.setDecimals(3)
        self.exp_factor_spin.setSingleStep(0.01)
        self.exp_factor_spin.setValue(4.0)
        form_layout.addRow("Expansion factor:", self.exp_factor_spin)

        self.show_napari_check = QCheckBox("Show Napari preview")
        self.show_napari_check.setChecked(True)
        form_layout.addRow("", self.show_napari_check)

        self.auto_roi_check = QCheckBox("Auto-ROI")
        self.auto_roi_check.setToolTip("Auto-ROI treats entire image as a single cell (Per-FOV metrics).")
        form_layout.addRow("", self.auto_roi_check)

        self.review_check = QCheckBox("Review each image")
        self.review_check.setToolTip("Pause after each image to approve or discard its segmentation before the next file.")
        form_layout.addRow("", self.review_check)

        report_row = QHBoxLayout()
        self.generate_reports_check = QCheckBox("Generate")
        self.generate_reports_check.setToolTip(
            "The original machine/audit CSV is always saved. This switch controls only additional reports."
        )
        self.report_options_button = QPushButton("Configure...")
        self.report_options_button.clicked.connect(self.open_report_options)
        report_row.addWidget(self.generate_reports_check)
        report_row.addWidget(self.report_options_button)
        form_layout.addRow("Reports:", report_row)

        self.plot_min_size_spin = QDoubleSpinBox()
        self.plot_min_size_spin.setRange(0.0, 100000.0)
        self.plot_min_size_spin.setDecimals(5)
        self.plot_min_size_spin.setToolTip(
            "Lower calibrated biological-size bound used by primary statistics, clean CSVs and plots. "
            "The original raw audit CSV remains unfiltered."
        )
        form_layout.addRow("Min size:", self.plot_min_size_spin)

        self.plot_max_size_spin = QDoubleSpinBox()
        self.plot_max_size_spin.setRange(0.0001, 100000.0)
        self.plot_max_size_spin.setDecimals(5)
        self.plot_max_size_spin.setToolTip(
            "Upper calibrated biological-size bound used by primary statistics, clean CSVs and plots. "
            "The original raw audit CSV remains unfiltered."
        )
        form_layout.addRow("Max size:", self.plot_max_size_spin)

        self.size_preview_button = QPushButton("Preview sizes...")
        self.size_preview_button.setToolTip(
            "Scan the supplied masks before analysis and choose the biological-size range interactively. "
            "The preview is calculated before manual ROI."
        )
        self.size_preview_button.clicked.connect(self.open_size_preview)
        form_layout.addRow("Size range:", self.size_preview_button)

        self.left_panel_layout.addLayout(form_layout)
        self.left_panel_layout.addStretch()
        self.btn_run = QPushButton("Start analysis")
        self.btn_run.setStyleSheet("background-color: #2a82da; color: white; padding: 10px; font-weight: bold;")
        self.left_panel_layout.addWidget(self.btn_run)

        self.btn_confirm_roi = QPushButton("Confirm ROI and Continue")
        self.btn_confirm_roi.setStyleSheet("background-color: #28a745; color: white; padding: 10px; font-weight: bold;")
        self.btn_confirm_roi.hide()
        self.btn_confirm_roi.clicked.connect(self.confirm_roi)
        self.left_panel_layout.addWidget(self.btn_confirm_roi)

        self.btn_next_image = QPushButton("Approve & Next Image")
        self.btn_next_image.setStyleSheet("background-color: #ff9800; color: white; padding: 10px; font-weight: bold;")
        self.btn_next_image.hide()
        self.btn_next_image.clicked.connect(self.next_image_confirmed)
        self.left_panel_layout.addWidget(self.btn_next_image)

        self.btn_stop_review = QPushButton("Stop and Discard This Image")
        self.btn_stop_review.setStyleSheet("background-color: #d9534f; color: white; padding: 10px; font-weight: bold;")
        self.btn_stop_review.setToolTip("Stops the batch here and discards the image currently shown. Previously approved images are still saved.")
        self.btn_stop_review.hide()
        self.btn_stop_review.clicked.connect(self.stop_and_discard)
        self.left_panel_layout.addWidget(self.btn_stop_review)

        self.align_progress_bar = QProgressBar()
        self.align_progress_bar.setRange(0, 100)
        self.align_progress_bar.setValue(0)
        self.align_progress_bar.setTextVisible(True)
        self.align_progress_bar.setStyleSheet("QProgressBar { text-align: center; font-weight: bold; }")
        self.align_progress_bar.hide()
        self.left_panel_layout.addWidget(self.align_progress_bar)

        self.status_label = QLabel("Prepared")
        self.status_label.setStyleSheet("color: gray; font-weight: bold; margin-top: 10px;")
        # Long messages (e.g. the ROI drawing help) must wrap, otherwise the label forces the
        # whole left panel to the width of the text and squeezes the Napari viewer.
        self.status_label.setWordWrap(True)
        self.left_panel_layout.addWidget(self.status_label)

        self.viewer = napari.Viewer(show=False)
        left_panel = QWidget()
        left_panel.setLayout(self.left_panel_layout)

        # Adjustable split between the controls and the viewer; the position is remembered.
        self.main_splitter = QSplitter(Qt.Horizontal)
        self.main_splitter.addWidget(left_panel)
        self.main_splitter.addWidget(self.viewer.window._qt_window)
        self.main_splitter.setChildrenCollapsible(False)
        self.main_splitter.setStretchFactor(0, 1)
        self.main_splitter.setStretchFactor(1, 3)
        self.main_splitter.setSizes([340, 1460])
        self.setCentralWidget(self.main_splitter)

        self.settings = QSettings("MyLab", "ExQt")
        #TIFF calibration is applied only when the user explicitly chooses it (calibration_source).
        self.detected_metadata_by_file = {}
        self.calibration_source = CALIBRATION_SOURCE_GUI
        self.calibration_warning = ""
        self.metadata_scanned_folder = None
        self.load_settings()

        self.btn_run.clicked.connect(self.start_analysis)
        self.mode_combo.currentTextChanged.connect(self.on_mode_changed)
        self.on_mode_changed(self.mode_combo.currentText()) 

    def start_alignment(self):
        folder_path = self.folder_input.text().strip()
        if not folder_path or not Path(folder_path).is_dir():
            QMessageBox.warning(self, "Alignment", "Please select a valid input folder first.")
            return

        source_files = _source_tiff_files(folder_path)
        if not source_files:
            QMessageBox.warning(self, "Alignment", "No source TIFF stacks were found in the selected folder.")
            return

        dialog = AlignmentOptionsDialog(folder_path, self)
        if dialog.exec() != QDialog.Accepted:
            return
        output_text, reference_channel, align_masks, max_step_px, max_fail_pct, expand_canvas = dialog.values()
        if not output_text:
            QMessageBox.warning(self, "Alignment settings", "Please choose an output folder.")
            return
        output_folder = Path(output_text).expanduser().resolve()
        if output_folder == Path(folder_path).resolve():
            QMessageBox.warning(
                self,
                "Alignment settings",
                "The alignment output folder must be different from the input folder. Source files are protected from overwrite.",
            )
            return

        alignment_config = AlignmentConfig(
            max_step_px=float(max_step_px),
            max_fail_fraction=max_fail_pct / 100.0,
            expand_canvas=expand_canvas,
        )
        self.align_stacks_action.setEnabled(False)
        self.btn_run.setEnabled(False)
        self.align_progress_bar.setValue(0)
        self.align_progress_bar.show()
        self.status_label.setText("Preparing stack alignment...")
        self.alignment_worker = AlignmentWorker(
            folder_path,
            str(output_folder),
            reference_channel,
            align_masks,
            config=alignment_config,
        )
        self.alignment_worker.progress.connect(self.update_alignment_status)
        self.alignment_worker.completed.connect(self.alignment_completed)
        self.alignment_worker.failed.connect(self.alignment_failed)
        self.alignment_worker.start()

    def update_alignment_status(self, pct, text):
        print(f"ALIGNMENT LOG [{pct}%]: {text}")
        self.align_progress_bar.setValue(pct)
        self.status_label.setText(text)

    def alignment_completed(self, batch_result):
        self.align_stacks_action.setEnabled(True)
        self.btn_run.setEnabled(True)
        self.align_progress_bar.setValue(100)
        self.align_progress_bar.hide()
        records = batch_result.records
        status_counts = records["status"].value_counts().to_dict()
        pass_count = int(status_counts.get("PASS", 0))
        review_count = int(status_counts.get("REVIEW", 0))
        fail_count = int(status_counts.get("FAIL", 0))
        error_count = int(status_counts.get("ERROR", 0))
        exported_count = pass_count + review_count
        self.status_label.setText(
            f"Alignment finished: {pass_count} passed, {review_count} review, "
            f"{fail_count} failed, {error_count} errors."
        )
        message = (
            f"Stacks processed: {len(records)}\n"
            f"Passed (fully reliable): {pass_count}\n"
            f"Review (exported — inspect drift CSV): {review_count}\n"
            f"Failed (not exported): {fail_count}\n"
            f"Errors: {error_count}\n\n"
            f"Batch summary: {batch_result.summary_csv}\n"
            f"The output folder contains aligned TIFFs and drift CSV/PNG files."
        )
        if review_count or fail_count or error_count:
            QMessageBox.warning(self, "Alignment finished with warnings", message)
        else:
            QMessageBox.information(self, "Alignment complete", message)

        if exported_count:
            aligned_folder = str(batch_result.summary_csv.parent)
            switch = QMessageBox.question(
                self,
                "Use aligned data for ExQt?",
                f"Alignment produced {exported_count} usable stack(s).\n\n"
                f"Switch the ExQt input folder to:\n{aligned_folder}\n\n"
                "This prevents accidentally running statistics on the original, unaligned TIFFs.",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.Yes,
            )
            if switch == QMessageBox.Yes:
                self.folder_input.setText(aligned_folder)
                self.metadata_scanned_folder = ""
                self._try_auto_fill_metadata(aligned_folder)
                self.status_label.setText(f"Using aligned data: {aligned_folder}")

    def alignment_failed(self, text):
        self.align_stacks_action.setEnabled(True)
        self.btn_run.setEnabled(True)
        self.align_progress_bar.hide()
        self.status_label.setText(f"Alignment failed: {text}")
        QMessageBox.warning(self, "Alignment failed", text)

    def on_mode_changed(self, mode):
        # Size bounds are volumes in 3D and areas in 2D / single slice.
        unit = " µm³" if mode == "3d" else " µm²"
        self.plot_min_size_spin.setSuffix(unit)
        self.plot_max_size_spin.setSuffix(unit)

    def choose_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Choose folder with data")
        if folder:
            self.folder_input.setText(folder)
            self._try_auto_fill_metadata(folder)

    #Inspect every source TIFF, but never silently overwrite the user's calibration.
    def _try_auto_fill_metadata(self, folder_path):
        folder = Path(folder_path)
        self.metadata_scanned_folder = str(folder.resolve())
        files = _source_tiff_files(folder)
        self.detected_metadata_by_file = {}
        self.calibration_source = CALIBRATION_SOURCE_GUI
        self.calibration_warning = ""
        implausible = []

        for tif_path in files:
            meta = get_metadata_from_tif(tif_path)
            if meta and "pixel_size" in meta and "z_step" in meta:
                try:
                    validate_calibration(meta["pixel_size"], meta["z_step"])
                except ValueError as error:
                    implausible.append(f"{tif_path.name}: {error}")
                    continue
                self.detected_metadata_by_file[tif_path.name] = {
                    "pixel_size_nm": float(meta["pixel_size"]),
                    "z_step_nm": float(meta["z_step"]),
                    "sources": meta.get("sources", {}),
                }

        if implausible:
            QMessageBox.warning(
                self,
                "Implausible TIFF calibration ignored",
                "These files report a calibration outside the supported range and it will NOT be used:\n\n"
                + "\n".join(implausible)
                + "\n\nThe calibration from Settings applies to them.",
            )

        gui_pixel_size_nm = _safe_float(self.settings.value("adv_pixel_size"), DEFAULT_SETTINGS["adv_pixel_size"])
        gui_z_step_nm = _safe_float(self.settings.value("adv_z_step"), DEFAULT_SETTINGS["adv_z_step"])
        expansion_factor = self.exp_factor_spin.value()
        summary = summarize_calibrations(self.detected_metadata_by_file)
        pre_expansion_note = ""
        if any(
            looks_pre_expansion(xy, z, gui_pixel_size_nm, gui_z_step_nm, expansion_factor)
            for xy, z in summary["calibrations"]
        ):
            pre_expansion_note = (
                "\n\nWARNING: the TIFF calibration ≈ current Settings / ExF "
                f"(XY {gui_pixel_size_nm:g} / {expansion_factor:g} nm). The files seem to be calibrated to the "
                "PRE-EXPANSION size; using it would divide by ExF twice."
            )
        if summary["has_mismatch"]:
            calib_details = "\n".join(f"  • XY={xy:g} nm, Z={z:g} nm" for xy, z in summary["calibrations"])
            reply = QMessageBox.question(
                self,
                "Different TIFF calibrations",
                (
                    f"Input files contain {summary['calibration_count']} different physical calibrations:\n\n"
                    f"{calib_details}\n\n"
                    "Yes: apply each file's own TIFF calibration (Settings only for files without it).\n"
                    f"No: apply Settings (XY={gui_pixel_size_nm:g} nm, Z={gui_z_step_nm:g} nm) to all files."
                    f"{pre_expansion_note}"
                ),
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No if pre_expansion_note else QMessageBox.Yes,
            )
            if reply == QMessageBox.Yes:
                self.calibration_source = CALIBRATION_SOURCE_TIFF
                self.status_label.setText(f"Per-file TIFF calibration selected ({summary['calibration_count']} calibrations).")
            else:
                self.status_label.setText("TIFF calibrations ignored; Settings apply to all files.")
        elif summary["calibration_count"] == 1:

            pixel_size_nm, z_step_nm = summary["calibrations"][0]
            reply = QMessageBox.question(
                self,
                "Use detected TIFF calibration?",
                (
                    "All readable input files report the same calibration:\n\n"
                    f"XY pixel size: {pixel_size_nm:g} nm\n"
                    f"Z-step: {z_step_nm:g} nm\n\n"
                    "Use these values in ExQt Settings? They will still be shown "
                    "for confirmation before analysis starts."
                    f"{pre_expansion_note}"
                ),
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No if pre_expansion_note else QMessageBox.Yes,
            )
            # Either answer leaves calibration_source="gui": Settings are what gets applied.
            if reply == QMessageBox.Yes:
                self.settings.setValue("adv_pixel_size", pixel_size_nm)
                self.settings.setValue("adv_z_step", z_step_nm)
                self.status_label.setText(
                    f"TIFF calibration selected: XY={pixel_size_nm:g} nm, Z={z_step_nm:g} nm."
                )
            else:
                self.status_label.setText(
                    f"Detected XY={pixel_size_nm:g} nm, Z={z_step_nm:g} nm; current manual settings retained."
                )
        elif files:
            self.status_label.setText(
                "No usable physical metadata detected; enter calibration manually in Settings."
            )
            
    def create_menu(self):
        menu_bar = self.menuBar()

        tools_menu = menu_bar.addMenu("Tools")
        self.align_stacks_action = QAction("Align Z-stacks...", self)
        self.align_stacks_action.setToolTip(
            "Create aligned copies of source TIFF stacks, optional shifted masks, and drift QC files in a new folder."
        )
        tools_menu.addAction(self.align_stacks_action)
        self.align_stacks_action.triggered.connect(self.start_alignment)

        self.merge_runs_action = QAction("Merge existing runs...", self)
        self.merge_runs_action.setToolTip(
            "Recursively merge compatible ExQt batch CSVs from one folder tree."
        )
        tools_menu.addAction(self.merge_runs_action)
        self.merge_runs_action.triggered.connect(self.merge_existing_runs)

        # A top-level action (no sub-menu): one click opens the settings window.
        self.settings_action = QAction("Settings", self)
        menu_bar.addAction(self.settings_action)
        self.settings_action.triggered.connect(self.open_settings)

    def open_settings(self):
        dialog = AdvancedSettingsDialog(self, mode=self.mode_combo.currentText())
        dialog.exec()

    def open_report_options(self):
        ReportOptionsDialog(self).exec()

    def merge_existing_runs(self):
        root = QFileDialog.getExistingDirectory(
            self,
            "Choose folder containing ExQt run CSVs",
            self.output_path_edit.text() or self.folder_input.text(),
        )
        if not root:
            return
        output_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save merged statistics",
            str(Path(root) / "Merged_Stats.xlsx"),
            "Excel workbook (*.xlsx)",
        )
        if not output_path:
            return
        if not output_path.lower().endswith(".xlsx"):
            output_path += ".xlsx"
        include_raw = QMessageBox.question(
            self,
            "Include merged raw audit?",
            (
                "Include every raw connected component in an additional workbook sheet?\n\n"
                "Recommended: No. Primary and QC-excluded objects plus per-run summaries are normally sufficient."
            ),
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No,
        ) == QMessageBox.Yes
        try:
            result = merge_statistics_folder(root, output_path, include_raw=include_raw)
        except QCPolicyMismatchError as exc:
            QMessageBox.critical(
                self,
                "Runs cannot be merged",
                str(exc),
            )
            return
        except Exception as exc:
            QMessageBox.critical(self, "Merge failed", str(exc))
            return
        QMessageBox.information(
            self,
            "Merge completed",
            (
                f"Merged {result['run_count']} compatible runs.\n"
                f"Primary objects: {result['primary_count']}\n"
                f"QC-excluded objects: {result['excluded_count']}\n\n"
                f"Workbook:\n{result['excel']}\n\n"
                f"Merged graph:\n{result['plot']}"
            ),
        )

    def choose_output_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Choose output folder")
        if folder:
            self.output_path_edit.setText(folder)

    def open_size_preview(self):
        """Measure saved masks and let the user choose the analyzed size range."""
        folder_path = self.folder_input.text().strip()
        if not folder_path or not os.path.isdir(folder_path):
            QMessageBox.warning(self, "Size preview", "Please select a valid input folder first.")
            return

        folder = Path(folder_path)
        source_files = _source_tiff_files(folder)
        if not source_files:
            QMessageBox.warning(self, "Size preview", "No source TIF files were found in the selected folder.")
            return
        missing_masks = [
            path.name for path in source_files
            if not _matching_mask_path(path).exists()
        ]
        if missing_masks:
            QMessageBox.warning(
                self,
                "Size preview",
                "Missing matching masks for:\n" + "\n".join(missing_masks),
            )
            return

        if self.metadata_scanned_folder != str(folder.resolve()):
            self._try_auto_fill_metadata(folder_path)
        if self.calibration_warning:
            QMessageBox.warning(self, "Calibration mismatch", self.calibration_warning)
            return

        pixel_size_nm = _safe_float(
            self.settings.value("adv_pixel_size", DEFAULT_SETTINGS["adv_pixel_size"]),
            DEFAULT_SETTINGS["adv_pixel_size"],
        )
        z_step_nm = _safe_float(
            self.settings.value("adv_z_step", DEFAULT_SETTINGS["adv_z_step"]),
            DEFAULT_SETTINGS["adv_z_step"],
        )
        expansion_factor = self.exp_factor_spin.value()
        mode = self.mode_combo.currentText()
        min_voxels = int(
            self.settings.value("raw_min_voxels", DEFAULT_SETTINGS["raw_min_voxels"])
        )
        signal_channel = int(
            self.settings.value("adv_signal_ch", DEFAULT_SETTINGS["adv_signal_ch"])
        )

        self.size_preview_button.setEnabled(False)
        self.status_label.setText("Measuring mask-size distribution...")
        QApplication.processEvents()
        try:
            size_table = collect_size_preview(
                input_folder=folder_path,
                mode=mode,
                expansion_factor=expansion_factor,
                pixel_size_nm=pixel_size_nm,
                z_step_nm=z_step_nm,
                min_voxels=min_voxels,
                signal_channel=signal_channel,
                calibration_source=self.calibration_source,
                detected_metadata_by_file=self.detected_metadata_by_file,
            )
        except Exception as error:
            self.status_label.setText("Size preview failed.")
            QMessageBox.warning(self, "Size preview failed", str(error))
            return
        finally:
            self.size_preview_button.setEnabled(True)

        if size_table.empty:
            self.status_label.setText("No objects remain after the raw noise filter.")
            QMessageBox.information(
                self,
                "Size preview",
                "No connected components remain after the raw pixel/voxel noise filter.",
            )
            return

        unit = str(size_table["unit"].iloc[0])
        if mode == "3d":
            calibration_summary = (
                f"Preview calibration: XY={pixel_size_nm:g} nm, Z={z_step_nm:g} nm, "
                f"ExF={expansion_factor:g}x; raw floor={min_voxels} voxels; output={unit}."
            )
        else:
            calibration_summary = (
                f"Preview calibration: XY={pixel_size_nm:g} nm, ExF={expansion_factor:g}x; "
                f"raw floor={min_voxels} pixels; output={unit}."
            )
        dialog = SizePreviewDialog(
            size_table,
            self.plot_min_size_spin.value(),
            self.plot_max_size_spin.value(),
            calibration_summary=calibration_summary,
            parent=self,
        )
        range_was_set = {"value": False}

        def set_preview_range(minimum, maximum):
            """Persist a preview range without forcing the graph window to close."""
            self.plot_min_size_spin.setValue(minimum)
            self.plot_max_size_spin.setValue(maximum)
            self.settings.setValue("plot_min_size", minimum)
            self.settings.setValue("plot_max_size", maximum)
            selected_count = int(
                ((size_table["biological_size"] >= minimum)
                & (size_table["biological_size"] <= maximum)).sum()
            )
            range_was_set["value"] = True
            self.status_label.setText(
                f"Size range set: {minimum:g}-{maximum:g} {unit} "
                f"({selected_count}/{len(size_table)} pre-ROI objects)."
            )

        dialog.range_set.connect(set_preview_range)
        dialog.exec()
        if not range_was_set["value"]:
            self.status_label.setText("Size preview closed without changing the range.")

    def switch_same_folder(self, checked):
        self.output_path_edit.setDisabled(checked)
        self.btn_output_browse.setDisabled(checked)
        
        if checked:
            self.output_path_edit.setText(self.folder_input.text())
        else:
            self.output_path_edit.clear()

    def load_settings(self):
        dark_mode = self.settings.value("dark_mode")
        if dark_mode is not None:
            apply_theme(_safe_bool(dark_mode))
        splitter_state = self.settings.value("main_splitter_state")
        if splitter_state is not None:
            self.main_splitter.restoreState(splitter_state)
        self.folder_input.setText(self.settings.value("input_folder", ""))
        self.output_path_edit.setText(self.settings.value("output_folder", ""))
        
        self.plot_min_size_spin.setValue(_safe_float(self.settings.value("plot_min_size", DEFAULT_SETTINGS["plot_min_size"]), DEFAULT_SETTINGS["plot_min_size"]))
        self.plot_max_size_spin.setValue(_safe_float(self.settings.value("plot_max_size", DEFAULT_SETTINGS["plot_max_size"]), DEFAULT_SETTINGS["plot_max_size"]))
        self.generate_reports_check.setChecked(_safe_bool(self.settings.value("generate_reports", False)))
        
        same_folder_saved = self.settings.value("same_folder", "false")
        if str(same_folder_saved).lower() == "true":
            self.same_folder_checkbox.setChecked(True)

    def closeEvent(self, event):
        self.settings.setValue("main_splitter_state", self.main_splitter.saveState())
        self.settings.setValue("input_folder", self.folder_input.text())
        self.settings.setValue("output_folder", self.output_path_edit.text())
        self.settings.setValue("same_folder", self.same_folder_checkbox.isChecked())
        self.settings.setValue("plot_min_size", self.plot_min_size_spin.value())
        self.settings.setValue("plot_max_size", self.plot_max_size_spin.value())
        self.settings.setValue("generate_reports", self.generate_reports_check.isChecked())
        event.accept()

    def start_analysis(self):
        folder_path = self.folder_input.text()
        if not folder_path or not os.path.isdir(folder_path):
            QMessageBox.warning(self, "Error", "Please select a valid input folder first.")
            return

        folder = Path(folder_path)
        files = _source_tiff_files(folder)
        if not files:
            QMessageBox.warning(self, "Error", "In the selected folder, no source TIF files were found.")
            return
        missing_masks = [f.name for f in files if not _matching_mask_path(f).exists()]

        if missing_masks:
            QMessageBox.warning(self, "Missing Masks", f"For these files there are no corresponding mask files:\n{', '.join(missing_masks)}")
            return 

        if self.plot_min_size_spin.value() >= self.plot_max_size_spin.value():
            QMessageBox.warning(
                self,
                "Invalid Size Range",
                "Analyzed Size Min must be smaller than Analyzed Size Max.",
            )
            return

        if self.metadata_scanned_folder != str(folder.resolve()):
            self._try_auto_fill_metadata(folder_path)

        if self.calibration_warning:
            QMessageBox.warning(self, "Calibration mismatch", self.calibration_warning)
            return

        pixel_size_nm = _safe_float(
            self.settings.value("adv_pixel_size", DEFAULT_SETTINGS["adv_pixel_size"]),
            DEFAULT_SETTINGS["adv_pixel_size"],
        )
        z_step_nm = _safe_float(
            self.settings.value("adv_z_step", DEFAULT_SETTINGS["adv_z_step"]),
            DEFAULT_SETTINGS["adv_z_step"],
        )
        expansion_factor = self.exp_factor_spin.value()
        try:
            validate_calibration(pixel_size_nm, z_step_nm)
        except ValueError as error:
            QMessageBox.warning(self, "Invalid calibration", f"Settings calibration is not usable:\n\n{error}")
            return
        summary = summarize_calibrations(self.detected_metadata_by_file)

        if self.calibration_source == CALIBRATION_SOURCE_TIFF:
            calibration_text = (
                f"Per-file TIFF calibration is ACTIVE ({summary['calibration_count']} distinct calibrations found):\n\n"
            )
            for xy, z in summary["calibrations"]:
                calibration_text += (
                    f"  • XY={xy:g} nm, Z={z:g} nm "
                    f"(biological: Y/X={xy / expansion_factor:g} nm, Z={z / expansion_factor:g} nm)\n"
                )
            without_tiff = [f.name for f in files if f.name not in self.detected_metadata_by_file]
            calibration_text += (
                f"\nFallback calibration (Settings): XY={pixel_size_nm:g} nm, Z={z_step_nm:g} nm, "
                f"used for {len(without_tiff)} file(s) without TIFF calibration\n"
                f"Expansion factor: {expansion_factor:g}×"
            )
            if any(
                looks_pre_expansion(xy, z, pixel_size_nm, z_step_nm, expansion_factor)
                for xy, z in summary["calibrations"]
            ):
                calibration_text += (
                    "\n\nWARNING: a TIFF calibration ≈ Settings / ExF. The files seem to be "
                    "calibrated to the PRE-EXPANSION size; applying ExF again would make volumes "
                    f"{expansion_factor:g}³× too small."
                )
        elif self.mode_combo.currentText() == "3d":
            calibration_text = (
                f"Acquisition XY pixel size: {pixel_size_nm:g} nm\n"
                f"Acquisition Z-step: {z_step_nm:g} nm\n"
                f"Expansion factor: {expansion_factor:g}×\n\n"
                f"Effective biological sampling: Z={z_step_nm / expansion_factor:g} nm, "
                f"Y/X={pixel_size_nm / expansion_factor:g} nm"
            )
        else:
            calibration_text = (
                f"Acquisition XY pixel size: {pixel_size_nm:g} nm\n"
                f"Expansion factor: {expansion_factor:g}×\n\n"
                f"Effective biological XY sampling: {pixel_size_nm / expansion_factor:g} nm"
            )
        detector_offset_adu = _safe_float(
            self.settings.value("adv_detector_offset", DEFAULT_SETTINGS["adv_detector_offset"]),
            DEFAULT_SETTINGS["adv_detector_offset"],
        )
        calibration_text += (
            f"\n\nDetector offset subtracted for K_part: {detector_offset_adu:g} ADU"
            + (" (none)" if detector_offset_adu == 0 else "")
        )
        if self.auto_roi_check.isChecked():
            calibration_text += "\nAuto-ROI: K_part is not computed (no nucleus outline)."
        if _safe_bool(self.settings.value("mode_a_enabled", DEFAULT_SETTINGS["mode_a_enabled"])):
            calibration_text += (
                "\n\nRadial FA Profiling min. voxels per layer: "
                f"{int(self.settings.value('mode_a_min_core_voxels', DEFAULT_SETTINGS['mode_a_min_core_voxels']))}"
            )
        reply = QMessageBox.question(
            self,
            "Confirm calibration",
            calibration_text + "\n\nContinue with these values?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.Yes,
        )
        if reply != QMessageBox.Yes:
            self.status_label.setText("Analysis cancelled: calibration was not confirmed.")
            return


        #Snapshot all GUI values before starting the thread so worker execution is independent of subsequent widget changes
        params = {
            "input_folder": folder_path,
            "output_folder": self.output_path_edit.text(),
            "mode": self.mode_combo.currentText(),
            "expansion_factor": expansion_factor,
            "min_voxels": int(self.settings.value("raw_min_voxels", DEFAULT_SETTINGS["raw_min_voxels"])),
            "auto_roi": self.auto_roi_check.isChecked(),
            "review_each_image": self.review_check.isChecked(),
            "show_napari": self.show_napari_check.isChecked(),
            "generate_reports": self.generate_reports_check.isChecked(),
            **{key: _safe_bool(self.settings.value(key, default), default) for key, default in REPORT_DEFAULTS.items()},
            "plot_min_size": self.plot_min_size_spin.value(),
            "plot_max_size": self.plot_max_size_spin.value(),
            "pixel_size_nm": pixel_size_nm,
            "z_step_nm": z_step_nm,
            "calibration_source": self.calibration_source,
            "calibration_confirmation": "confirmed_by_user_at_run_start",
            "detected_metadata_by_file": self.detected_metadata_by_file,
            "signal_channel": int(self.settings.value("adv_signal_ch", DEFAULT_SETTINGS["adv_signal_ch"])),
            "dapi_channel": int(self.settings.value("adv_dapi_ch", DEFAULT_SETTINGS["adv_dapi_ch"])),
            "detector_offset_adu": detector_offset_adu,
            "detector_offset_source": _load_offset_source(self.settings),
            "mode_a_enabled": _safe_bool(self.settings.value("mode_a_enabled", DEFAULT_SETTINGS["mode_a_enabled"])),
            "mode_a_min_core_voxels": int(self.settings.value("mode_a_min_core_voxels", DEFAULT_SETTINGS["mode_a_min_core_voxels"])),
            "mode_a_exclude_split_slices": _safe_bool(self.settings.value("mode_a_exclude_split_slices", DEFAULT_SETTINGS["mode_a_exclude_split_slices"])),
            **{key: DEFAULT_SETTINGS[key] for key in MODE_A_KEYS if key.startswith("mode_a_z_split_")},
        }

        self.btn_run.setEnabled(False)
        self.btn_run.setText("Processing...")
        self.align_stacks_action.setEnabled(False)
        self.viewer.layers.clear()

        self.worker = AnalysisWorker(params)
        self.worker.layer_ready.connect(self.receive_layer)
        self.worker.progress.connect(self.update_button_text)
        self.worker.request_roi_signal.connect(self.prepare_manual_roi)
        self.worker.request_review_signal.connect(self.prepare_review)
        self.worker.start()
        
    def update_button_text(self, text):
        print(f"GUI LOG: {text}") 
        
        if text == "Done":
            self.btn_stop_review.hide()
            self.btn_run.setEnabled(True)
            self.btn_run.setText("Start analysis")
            self.align_stacks_action.setEnabled(True)
            
            if getattr(getattr(self, "worker", None), "abort_requested", False):
                self.status_label.setText("Stopped by user; the current image was discarded, earlier approved images were saved.")
            elif not self.status_label.text().startswith("Error") and not self.status_label.text().startswith("No data"):
                self.status_label.setText("Analysis completed successfully.")
                QMessageBox.information(self, "Done", "Analysis was completed successfully\n\nResults are saved in CSV.")
        elif text.startswith("Error") or text.startswith("No data"):
            self.status_label.setText(text)
            QMessageBox.warning(self, "Analysis Failed", text) 
        else:
            self.status_label.setText(text)
            self.btn_run.setText("Processing...")

    def receive_layer(self, layer_info):
        layer_type = layer_info.get("type", "image")

        if layer_type == "clear_layers":
            self.viewer.layers.clear()
            return

        name = layer_info.get("name", "Unknown Layer")
        data = layer_info.get("data")
        kwargs = layer_info.get("kwargs", {})

        if layer_type == "image":
            self.viewer.add_image(data, name=name, **kwargs)
        elif layer_type == "labels":
            self.viewer.add_labels(data, name=name, **kwargs)
        elif layer_type == "points":
            self.viewer.add_points(data, name=name, **kwargs)

    def prepare_manual_roi(self, info):
        self._current_roi_info = info
        shape = info["shape"]
        ndim = len(shape)
        first_layer = next(iter(self.viewer.layers), None)
        layer_scale = getattr(first_layer, "scale", None) if first_layer is not None else None

        for name in ("Draw ROI", "Paint ROI"):
            if name in self.viewer.layers:
                self.viewer.layers.remove(name)

        kwargs = {
            "name": "Draw ROI",
            "ndim": ndim,
            "opacity": 0.5,
            "face_color": "cyan",
            "edge_color": "blue",
            "edge_width": 2,
        }
        if layer_scale is not None:
            kwargs["scale"] = layer_scale

        shapes_layer = self.viewer.add_shapes(**kwargs)
        shapes_layer.mode = "add_ellipse"

        self.btn_run.hide()
        self.btn_confirm_roi.show()
        self.btn_stop_review.show()
        if hasattr(self, "status_label"):
            self.status_label.setText("Draw one shape per nucleus on 1 or more Z-slices (overlapping shapes on different slices = same nucleus, interpolated; outside its slices = 0; one slice = whole Z). Shapes must not touch on the same slice. Click 'Confirm ROI'...")

    def confirm_roi(self):
        info = getattr(self, "_current_roi_info", {})
        shape = info.get("shape", None)

        mask_data = None
        roi_layer_name = None
        if "Draw ROI" in self.viewer.layers:
            roi_layer_name = "Draw ROI"
            layer = self.viewer.layers["Draw ROI"]
            if hasattr(layer, "to_labels") and shape is not None:
                mask_data = layer.to_labels(labels_shape=shape)
            elif hasattr(layer, "to_masks") and shape is not None:
                masks = layer.to_masks(mask_shape=shape)
                mask_data = np.sum(masks, axis=0, dtype=int) if masks.ndim > len(shape) else masks.astype(int)
        elif "Paint ROI" in self.viewer.layers:
            roi_layer_name = "Paint ROI"
            mask_data = self.viewer.layers["Paint ROI"].data

        # Check the multi-nucleus grouping here so a drawing mistake can be fixed
        # instead of failing the whole image later in the worker.
        if mask_data is not None and info.get("is_3d") and mask_data.ndim == 3 and mask_data.max() > 0:
            try:
                _interpolate_or_extrude_roi(mask_data)
            except ValueError as error:
                QMessageBox.warning(self, "ROI cannot be used", f"{error}\n\nFix the drawn shapes and confirm again.")
                return
        if roi_layer_name is not None:
            self.viewer.layers.remove(roi_layer_name)

        roi_source = "manual"
        if mask_data is None or mask_data.max() == 0:
            # Never fall back to the whole FOV silently: ask, and record it in the CSV.
            reply = QMessageBox.question(
                self,
                "No ROI drawn",
                "No ROI was drawn for this image.\n\n"
                "Yes: analyse the whole field of view (recorded as roi_source = empty_fallback_fov).\n"
                "No: go back and draw an ROI.",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            if reply != QMessageBox.Yes:
                self.prepare_manual_roi(info)
                return
            roi_source = "empty_fallback_fov"
            first_layer = next(iter(self.viewer.layers), None)
            if first_layer is not None:
                mask_data = np.ones(first_layer.data.shape, dtype=int)
            elif shape is not None:
                mask_data = np.ones(shape, dtype=int)
            else:
                mask_data = np.ones((1, 1), dtype=int)

        self.worker.user_roi_data = mask_data
        self.worker.user_roi_source = roi_source
        self.btn_confirm_roi.hide()
        self.btn_stop_review.hide()
        self.btn_run.show()
        if hasattr(self, "status_label"):
            self.status_label.setText("ROI confirmed. Processing...")
        self.worker.roi_event.set()

    def prepare_review(self):
        self.btn_run.hide()
        self.btn_confirm_roi.hide()
        self.btn_next_image.show()
        self.btn_stop_review.show()
        if hasattr(self, "status_label"):
            self.status_label.setText("Waiting for user review...")

    def next_image_confirmed(self):
        self.btn_next_image.hide()
        self.btn_stop_review.hide()
        self.btn_run.show()
        if hasattr(self, "status_label"):
            self.status_label.setText("Processing next...")
        self.worker.review_event.set()

    def stop_and_discard(self):
        self.btn_next_image.hide()
        self.btn_stop_review.hide()
        self.btn_confirm_roi.hide() 
        self.btn_run.show()  
        self.viewer.layers.clear()  
        if hasattr(self, "status_label"):
            self.status_label.setText("Stopping...")
        self.worker.request_abort()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ExQt()
    window.showMaximized()
    sys.exit(app.exec())