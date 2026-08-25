import os
import json
import threading
from datetime import datetime
from Batch import (
    MODE_A_Z_SPLIT_MIN_COMPONENT_FRACTION,
    MODE_A_Z_SPLIT_MIN_COMPONENT_VOXELS,
    MODE_A_Z_SPLIT_PASS_FRACTION,
    MODE_A_Z_SPLIT_REVIEW_FRACTION,
    get_metadata_from_tif,
    process_condensates,
)
from rezim_a_metrics import MODE_A_LAYER_SCHEME
import sys
import napari
from PySide6.QtWidgets import (QApplication, QLabel, QMainWindow, QPushButton, 
                                QVBoxLayout, QHBoxLayout, QWidget, QFormLayout, 
                                QSpinBox, QDoubleSpinBox, QComboBox, QCheckBox, QLineEdit, QFileDialog, QDialog, QDialogButtonBox,
                                QMessageBox)
from PySide6.QtGui import QAction
import qdarktheme
from PySide6.QtCore import QSettings, QThread, Signal
import numpy as np
import tifffile
import pandas as pd
from pathlib import Path
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from postprocessing import (
    QCPolicyMismatchError,
    generate_excel_stats,
    generate_plots,
    generate_rezim_a_plots,
    merge_statistics_folder,
)
from calibration_policy import summarize_calibrations
from size_preview import collect_size_preview

DEFAULT_SETTINGS = {
    "adv_pixel_size": 58.0,
    "adv_z_step": 250.0,
    "adv_signal_ch": 1,
    "adv_dapi_ch": 0,
    "raw_min_voxels": 5,
    "mode_a_enabled": False,
    "mode_a_min_core_voxels": 20,
    "mode_a_exclude_split_slices": True,
    "mode_a_z_split_min_component_voxels": MODE_A_Z_SPLIT_MIN_COMPONENT_VOXELS,
    "mode_a_z_split_min_component_fraction": MODE_A_Z_SPLIT_MIN_COMPONENT_FRACTION,
    "mode_a_z_split_pass_fraction": MODE_A_Z_SPLIT_PASS_FRACTION,
    "mode_a_z_split_review_fraction": MODE_A_Z_SPLIT_REVIEW_FRACTION,
}

REPORT_DEFAULTS = {
    "report_excel": True,
    "report_primary_csv": True,
    "report_excluded_csv": True,
    "report_raw_audit_csv": False,
    "report_standard_plots": True,
    "report_mode_a_plots": True,
}

def _safe_float(value, default):
    if value is None:
        return default
    #Parse  numeric settings without breaking GUI startup.
    try:
        if isinstance(value, str):
            value = value.replace(',', '.')
        return float(value)
    except (ValueError, TypeError):
        return default


#Parse QSettings boolean values from either native or text forms.
def _safe_bool(value, default=False):
    if value is None:
        return default
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)

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
        self.abort_requested = False

    # Pause the worker until the GUI returns a painted ROI mask.
    def request_roi_callback(self, img_shape, is_3d):
        self.request_roi_signal.emit({"shape": img_shape, "is_3d": is_3d})
        self.roi_event.wait()
        self.roi_event.clear()
        return self.user_roi_data

    # Wake any pending ROI/review wait and mark the batch for shutdown.
    def request_abort(self):
        self.abort_requested = True
        self.review_event.set()
        self.roi_event.set() 

    # Validate the folder, delegate measurements to Batch, and write outputs.
    def run(self):
        try:
            folder_path = Path(self.params["input_folder"])
            if not folder_path.exists() or not folder_path.is_dir():
                self.progress.emit("Error: Input folder does not exist or is not valid.")
                self.progress.emit("Done")
                return

            #Pair each source TIFF with its convention-based *_Mask.tif file;
            raw_files = [f for f in sorted(folder_path.glob("*.tif")) if "Mask" not in f.name and "Final" not in f.name]
            all_dataframes = []

            for raw_tif in raw_files:
                if self.abort_requested:
                    self.progress.emit("Analysis stopped by user.")
                    break

                mask_file = raw_tif.with_name(f"{raw_tif.stem}_Mask.tif")
                if not mask_file.exists():
                    self.progress.emit(f"Skipping {raw_tif.name}: Mask not found.")
                    continue

                self.progress.emit(f"Processing: {raw_tif.name}")

                try:
                    #Batch owns image measurements; the worker only supplies
                    df_file = process_condensates(
                        tif_path=raw_tif,
                        mask_path=mask_file,
                        mode=self.params["mode"],
                        expansion_factor=self.params["expansion_factor"],
                        min_voxels=self.params.get("min_voxels", 5),
                        auto_roi=self.params["auto_roi"],
                        send_layer_func=self.layer_ready.emit if self.params.get("show_napari", True) else None,
                        request_roi_func=self.request_roi_callback if not self.params.get("auto_roi", False) else None,
                        pixel_size_nm=self.params["pixel_size_nm"],
                        z_step_nm=self.params["z_step_nm"],
                        signal_channel=self.params["signal_channel"],
                        dapi_channel=self.params["dapi_channel"],
                        mode_a_enabled=self.params.get("mode_a_enabled", False),
                        mode_a_min_core_voxels=self.params.get("mode_a_min_core_voxels", 20),
                        mode_a_exclude_split_slices=self.params.get("mode_a_exclude_split_slices", True),
                        mode_a_z_split_min_component_voxels=self.params.get("mode_a_z_split_min_component_voxels", MODE_A_Z_SPLIT_MIN_COMPONENT_VOXELS),
                        mode_a_z_split_min_component_fraction=self.params.get("mode_a_z_split_min_component_fraction", MODE_A_Z_SPLIT_MIN_COMPONENT_FRACTION),
                        mode_a_z_split_pass_fraction=self.params.get("mode_a_z_split_pass_fraction", MODE_A_Z_SPLIT_PASS_FRACTION),
                        mode_a_z_split_review_fraction=self.params.get("mode_a_z_split_review_fraction", MODE_A_Z_SPLIT_REVIEW_FRACTION),
                    )
                    appended_this_round = False
                    if df_file is not None and not df_file.empty:
                        all_dataframes.append(df_file)
                        appended_this_round = True

                    #Review happens after a file is computed
                    if self.params.get("review_each_image", False):
                        self.progress.emit(f"Review: {raw_tif.name}")
                        self.request_review_signal.emit()
                        self.review_event.wait()
                        self.review_event.clear()

                        if self.abort_requested:
                            if appended_this_round:
                                all_dataframes.pop() 
                            self.progress.emit(f"Discarded and stopped at: {raw_tif.name}")
                            break
                except Exception as e:
                    self.progress.emit(f"Error: {str(e)}")

            #Concatenate approved files once, then create the CSV 
            if all_dataframes:
                final_df = __import__("pandas").concat(all_dataframes, ignore_index=True)
                output_folder = Path(self.params["output_folder"]) if self.params["output_folder"] else folder_path
                output_folder.mkdir(parents=True, exist_ok=True)
                output_csv = output_folder / f"{folder_path.name}_Output_Batch_{self.params['mode']}.csv"
                final_df.to_csv(output_csv, index=False)
                self.progress.emit(f"CSV saved: {output_csv.name}")

                #Store calibration, mode, and interpretation limits in a JSON alongside the CSV for reproducibility and future reference.
                metadata = {
                    "timestamp": datetime.now().isoformat(),
                    "software": "ExQt",
                    "parameters": {
                        "mode": self.params["mode"],
                        "expansion_factor": self.params["expansion_factor"],
                        "min_voxels": self.params.get("min_voxels", 5),
                        "pixel_size_nm": self.params["pixel_size_nm"],
                        "z_step_nm": self.params["z_step_nm"],
                        # Keep the applied and detected calibration evidence together in metadata JSON.
                        "calibration_policy": self.params.get("calibration_policy", "one_explicit_calibration_per_batch"),
                        "calibration_confirmation": self.params.get("calibration_confirmation", "unknown"),
                        "detected_metadata_by_file": self.params.get("detected_metadata_by_file", {}),
                        "plot_min_size": self.params.get("plot_min_size", 0.0001),
                        "plot_max_size": self.params.get("plot_max_size", 2.0)
                    },
                    "channels": {
                        "signal_channel": self.params["signal_channel"],
                        "dapi_channel": self.params["dapi_channel"]
                    },
                    "mode_a": {
                        "enabled": self.params.get("mode_a_enabled", False),
                        "available_in_mode": "3d",
                        "min_core_voxels": self.params.get("mode_a_min_core_voxels", 20),
                        "minimum_voxel_policy": "same minimum applied to object, shell, middle, and core FA validity; also used as minimum core size",
                        "exclude_split_slices_from_primary": self.params.get("mode_a_exclude_split_slices", True),
                        "require_z_topology_pass_for_primary": self.params.get("mode_a_exclude_split_slices", True),
                        "z_split_policy": "substantial_component_fraction_v1",
                        "z_split_min_component_voxels": self.params.get("mode_a_z_split_min_component_voxels", MODE_A_Z_SPLIT_MIN_COMPONENT_VOXELS),
                        "z_split_min_component_fraction": self.params.get("mode_a_z_split_min_component_fraction", MODE_A_Z_SPLIT_MIN_COMPONENT_FRACTION),
                        "z_split_pass_fraction": self.params.get("mode_a_z_split_pass_fraction", MODE_A_Z_SPLIT_PASS_FRACTION),
                        "z_split_review_fraction": self.params.get("mode_a_z_split_review_fraction", MODE_A_Z_SPLIT_REVIEW_FRACTION),
                        "layer_scheme": MODE_A_LAYER_SCHEME,
                        "sampling_order": "Z,Y,X",
                        "metrics": [
                            "A_object", "A_shell", "A_middle", "A_core",
                            "Delta_A_middle_shell", "Delta_A_core_middle", "Delta_A_core_shell",
                        ],
                        "interpretation": "Geometric structural-response metrics; not direct stiffness, liquidity, or viscosity measurements."
                    },
                    "reporting": {
                        "enabled": self.params.get("generate_reports", False),
                        "excel": self.params.get("report_excel", True),
                        "primary_csv": self.params.get("report_primary_csv", True),
                        "excluded_csv": self.params.get("report_excluded_csv", True),
                        "raw_audit_csv": self.params.get("report_raw_audit_csv", False),
                        "standard_plots": self.params.get("report_standard_plots", True),
                        "mode_a_plots": self.params.get("report_mode_a_plots", True),
                    },
                }
                meta_path = output_folder / f"{folder_path.name}_Output_Batch_{self.params['mode']}_metadata.json"
                try:
                    with open(meta_path, "w", encoding="utf-8") as f:
                        json.dump(metadata, f, indent=4)
                    self.progress.emit("Metadata saved.")
                except Exception as e:
                    self.progress.emit(f"Error saving metadata: {str(e)}")

                report_tables_requested = self.params.get("generate_reports", False) and any((
                    self.params.get("report_excel", True),
                    self.params.get("report_primary_csv", True),
                    self.params.get("report_excluded_csv", True),
                    self.params.get("report_raw_audit_csv", False),
                ))
                if report_tables_requested:
                    try:
                        report_result = generate_excel_stats(
                            str(output_csv),
                            min_size=self.params.get("plot_min_size", 0.0001),
                            max_size=self.params.get("plot_max_size", 2.0),
                            generate_excel=self.params.get("report_excel", True),
                            generate_primary_csv=self.params.get("report_primary_csv", True),
                            generate_excluded_csv=self.params.get("report_excluded_csv", True),
                            generate_raw_audit_csv=self.params.get("report_raw_audit_csv", False),
                        )
                        generated = [name for name in ("excel", "primary_csv", "excluded_csv", "all_objects_csv") if report_result.get(name)]
                        self.progress.emit(f"Selected tables generated: {', '.join(generated)}")
                    except Exception as e:
                        self.progress.emit(f"Error generating report tables: {str(e)}")

                if self.params.get("generate_reports", False) and self.params.get("report_standard_plots", True):
                    try:
                        generate_plots(
                            str(output_csv),
                            min_size=self.params.get("plot_min_size", 0.0001),
                            max_size=self.params.get("plot_max_size", 2.0),
                        )
                        self.progress.emit("Graphs were generated.")
                    except Exception as e:
                        self.progress.emit(f"Error generating graphs: {str(e)}")

                # Radial FA Profiling creates an additional audit report. Standard ExQt plots above remain unchanged and do not contain radial FA values.
                if (
                    self.params.get("generate_reports", False)
                    and self.params.get("report_mode_a_plots", True)
                    and self.params.get("mode_a_enabled")
                    and self.params["mode"] == "3d"
                ):
                    try:
                        generate_rezim_a_plots(
                            str(output_csv),
                            min_size=self.params.get("plot_min_size"),
                            max_size=self.params.get("plot_max_size"),
                        )
                        self.progress.emit("Radial FA Profiling plots were generated.")
                    except Exception as e:
                        self.progress.emit(f"Error generating Radial FA Profiling plots: {str(e)}")
            else:
                self.progress.emit("No data available for processing.")

            self.progress.emit("Done")
        except Exception as e:
            self.progress.emit(f"Error when loading: {str(e)}")
            self.progress.emit("Done")

class ReportOptionsDialog(QDialog):
    #Choose derived outputs without changing the source audit CSV

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
            # Preserve the original quick-click behaviour away from a handle.
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
    #Edit and persist calibration, channel, and Radial FA Profiling settings.
    def __init__(self, parent=None, mode="3d"):
        super().__init__(parent)
        self.setWindowTitle("Advanced Settings")
        self.setMinimumWidth(300)

        self.settings = QSettings("MyLab", "ExQt")

        layout = QVBoxLayout()
        form = QFormLayout()

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

        self.raw_min_voxels_spin = QSpinBox()
        self.raw_min_voxels_spin.setRange(1, 10000)
        self.raw_min_voxels_spin.setToolTip(
            "Early connected-component noise floor before calibrated biological-size filtering. "
            "This is separate from the analyzed µm²/µm³ range and from Radial FA Profiling layer validity."
        )
        form.addRow("Raw noise filter (pixels/voxels):", self.raw_min_voxels_spin)

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
        form.addRow("Radial FA Profiling min. voxels per layer:", self.mode_a_min_core_spin)
        self.mode_a_label = form.labelForField(self.mode_a_min_core_spin)

        self.mode_a_exclude_split_check = QCheckBox(
            "Require Z-topology PASS for primary comparison"
        )
        self.mode_a_exclude_split_check.setToolTip(
            "Ignores tiny detached components, grades persistent substantial splitting as PASS/REVIEW/FAIL, "
            "and keeps REVIEW/FAIL objects in the CSV outside the primary comparison."
        )
        form.addRow("", self.mode_a_exclude_split_check)

        layout.addLayout(form)

        self.load_adv_settings()
        self._apply_mode_state(mode)

        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)

        self.setLayout(layout)

    #Enable controls only when the selected processing mode supports them.
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

    #Populate controls from QSettings, applying safe defaults.
    def load_adv_settings(self):
        self.pixel_size_spin.setValue(_safe_float(self.settings.value("adv_pixel_size", DEFAULT_SETTINGS["adv_pixel_size"]), DEFAULT_SETTINGS["adv_pixel_size"]))
        self.z_step_spin.setValue(_safe_float(self.settings.value("adv_z_step", DEFAULT_SETTINGS["adv_z_step"]), DEFAULT_SETTINGS["adv_z_step"]))
        self.signal_spin.setValue(int(self.settings.value("adv_signal_ch", DEFAULT_SETTINGS["adv_signal_ch"])))
        self.dapi_spin.setValue(int(self.settings.value("adv_dapi_ch", DEFAULT_SETTINGS["adv_dapi_ch"])))
        self.raw_min_voxels_spin.setValue(int(self.settings.value("raw_min_voxels", DEFAULT_SETTINGS["raw_min_voxels"])))
        self.mode_a_enabled_check.setChecked(_safe_bool(self.settings.value("mode_a_enabled", DEFAULT_SETTINGS["mode_a_enabled"])))
        self.mode_a_min_core_spin.setValue(int(self.settings.value("mode_a_min_core_voxels", DEFAULT_SETTINGS["mode_a_min_core_voxels"])))
        self.mode_a_exclude_split_check.setChecked(_safe_bool(self.settings.value("mode_a_exclude_split_slices", DEFAULT_SETTINGS["mode_a_exclude_split_slices"])))

    #Persist dialog values before closing successfully.
    def accept(self):
        self.settings.setValue("adv_pixel_size", self.pixel_size_spin.value())
        self.settings.setValue("adv_z_step", self.z_step_spin.value())
        self.settings.setValue("adv_signal_ch", self.signal_spin.value())
        self.settings.setValue("adv_dapi_ch", self.dapi_spin.value())
        self.settings.setValue("raw_min_voxels", self.raw_min_voxels_spin.value())
        self.settings.setValue("mode_a_enabled", self.mode_a_enabled_check.isChecked())
        self.settings.setValue("mode_a_min_core_voxels", self.mode_a_min_core_spin.value())
        self.settings.setValue("mode_a_exclude_split_slices", self.mode_a_exclude_split_check.isChecked())
        super().accept()

class ExQt(QMainWindow): 
    #Main window connecting user controls, Napari, and AnalysisWorker.
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

        self.same_folder_checkbox = QCheckBox("Save to the same folder as input")
        self.same_folder_checkbox.toggled.connect(self.switch_same_folder)
        self.left_panel_layout.addWidget(self.same_folder_checkbox)
        
        form_layout = QFormLayout()

        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["3d", "2d", "single_slice"])
        form_layout.addRow("Process mode:", self.mode_combo)

        self.exp_factor_spin = QDoubleSpinBox()
        self.exp_factor_spin.setMinimum(0.1)
        self.exp_factor_spin.setValue(4.0)
        self.exp_factor_spin.setDecimals(1)
        form_layout.addRow("Expansion Factor:", self.exp_factor_spin)

        self.show_napari_check = QCheckBox("Show Napari preview")
        self.show_napari_check.setChecked(True)
        form_layout.addRow("", self.show_napari_check)

        self.auto_roi_check = QCheckBox("Auto-ROI")
        self.auto_roi_check.setToolTip("Auto-ROI treats entire image as a single cell (Per-FOV metrics).")
        form_layout.addRow("", self.auto_roi_check)

        self.review_check = QCheckBox("Pause and review segmentations")
        form_layout.addRow("", self.review_check)

        report_row = QHBoxLayout()
        self.generate_reports_check = QCheckBox("Generate selected reports")
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
        form_layout.addRow("Analyzed Size Min (µm²/µm³):", self.plot_min_size_spin)

        self.plot_max_size_spin = QDoubleSpinBox()
        self.plot_max_size_spin.setRange(0.0001, 100000.0)
        self.plot_max_size_spin.setDecimals(5)
        self.plot_max_size_spin.setToolTip(
            "Upper calibrated biological-size bound used by primary statistics, clean CSVs and plots. "
            "The original raw audit CSV remains unfiltered."
        )
        form_layout.addRow("Analyzed Size Max (µm²/µm³):", self.plot_max_size_spin)

        self.size_preview_button = QPushButton("Preview size distribution...")
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

        self.status_label = QLabel("Prepared")
        self.status_label.setStyleSheet("color: gray; font-weight: bold; margin-top: 10px;")
        self.left_panel_layout.addWidget(self.status_label)

        master_layout = QHBoxLayout()

        self.viewer = napari.Viewer(show=False) 
        master_layout.addLayout(self.left_panel_layout, stretch=1)
        master_layout.addWidget(self.viewer.window._qt_window, stretch=3)

        central_widget = QWidget()
        central_widget.setLayout(master_layout) 
        self.setCentralWidget(central_widget)

        self.settings = QSettings("MyLab", "ExQt")
        # TIFF metadata is advisory only; the user controls the calibration in Advanced Settings.
        self.detected_metadata_by_file = {}
        self.calibration_warning = ""
        self.metadata_scanned_folder = None
        self.load_settings()

        self.btn_run.clicked.connect(self.start_analysis)
        self.mode_combo.currentTextChanged.connect(self.on_mode_changed)
        self.on_mode_changed(self.mode_combo.currentText()) 

    #Update controls when the processing mode changes; Radial FA Profiling is 3D-only.
    def on_mode_changed(self, mode):
        pass

    #Select the input folder and opportunistically load TIFF calibration.
    def choose_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Choose folder with data")
        if folder:
            self.folder_input.setText(folder)
            self._try_auto_fill_metadata(folder)

    #Inspect every source TIFF, but never silently overwrite the user's calibration.
    def _try_auto_fill_metadata(self, folder_path):
        folder = Path(folder_path)
        self.metadata_scanned_folder = str(folder.resolve())
        files = sorted(
            f for f in folder.glob("*.tif")
            if "Mask" not in f.name and "Final" not in f.name
        )
        self.detected_metadata_by_file = {}
        self.calibration_warning = ""

        for tif_path in files:
            meta = get_metadata_from_tif(tif_path)
            if meta and "pixel_size" in meta and "z_step" in meta:
                self.detected_metadata_by_file[tif_path.name] = {
                    "pixel_size_nm": float(meta["pixel_size"]),
                    "z_step_nm": float(meta["z_step"]),
                    "sources": meta.get("sources", {}),
                }

        summary = summarize_calibrations(self.detected_metadata_by_file)
        if summary["has_mismatch"]:
            self.calibration_warning = (
                "Input files contain different physical calibrations. "
                "Set one explicit calibration for this batch or split the batch."
            )
            self.status_label.setText("Calibration mismatch detected — review Advanced Settings.")
            QMessageBox.warning(self, "Calibration mismatch", self.calibration_warning)
        elif summary["calibration_count"] == 1:
            pixel_size_nm, z_step_nm = summary["calibrations"][0]
            reply = QMessageBox.question(
                self,
                "Use detected TIFF calibration?",
                (
                    "All readable input files report the same calibration:\n\n"
                    f"XY pixel size: {pixel_size_nm:g} nm\n"
                    f"Z-step: {z_step_nm:g} nm\n\n"
                    "Use these values in ExQt Advanced Settings? They will still be shown "
                    "for confirmation before analysis starts."
                ),
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.Yes,
            )
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
                "No usable physical metadata detected; enter calibration manually in Advanced Settings."
            )
            
    def create_menu(self):
        menu_bar = self.menuBar()

        tools_menu = menu_bar.addMenu("Tools")
        self.merge_runs_action = QAction("Merge existing runs...", self)
        self.merge_runs_action.setToolTip(
            "Recursively merge compatible ExQt batch CSVs from one folder tree."
        )
        tools_menu.addAction(self.merge_runs_action)
        self.merge_runs_action.triggered.connect(self.merge_existing_runs)

        settings_menu = menu_bar.addMenu("Settings")
        self.dark_mode_action = QAction("Dark Mode", self)
        self.dark_mode_action.setCheckable(True)
        settings_menu.addAction(self.dark_mode_action)
        self.dark_mode_action.toggled.connect(self.switch_theme)

        settings_menu.addSeparator() 

        self.advanced = QAction("Advanced...", self)
        settings_menu.addAction(self.advanced)
        self.advanced.triggered.connect(self.open_advanced_settings)

    #Apply the selected qdarktheme stylesheet to the Qt application.
    def switch_theme(self, active):
        app = QApplication.instance() 
        if active:
            app.setStyleSheet(qdarktheme.load_stylesheet("dark"))
        else:
            app.setStyleSheet(qdarktheme.load_stylesheet("light"))

    #Open calibration and Radial FA Profiling settings for the current process mode.
    def open_advanced_settings(self):
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

    #Select the destination used for generated CSV and derived outputs.
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
        source_files = sorted(
            path for path in folder.glob("*.tif")
            if "Mask" not in path.name and "Final" not in path.name
        )
        if not source_files:
            QMessageBox.warning(self, "Size preview", "No source TIF files were found in the selected folder.")
            return
        missing_masks = [
            path.name for path in source_files
            if not (folder / f"{path.stem}_Mask.tif").exists()
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

    #Synchronize output selection with the input folder when requested.
    def switch_same_folder(self, checked):
        self.output_path_edit.setDisabled(checked)
        self.btn_output_browse.setDisabled(checked)
        
        if checked:
            self.output_path_edit.setText(self.folder_input.text())
        else:
            self.output_path_edit.clear()

    #Restore folder and plotting preferences saved by the last session.
    def load_settings(self):
        self.folder_input.setText(self.settings.value("input_folder", ""))
        self.output_path_edit.setText(self.settings.value("output_folder", ""))
        
        self.plot_min_size_spin.setValue(_safe_float(self.settings.value("plot_min_size", 0.0001), 0.0001))
        self.plot_max_size_spin.setValue(_safe_float(self.settings.value("plot_max_size", 2.0), 2.0))
        self.generate_reports_check.setChecked(_safe_bool(self.settings.value("generate_reports", False)))
        
        same_folder_saved = self.settings.value("same_folder", "false")
        if str(same_folder_saved).lower() == "true":
            self.same_folder_checkbox.setChecked(True)

    #Persist user-facing paths and plot limits before the window closes.
    def closeEvent(self, event):
        self.settings.setValue("input_folder", self.folder_input.text())
        self.settings.setValue("output_folder", self.output_path_edit.text())
        self.settings.setValue("same_folder", self.same_folder_checkbox.isChecked())
        self.settings.setValue("plot_min_size", self.plot_min_size_spin.value())
        self.settings.setValue("plot_max_size", self.plot_max_size_spin.value())
        self.settings.setValue("generate_reports", self.generate_reports_check.isChecked())
        event.accept()

    #Validate selections, collect settings, and start the worker thread.
    def start_analysis(self):
        folder_path = self.folder_input.text()
        if not folder_path or not os.path.isdir(folder_path):
            QMessageBox.warning(self, "Error", "Please select a valid input folder first.")
            return

        folder = Path(folder_path)
        files = [f for f in folder.glob("*.tif") if "Mask" not in f.name and "Final" not in f.name]
        if not files:
            QMessageBox.warning(self, "Error", "In the selected folder, no source TIF files were found.")
            return
        missing_masks = [f.name for f in files if not (folder / f"{f.stem}_Mask.tif").exists()]

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
        if self.mode_combo.currentText() == "3d":
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
            QMessageBox.No,
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
            "report_excel": _safe_bool(self.settings.value("report_excel", REPORT_DEFAULTS["report_excel"]), REPORT_DEFAULTS["report_excel"]),
            "report_primary_csv": _safe_bool(self.settings.value("report_primary_csv", REPORT_DEFAULTS["report_primary_csv"]), REPORT_DEFAULTS["report_primary_csv"]),
            "report_excluded_csv": _safe_bool(self.settings.value("report_excluded_csv", REPORT_DEFAULTS["report_excluded_csv"]), REPORT_DEFAULTS["report_excluded_csv"]),
            "report_raw_audit_csv": _safe_bool(self.settings.value("report_raw_audit_csv", REPORT_DEFAULTS["report_raw_audit_csv"]), REPORT_DEFAULTS["report_raw_audit_csv"]),
            "report_standard_plots": _safe_bool(self.settings.value("report_standard_plots", REPORT_DEFAULTS["report_standard_plots"]), REPORT_DEFAULTS["report_standard_plots"]),
            "report_mode_a_plots": _safe_bool(self.settings.value("report_mode_a_plots", REPORT_DEFAULTS["report_mode_a_plots"]), REPORT_DEFAULTS["report_mode_a_plots"]),
            "plot_min_size": self.plot_min_size_spin.value(),
            "plot_max_size": self.plot_max_size_spin.value(),
            "pixel_size_nm": pixel_size_nm,
            "z_step_nm": z_step_nm,
            #The normal GUI workflow uses one deliberate calibration per batch.
            "calibration_policy": "one_explicit_calibration_per_batch",
            "calibration_confirmation": "confirmed_by_user_at_run_start",
            "detected_metadata_by_file": self.detected_metadata_by_file,
            "signal_channel": int(self.settings.value("adv_signal_ch", DEFAULT_SETTINGS["adv_signal_ch"])),
            "dapi_channel": int(self.settings.value("adv_dapi_ch", DEFAULT_SETTINGS["adv_dapi_ch"])),
            "mode_a_enabled": _safe_bool(self.settings.value("mode_a_enabled", DEFAULT_SETTINGS["mode_a_enabled"])),
            "mode_a_min_core_voxels": int(self.settings.value("mode_a_min_core_voxels", DEFAULT_SETTINGS["mode_a_min_core_voxels"])),
            "mode_a_exclude_split_slices": _safe_bool(self.settings.value("mode_a_exclude_split_slices", DEFAULT_SETTINGS["mode_a_exclude_split_slices"])),
            "mode_a_z_split_min_component_voxels": DEFAULT_SETTINGS["mode_a_z_split_min_component_voxels"],
            "mode_a_z_split_min_component_fraction": DEFAULT_SETTINGS["mode_a_z_split_min_component_fraction"],
            "mode_a_z_split_pass_fraction": DEFAULT_SETTINGS["mode_a_z_split_pass_fraction"],
            "mode_a_z_split_review_fraction": DEFAULT_SETTINGS["mode_a_z_split_review_fraction"],
        }

        self.btn_run.setEnabled(False)
        self.btn_run.setText("Processing...")
        self.viewer.layers.clear()

        self.worker = AnalysisWorker(params)
        self.worker.layer_ready.connect(self.receive_layer)
        self.worker.progress.connect(self.update_button_text)
        self.worker.request_roi_signal.connect(self.prepare_manual_roi)
        self.worker.request_review_signal.connect(self.prepare_review)
        self.worker.start()
        
    #Map worker progress signals onto status text and button state.
    def update_button_text(self, text):
        print(f"GUI LOG: {text}") 
        
        if text == "Done":
            self.btn_stop_review.hide()
            self.btn_run.setEnabled(True)
            self.btn_run.setText("Start analysis")
            
            if not self.status_label.text().startswith("Error") and not self.status_label.text().startswith("No data"):
                self.status_label.setText("Analysis completed successfully.")
                QMessageBox.information(self, "Done", "Analysis was completed successfully\n\nResults are saved in CSV.")
        elif text.startswith("Error") or text.startswith("No data"):
            self.status_label.setText(text)
            QMessageBox.warning(self, "Analysis Failed", text) 
        else:
            self.status_label.setText(text)
            self.btn_run.setText("Processing...")

    #Translate Batch preview messages into Napari layers.
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

    #Present a paintable ROI layer while the worker waits for input.
    def prepare_manual_roi(self, info):
        shape = info["shape"]
        empty_mask = np.zeros(shape, dtype=int)
        self.viewer.add_labels(empty_mask, name="Paint ROI", opacity=0.5)
        self.viewer.layers["Paint ROI"].mode = 'paint'
        self.btn_run.hide()
        self.btn_confirm_roi.show()
        self.btn_stop_review.show() 

    #Return the painted ROI to the worker and resume processing.
    def confirm_roi(self):
        if "Paint ROI" in self.viewer.layers:
            mask_data = self.viewer.layers["Paint ROI"].data
            self.worker.user_roi_data = mask_data
            self.viewer.layers.remove("Paint ROI")
        else:
            first_layer = next(iter(self.viewer.layers), None)
            if first_layer is not None:
                self.worker.user_roi_data = np.ones(first_layer.data.shape, dtype=int)
            else:
                self.worker.user_roi_data = np.ones((1, 1), dtype=int)

        self.btn_confirm_roi.hide()
        self.btn_stop_review.hide()
        self.btn_run.show()
        self.worker.roi_event.set()

    #Pause between images so the user can approve the current preview.
    def prepare_review(self):
        self.btn_run.hide()
        self.btn_confirm_roi.hide()
        self.btn_next_image.show()
        self.btn_stop_review.show()
        if hasattr(self, "status_label"):
            self.status_label.setText("Waiting for user review...")

    #Approve the current image and release the worker for the next one.
    def next_image_confirmed(self):
        self.btn_next_image.hide()
        self.btn_stop_review.hide()
        self.btn_run.show()
        if hasattr(self, "status_label"):
            self.status_label.setText("Processing next...")
        self.worker.review_event.set()

    #Stop the batch and discard the image currently under review.
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
    #Keep Qt application startup in the script entry point so imports remain
    app = QApplication(sys.argv)
    window = ExQt()
    window.showMaximized()
    sys.exit(app.exec())
