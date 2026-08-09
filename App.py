import os
import json
import h5py
import threading
from datetime import datetime
from Batch import process_condensates_h5
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
from postprocessing import generate_excel_stats, generate_plots

DEFAULT_SETTINGS = {
    "adv_sigma": 1.0,
    "adv_pixel_size": 58.0,
    "adv_z_step": 250.0,
    "adv_prob_ch": 0,
    "adv_signal_ch": 1,
    "adv_dapi_ch": 0,
}

class AnalysisWorker(QThread):
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

    def request_roi_callback(self, img_shape, is_3d):
        self.request_roi_signal.emit({"shape": img_shape, "is_3d": is_3d})
        self.roi_event.wait()
        self.roi_event.clear()
        return self.user_roi_data

    def run(self):
        try:
            folder_path = Path(self.params["input_folder"])
            if not folder_path.exists() or not folder_path.is_dir():
                self.progress.emit("Error: Input folder does not exist or is not valid.")
                self.progress.emit("Done_Error")
                return

            raw_files = [f for f in sorted(folder_path.glob("*.tif")) if "Probabilities" not in f.name and "Final" not in f.name]
            all_dataframes = []
            processed_pairs = 0
            had_error = False

            for raw_tif in raw_files:
                h5_file = raw_tif.with_name(f"{raw_tif.stem}_Probabilities.h5")
                if not h5_file.exists():
                    self.progress.emit(f"Skipping {raw_tif.name}: matching H5 not found ({h5_file.name}).")
                    continue

                processed_pairs += 1
                self.progress.emit(f"Processing: {raw_tif.name}")

                try:
                    df_file = process_condensates_h5(
                        tif_path=raw_tif,
                        h5_path=h5_file,
                        mode=self.params["mode"],
                        expansion_factor=self.params["expansion_factor"],
                        prob_threshold=self.params["prob_threshold"],
                        sigma=self.params["sigma"],
                        min_voxels=self.params.get("min_voxels", 5),
                        auto_roi=self.params["auto_roi"],
                        send_layer_func=self.layer_ready.emit if self.params.get("show_napari", True) else None,
                        request_roi_func=self.request_roi_callback if not self.params.get("auto_roi", False) else None,
                        pixel_size_nm=self.params["pixel_size_nm"],
                        z_step_nm=self.params["z_step_nm"],
                        prob_channel=self.params["prob_channel"],
                        signal_channel=self.params["signal_channel"],
                        dapi_channel=self.params["dapi_channel"]
                    )
                    if df_file is not None and not df_file.empty:
                        all_dataframes.append(df_file)
                    else:
                        self.progress.emit(f"No detectable objects found in {raw_tif.name}.")

                    if self.params.get("review_each_image", False):
                        self.progress.emit(f"Review: {raw_tif.name}")
                        self.request_review_signal.emit()
                        self.review_event.wait()
                        self.review_event.clear()
                except Exception as e:
                    self.progress.emit(f"Error processing {raw_tif.name}: {str(e)}")
                    had_error = True

            if not raw_files:
                self.progress.emit("No TIFF files found in input folder.")
                had_error = True
            elif processed_pairs == 0:
                self.progress.emit("No matching TIFF/H5 pairs were processed.")
                had_error = True

            if all_dataframes:
                final_df = __import__("pandas").concat(all_dataframes, ignore_index=True)
                output_folder = Path(self.params["output_folder"]) if self.params["output_folder"] else folder_path
                output_folder.mkdir(parents=True, exist_ok=True)
                output_csv = output_folder / f"{folder_path.name}_Output_Batch_{self.params['mode']}.csv"
                final_df.to_csv(output_csv, index=False)
                self.progress.emit(f"CSV saved: {output_csv.name}")

                metadata = {
                    "timestamp": datetime.now().isoformat(),
                    "software": "ExQt v1.0",
                    "parameters": {
                        "mode": self.params["mode"],
                        "expansion_factor": self.params["expansion_factor"],
                        "prob_threshold": self.params["prob_threshold"],
                        "min_voxels": self.params.get("min_voxels", 5),
                        "gaussian_sigma": self.params["sigma"],
                        "pixel_size_nm": self.params["pixel_size_nm"],
                        "z_step_nm": self.params["z_step_nm"]
                    },
                    "channels": {
                        "prob_channel": self.params["prob_channel"],
                        "signal_channel": self.params["signal_channel"],
                        "dapi_channel": self.params["dapi_channel"]
                    }
                }
                meta_path = output_folder / f"{folder_path.name}_Output_Batch_{self.params['mode']}_metadata.json"
                try:
                    with open(meta_path, "w", encoding="utf-8") as f:
                        json.dump(metadata, f, indent=4)
                    self.progress.emit("Metadata saved.")
                except Exception as e:
                    self.progress.emit(f"Error saving metadata: {str(e)}")

                if self.params.get("gen_excel"):
                    try:
                        generate_excel_stats(str(output_csv))
                        self.progress.emit("Excel was generated.")
                    except Exception as e:
                        self.progress.emit(f"Error generating Excel: {str(e)}")
                        had_error = True

                if self.params.get("gen_plots"):
                    try:
                        generate_plots(
                            str(output_csv), 
                            min_size=self.params.get("plot_min_size", 0.0001), 
                            max_size=self.params.get("plot_max_size", 2.0)
                        )
                        self.progress.emit("Graphs were generated.")
                    except Exception as e:
                        self.progress.emit(f"Graph error: {str(e)}")
                        print(f"DEBUG PLOT ERROR: {str(e)}")
                        return

                # If everything passed (including plots), report success
                self.progress.emit("Done_Success")
            else:
                if not had_error:
                    self.progress.emit("No data available for processing.")
                self.progress.emit("Done_Error")
        except Exception as e:
            self.progress.emit(f"Error when loading: {str(e)}")
            self.progress.emit("Done_Error")

class AdvancedSettingsDialog(QDialog):
    def __init__(self, parent=None, mode="3d"):
        super().__init__(parent)
        self.setWindowTitle("Advanced Settings")
        self.setMinimumWidth(300)

        self.settings = QSettings("MyLab", "ExQt")

        layout = QVBoxLayout()
        form = QFormLayout()

        self.sigma_spin = QDoubleSpinBox()
        self.sigma_spin.setSingleStep(0.5)
        form.addRow("Gaussian Sigma:", self.sigma_spin)

        self.pixel_size_spin = QDoubleSpinBox()
        self.pixel_size_spin.setRange(10.0, 2000.0)
        self.pixel_size_spin.setSingleStep(0.1)
        self.pixel_size_spin.setDecimals(1)
        form.addRow("Pixel Size XY (nm):", self.pixel_size_spin)

        self.z_step_spin = QDoubleSpinBox()
        self.z_step_spin.setRange(10.0, 5000.0)
        self.z_step_spin.setSingleStep(0.1)
        self.z_step_spin.setDecimals(1)
        form.addRow("Z-step (nm):", self.z_step_spin)
        self.z_step_label = form.labelForField(self.z_step_spin)

        self.prob_spin = QSpinBox()
        form.addRow("Probability Channel:", self.prob_spin)

        self.signal_spin = QSpinBox()
        form.addRow("Signal Channel:", self.signal_spin)

        self.dapi_spin = QSpinBox()
        form.addRow("DAPI Channel:", self.dapi_spin)

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
        tooltip = "" if is_3d else "Z-step is only used in 3D mode and is ignored for the current Process mode."
        self.z_step_spin.setToolTip(tooltip)

    def load_adv_settings(self):
        self.sigma_spin.setValue(float(self.settings.value("adv_sigma", DEFAULT_SETTINGS["adv_sigma"])))
        self.pixel_size_spin.setValue(float(self.settings.value("adv_pixel_size", DEFAULT_SETTINGS["adv_pixel_size"])))
        self.z_step_spin.setValue(float(self.settings.value("adv_z_step", DEFAULT_SETTINGS["adv_z_step"])))
        self.prob_spin.setValue(int(self.settings.value("adv_prob_ch", DEFAULT_SETTINGS["adv_prob_ch"])))
        self.signal_spin.setValue(int(self.settings.value("adv_signal_ch", DEFAULT_SETTINGS["adv_signal_ch"])))
        self.dapi_spin.setValue(int(self.settings.value("adv_dapi_ch", DEFAULT_SETTINGS["adv_dapi_ch"])))

    def accept(self):
        self.settings.setValue("adv_sigma", self.sigma_spin.value())
        self.settings.setValue("adv_pixel_size", self.pixel_size_spin.value())
        self.settings.setValue("adv_z_step", self.z_step_spin.value())
        self.settings.setValue("adv_prob_ch", self.prob_spin.value())
        self.settings.setValue("adv_signal_ch", self.signal_spin.value())
        self.settings.setValue("adv_dapi_ch", self.dapi_spin.value())
        super().accept()

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
        self.folder_input.setPlaceholderText("Choose folder with TIF and H5 files...")
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

        self.threshold_spin = QDoubleSpinBox()
        self.threshold_spin.setRange(0.0, 1.0)
        self.threshold_spin.setSingleStep(0.05)
        self.threshold_spin.setValue(0.3)
        form_layout.addRow("Probability Threshold:", self.threshold_spin)

        self.min_voxel_spinbox = QSpinBox()
        self.min_voxel_spinbox.setRange(1, 10000)
        self.min_voxel_spinbox.setValue(5)
        self.min_voxel_spinbox.setToolTip("Ignore objects smaller than this number of voxels/pixels (noise filtering).")
        form_layout.addRow("Min. size (voxels):", self.min_voxel_spinbox)
        self.min_voxel_label = form_layout.labelForField(self.min_voxel_spinbox)

        self.show_napari_check = QCheckBox("Show Napari preview")
        self.show_napari_check.setChecked(True)
        form_layout.addRow("", self.show_napari_check)

        self.auto_roi_check = QCheckBox("Auto-ROI")
        self.auto_roi_check.setToolTip("Auto-ROI treats entire image as a single cell (Per-FOV metrics).")
        form_layout.addRow("", self.auto_roi_check)

        self.review_check = QCheckBox("Pause and review segmentations")
        form_layout.addRow("", self.review_check)

        self.generate_excel_check = QCheckBox("Generate Excel stats")
        form_layout.addRow("", self.generate_excel_check)

        self.generate_plots_check = QCheckBox("Generate plots")
        form_layout.addRow("", self.generate_plots_check)

        self.plot_min_spin = QDoubleSpinBox()
        self.plot_min_spin.setRange(0.0, 10000.0)
        self.plot_min_spin.setDecimals(4)
        self.plot_min_spin.setValue(0.0001)
        
        self.plot_max_spin = QDoubleSpinBox()
        self.plot_max_spin.setRange(0.0, 100000.0)
        self.plot_max_spin.setDecimals(4)
        self.plot_max_spin.setValue(2.0)

        plot_range_layout = QHBoxLayout()
        plot_range_layout.setContentsMargins(0, 0, 0, 0)
        plot_range_layout.addWidget(QLabel("Min:"))
        plot_range_layout.addWidget(self.plot_min_spin)
        plot_range_layout.addWidget(QLabel("Max:"))
        plot_range_layout.addWidget(self.plot_max_spin)

        self.plot_range_widget = QWidget()
        self.plot_range_widget.setLayout(plot_range_layout)
        self.plot_range_widget.hide()

        form_layout.addRow("Plot limits (μm):", self.plot_range_widget)
        self.plot_range_label = form_layout.labelForField(self.plot_range_widget)
        self.plot_range_label.hide()
        self.generate_plots_check.toggled.connect(self.toggle_plot_range)
  
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
        self.load_settings()

        self.btn_run.clicked.connect(self.start_analysis)
        self.mode_combo.currentTextChanged.connect(self.on_mode_changed)
        self.on_mode_changed(self.mode_combo.currentText())  #

    def on_mode_changed(self, mode):
        is_3d = mode == "3d"
        unit_label = "voxels" if is_3d else "pixels"
        if self.min_voxel_label is not None:
            self.min_voxel_label.setText(f"Min. size ({unit_label}):")

    def choose_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Choose folder with data")
        if folder:
            self.folder_input.setText(folder)
    def create_menu(self):
        menu_bar = self.menuBar()

        settings_menu = menu_bar.addMenu("Settings")
        self.dark_mode_action = QAction("Dark Mode", self)
        self.dark_mode_action.setCheckable(True)
        settings_menu.addAction(self.dark_mode_action)
        self.dark_mode_action.toggled.connect(self.switch_theme)

        settings_menu.addSeparator() 

        self.advanced = QAction("Advanced...", self)
        settings_menu.addAction(self.advanced)
        self.advanced.triggered.connect(self.open_advanced_settings)

    def switch_theme(self, active):
        app = QApplication.instance() 
        if active:
            app.setStyleSheet(qdarktheme.load_stylesheet("dark"))
        else:
            app.setStyleSheet(qdarktheme.load_stylesheet("light"))

    def open_advanced_settings(self):
        dialog = AdvancedSettingsDialog(self, mode=self.mode_combo.currentText())
        dialog.exec()
        

    def choose_output_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Choose output folder")
        if folder:
            self.output_path_edit.setText(folder)

    def switch_same_folder(self, checked):
        self.output_path_edit.setDisabled(checked)
        self.btn_output_browse.setDisabled(checked)
        
        if checked:
            self.output_path_edit.setText(self.folder_input.text())
        else:
            self.output_path_edit.clear()

    def load_settings(self):
        self.folder_input.setText(self.settings.value("input_folder", ""))
        self.output_path_edit.setText(self.settings.value("output_folder", ""))
        
        same_folder_saved = self.settings.value("same_folder", "false")
        if str(same_folder_saved).lower() == "true":
            self.same_folder_checkbox.setChecked(True)

    def closeEvent(self, event):
        self.settings.setValue("input_folder", self.folder_input.text())
        self.settings.setValue("output_folder", self.output_path_edit.text())
        self.settings.setValue("same_folder", self.same_folder_checkbox.isChecked())
        event.accept()

    def start_analysis(self):
        params = {
            "input_folder": self.folder_input.text(),
            "output_folder": self.output_path_edit.text(),
            "mode": self.mode_combo.currentText(),
            "expansion_factor": self.exp_factor_spin.value(),
            "prob_threshold": self.threshold_spin.value(),
            "min_voxels": self.min_voxel_spinbox.value(),
            "sigma": float(self.settings.value("adv_sigma", DEFAULT_SETTINGS["adv_sigma"])),
            "auto_roi": self.auto_roi_check.isChecked(),
            "review_each_image": self.review_check.isChecked(),
            "show_napari": self.show_napari_check.isChecked(),
            "gen_excel": self.generate_excel_check.isChecked(),
            "gen_plots": self.generate_plots_check.isChecked(),
            "pixel_size_nm": float(self.settings.value("adv_pixel_size", DEFAULT_SETTINGS["adv_pixel_size"])),
            "z_step_nm": float(self.settings.value("adv_z_step", DEFAULT_SETTINGS["adv_z_step"])),
            "prob_channel": int(self.settings.value("adv_prob_ch", DEFAULT_SETTINGS["adv_prob_ch"])),
            "signal_channel": int(self.settings.value("adv_signal_ch", DEFAULT_SETTINGS["adv_signal_ch"])),
            "dapi_channel": int(self.settings.value("adv_dapi_ch", DEFAULT_SETTINGS["adv_dapi_ch"])),
            "plot_min_size": self.plot_min_spin.value(),
            "plot_max_size": self.plot_max_spin.value()
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
        
    def update_button_text(self, text):
        if text == "Done_Success":
            self.status_label.setText("Analysis completed successfully.")
            self.btn_run.setEnabled(True)
            self.btn_run.setText("Start analysis")
            QMessageBox.information(self, "Done", "Analysis was completed successfully\n\nResults are saved in CSV.")
        elif text == "Done_Error":
            self.btn_run.setEnabled(True)
            self.btn_run.setText("Start analysis")
            QMessageBox.warning(
                self,
                "Error / No data",
                "Processing finished, but no results were obtained!\n\nPlease check:\n1. The image is 3D if you selected mode '3d'.\n2. The TIF and H5 files have exactly matching names."
            )
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
        shape = info["shape"]
        empty_mask = np.zeros(shape, dtype=int)
        self.viewer.add_labels(empty_mask, name="Paint ROI", opacity=0.5)
        self.viewer.layers["Paint ROI"].mode = 'paint'
        self.btn_run.hide()
        self.btn_confirm_roi.show()

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
        self.btn_run.show()
        self.worker.roi_event.set()

    def prepare_review(self):
        self.btn_run.hide()
        self.btn_confirm_roi.hide()
        self.btn_next_image.show()
        if hasattr(self, "status_label"):
            self.status_label.setText("Waiting for user review...")

    def toggle_plot_range(self, checked):
        self.plot_range_widget.setVisible(checked)
        if self.plot_range_label:
            self.plot_range_label.setVisible(checked)

    def next_image_confirmed(self):
        self.btn_next_image.hide()
        self.btn_run.show()
        if hasattr(self, "status_label"):
            self.status_label.setText("Processing next...")
        self.worker.review_event.set()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ExQt()
    window.showMaximized()
    sys.exit(app.exec())
