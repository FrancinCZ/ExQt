import os
os.environ["LIBGL_ALWAYS_SOFTWARE"] = "1"
import sys
import napari
from PySide6.QtWidgets import (QApplication, QMainWindow, QPushButton, 
                               QVBoxLayout, QHBoxLayout, QWidget, QFormLayout, 
                               QSpinBox, QDoubleSpinBox, QComboBox, QCheckBox, QLineEdit, QFileDialog)

class ExQt(QMainWindow): 
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ExQt: Analysis of nuclear condensates")
        self.resize(1200, 800) 

        
        left_panel_layout = QVBoxLayout() 
        folder_layout = QHBoxLayout()
        self.folder_input = QLineEdit("Choose folder with TIF and H5 files")
        self.folder_input.setReadOnly(True)
        self.btn_browse = QPushButton("Browse")
        
        self.btn_browse.clicked.connect(self.choose_folder)
        
        folder_layout.addWidget(self.folder_input)
        folder_layout.addWidget(self.btn_browse)
        left_panel_layout.addLayout(folder_layout)

        form_layout = QFormLayout()

        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["3d", "2d", "single_slice"])
        form_layout.addRow("Process mode:", self.mode_combo)

        self.exp_factor_spin = QDoubleSpinBox()
        self.exp_factor_spin.setValue(4.0)
        self.exp_factor_spin.setDecimals(1)
        form_layout.addRow("Expansion Factor:", self.exp_factor_spin)

        self.threshold_spin = QDoubleSpinBox()
        self.threshold_spin.setRange(0.0, 1.0)
        self.threshold_spin.setSingleStep(0.05)
        self.threshold_spin.setValue(0.3)
        form_layout.addRow("Ilastik Threshold:", self.threshold_spin)

        self.show_napari_check = QCheckBox("Show Napari preview")
        self.show_napari_check.setChecked(True)
        form_layout.addRow("", self.show_napari_check)

        left_panel_layout.addLayout(form_layout)

        self.btn_run = QPushButton("Spustit Batch Zpracování")
        self.btn_run.setStyleSheet("background-color: #2a82da; color: white; padding: 10px; font-weight: bold;")
        left_panel_layout.addWidget(self.btn_run)

        master_layout = QHBoxLayout()

        viewer = napari.Viewer(show=False) 
        master_layout.addWidget(viewer.window._qt_viewer)
        central_widget = QWidget()
        central_widget.setLayout(master_layout) 
        
        self.setCentralWidget(central_widget)

    def choose_folder(self):
        slozka = QFileDialog.getExistingDirectory(self, "Vyberte složku s daty")
        if slozka:
            self.folder_input.setText(slozka)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ExQt()
    window.show()
    sys.exit(app.exec())