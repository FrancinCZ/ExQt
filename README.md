# ExQt: Expansion Microscopy Quantification Tool

**ExQt** (*Expansion Microscopy Quantification Tool*) is a desktop application for measuring objects in 3D confocal and Expansion Microscopy (ExM) images of the cell nucleus. It does not segment images itself: it takes existing segmentation masks (e.g. from Labkit or ilastik) and measures the size, enrichment (K_part) and shape (FA) of every object, optionally also for three concentric layers. Each run records its settings, code version and input files.

 **[📖 Read the Documentation](https://francincz.github.io/ExQt/)**

---

## Quick Installation

ExQt needs **Python 3.10 or newer** (tested with Python 3.14) and runs on Windows, Linux and macOS. Open a terminal (Command Prompt on Windows) and run:

```bash
git clone https://github.com/FrancinCZ/ExQt.git
cd ExQt
python -m pip install -r requirements.txt
python App.py
```

Without Git: on GitHub click **Code → Download ZIP**, unzip it, open a terminal in the unzipped folder and run the last two commands.

With Conda instead of pip: `conda env create -f environment.yml`, then `conda activate exqt-env` and `python App.py`.

---

## User Guide: 4 Steps to Run an Analysis

### Step 1: Prepare Your Input Files
Place your raw 3D fluorescence stacks and their corresponding segmentation masks in the same folder. ExQt pairs them automatically by filename:

```text
My_Experiment/
├── cell_01.tif          <-- Raw 3D fluorescence image
└── cell_01_Mask.tif     <-- Binary or labeled mask (from ilastik, Labkit, etc.)
```

### Step 2: Set Calibration
1. Select your **Input folder** and **Output folder**.
2. Click **Settings** in the menu bar and enter the microscope calibration (values before expansion correction):
   - **Pixel Size XY (nm)**, e.g. `65.0`, and **Z-step (nm)**, e.g. `200.0` (ExQt offers the calibration stored in OME/ImageJ TIFFs);
   - **Detector offset (ADU)**: `0` for photon-counting detectors; **From .lif…** reads it from a Leica file.
3. In the main window set the **Expansion factor** (`1.0` for classical confocal) and the **Min size / Max size** range used for the primary set (**Preview sizes…** shows the distribution first).

### Step 3: Define the Nuclear ROI (Napari)
K_part needs the nucleoplasm around the condensates, so the nucleus is outlined by hand:
- When prompted, draw one shape per nucleus in the integrated **Napari** viewer, on one or more Z-slices.
- Shapes of one nucleus on several slices are interpolated between them; a nucleus drawn on a single slice is extended through the whole stack.
- With **Auto-ROI** the whole field of view is used and K_part is not computed.

### Step 4: Run & Inspect Results
Click **Start analysis**. ExQt analyzes all matched pairs one by one and writes the results into the output folder:

- **`<run>_Primary_Condensates.csv`** — Objects that passed all QC criteria (the primary set).
- **`<run>_QC_Excluded.csv`** — Size-eligible objects excluded by QC, with explicit rejection reasons.
- **`<folder>_Output_Batch_3d.csv`** — All measured objects with every parameter (the source table).
- **`<folder>_Output_Batch_3d_metadata.json`** — Settings, code version and input file hashes needed to reproduce the run.
- **`<run>_Detailed_Stats.xlsx`** — Excel workbook with summary, primary and excluded objects.
- **`<run>_partitioning_analysis.png`** and **`<run>_ExQt_Report.pdf`** — Four-panel summary figure and report.
- **`<stem>_ROI.tif`** — The ROI used for each image.

---

## Detailed Documentation

Definitions, calculations and limitations are described in the documentation:

| Topic | Link |
| :--- | :--- |
| **GUI & Calibration** | [User Guide → GUI Overview](https://francincz.github.io/ExQt/user-guide/gui-overview/) |
| **Nuclear ROI & Napari** | [User Guide → ROI & Napari](https://francincz.github.io/ExQt/user-guide/roi-and-napari/) |
| **Stack Aligner** | [User Guide → Stack Aligner](https://francincz.github.io/ExQt/user-guide/stack-aligner/) |
| **Dual-Scale Volume ($V_{\text{bio}}$ vs $V_{\text{gel}}$)** | [Methodology → Dual-Scale](https://francincz.github.io/ExQt/methodology/dual-scale/) |
| **Partition Coefficient ($K_{\text{part}}$)** | [Methodology → Partitioning](https://francincz.github.io/ExQt/methodology/partitioning/) |
| **Radial Profiling (layers)** | [Methodology → Radial Profiling](https://francincz.github.io/ExQt/methodology/radial-profiling/) |
| **3D Fractional Anisotropy (FA)** | [Methodology → 3D Anisotropy](https://francincz.github.io/ExQt/methodology/anisotropy/) |
| **Geometric Null Model** | [Methodology → Null Model](https://francincz.github.io/ExQt/methodology/null-model/) |
| **Quality Control** | [Quality Control → QC Pipeline](https://francincz.github.io/ExQt/quality-control/qc-pipeline/) |
| **Exclusion Rules** | [Quality Control → Exclusion Rules](https://francincz.github.io/ExQt/quality-control/exclusion-rules/) |
| **CSV Data Dictionary** | [Data Dictionary (All Columns)](https://francincz.github.io/ExQt/data-dictionary/) |
| **Reports** | [Reports & Visualizations](https://francincz.github.io/ExQt/reporting/) |

---

## Key Libraries & Dependencies

- **[Napari](https://napari.org/):** image viewer and ROI drawing.
- **[PySide6](https://doc.qt.io/qtforpython/):** user interface.
- **[scikit-image](https://scikit-image.org/):** connected components, object measurements, phase cross-correlation.
- **[SciPy](https://scipy.org/):** distance transform, image shifts, statistical test.
- **[NumPy](https://numpy.org/) & [pandas](https://pandas.pydata.org/):** calculations and tables.
- **[tifffile](https://github.com/cgohlke/tifffile):** reading and writing TIFF files and their calibration.
- **[matplotlib](https://matplotlib.org/) & [seaborn](https://seaborn.pydata.org/):** figures.
- **[openpyxl](https://openpyxl.readthedocs.io/):** Excel output.

---

## Acknowledgments

- **Z-stack aligner:** The XY drift correction (`stack_aligner.py`) is inspired by and adapted from the MATLAB **[3D-Aligner](https://github.com/suzukilabmcardle/3D-Aligner)** developed by the Suzuki Lab at the McArdle Laboratory for Cancer Research. It is re-implemented in Python.

---

## License
ExQt is open-source software released under the [MIT License](LICENSE).
