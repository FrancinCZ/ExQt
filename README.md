# ExQt: Expansion Microscopy Quantification Tool

**ExQt** (*Expansion Microscopy Quantification Tool*) is a desktop application for measuring biomolecular condensates in 3D confocal and Expansion Microscopy (ExM) images. It does not segment images itself: it takes your existing segmentation masks (e.g. from Labkit or ilastik) and measures size, enrichment (K_part), shape (FA) and the core–shell profile of every object, with quality control and a full record of how each result was produced.

 **[📖 Read the Documentation](https://francincz.github.io/ExQt/)**

---

## Quick Installation

ExQt needs **Python 3.10 or newer** (tested with Python 3.14) on Windows, Linux or macOS. Conda is not required.

```bash
# 1. Clone the repository
git clone https://github.com/FrancinCZ/ExQt.git
cd ExQt

# 2. Install dependencies (a virtual environment is recommended)
pip install -r requirements.txt

# 3. Launch the application
python App.py
```

Alternatively, with Conda: `conda env create -f environment.yml`, then `conda activate exqt-env` and `python App.py`.

To check an installation, run the test suite: `python -m pytest tests`.

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

For full mathematical definitions, algorithm walkthroughs, and data interpretation, visit our documentation:

| Topic | Link |
| :--- | :--- |
| **GUI & Calibration** | [User Guide → GUI Overview](https://francincz.github.io/ExQt/user-guide/gui-overview/) |
| **Nuclear ROI & Napari** | [User Guide → ROI & Napari](https://francincz.github.io/ExQt/user-guide/roi-and-napari/) |
| **Stack Aligner** | [User Guide → Stack Aligner](https://francincz.github.io/ExQt/user-guide/stack-aligner/) |
| **Dual-Scale Volume ($V_{\text{bio}}$ vs $V_{\text{gel}}$)** | [Methodology → Dual-Scale](https://francincz.github.io/ExQt/methodology/dual-scale/) |
| **Partition Coefficient ($K_{\text{part}}$)** | [Methodology → Partitioning](https://francincz.github.io/ExQt/methodology/partitioning/) |
| **Core–Shell Radial Profiling** | [Methodology → Radial Profiling](https://francincz.github.io/ExQt/methodology/radial-profiling/) |
| **3D Fractional Anisotropy (FA)** | [Methodology → 3D Anisotropy](https://francincz.github.io/ExQt/methodology/anisotropy/) |
| **Geometric Null Model** | [Methodology → Null Model](https://francincz.github.io/ExQt/methodology/null-model/) |
| **Quality Control Tiers** | [Quality Control → QC Pipeline](https://francincz.github.io/ExQt/quality-control/qc-pipeline/) |
| **Exclusion Rules** | [Quality Control → Exclusion Rules](https://francincz.github.io/ExQt/quality-control/exclusion-rules/) |
| **CSV Data Dictionary** | [Data Dictionary (All Columns)](https://francincz.github.io/ExQt/data-dictionary/) |
| **Plot Interpretation** | [Reports & Visualizations](https://francincz.github.io/ExQt/reporting/) |

---

## Key Libraries & Dependencies

ExQt is built in Python, integrating core scientific and bio-imaging libraries:

- **[Napari](https://napari.org/):** Multi-dimensional 3D viewer embedded directly in the main window for interactive stack navigation and manual nuclear ROI delineation.
- **[PySide6 (Qt for Python)](https://doc.qt.io/qtforpython/):** Native desktop graphical user interface framework powering all dialogs, interactive controls, and background batch worker threads.
- **[scikit-image](https://scikit-image.org/):** 3D connected-component labeling, volumetric feature extraction (`measure.regionprops`), and Fourier phase cross-correlation.
- **[SciPy](https://scipy.org/):** Exact 3D Euclidean Distance Transform (`scipy.ndimage.distance_transform_edt`) for Core–Middle–Shell radial peeling, and subpixel n-dimensional image shifts.
- **[NumPy](https://numpy.org/) & [pandas](https://pandas.pydata.org/):** Spatial covariance tensor decomposition (Fractional Anisotropy), geometric null model simulations, and high-throughput tabular data management.
- **[tifffile](https://github.com/cgohlke/tifffile):** Robust reading and writing of multi-channel 3D OME-TIFF files while preserving microscope physical calibration tags.
- **[matplotlib](https://matplotlib.org/) & [seaborn](https://seaborn.pydata.org/):** Summary figures (300 DPI) and the interactive size preview histogram.
- **[openpyxl](https://openpyxl.readthedocs.io/):** Multi-sheet Excel workbook (`_Detailed_Stats.xlsx`).

---

## Acknowledgments & Inspiration

- **Z-stack aligner:** The XY drift correction (`stack_aligner.py`) is inspired by and adapted from the MATLAB **[3D-Aligner](https://github.com/suzukilabmcardle/3D-Aligner)** developed by the Suzuki Lab at the McArdle Laboratory for Cancer Research. It has been re-implemented in Python with sub-pixel phase cross-correlation, Hanning-window preconditioning, canvas expansion (the added border is marked and excluded from statistics) and per-step quality checks.
- **Expansion Microscopy Community:** Designed to support quantitative biophysical analysis of biomolecular condensates and nuclear assemblies across classical confocal and high-expansion microscopy (ExM).

---

## License
ExQt is open-source software released under the [MIT License](LICENSE).
