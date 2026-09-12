# ExQt: Expansion & Confocal Quantitative Toolkit

**ExQt** is an automated 3D image-analysis platform for the segmentation, biophysical quantification, and morphological profiling of biomolecular condensates in confocal and Expansion Microscopy (ExM).

 **[📖 Read the Documentation](https://francincz.github.io/ExQt/)**

---

## Quick Installation

ExQt requires **Python 3.10–3.14** on Windows, Linux, or macOS.

```bash
#1. Clone the repository
git clone https://github.com/FrancinCZ/ExQt.git
cd ExQt

#2. Install dependencies
pip install -r requirements.txt

#3. Launch the application
python App.py
```

---

## User Guide: 4 Steps to Run an Analysis

### Step 1: Prepare Your Input Files
Place your raw 3D fluorescence stacks and their corresponding segmentation masks in the same folder. ExQt pairs them automatically by filename:

```text
My_Experiment/
├── cell_01.tif          <-- Raw 3D fluorescence image
└── cell_01_Mask.tif     <-- Binary or labeled mask (from ilastik, Labkit, etc.)
```

### Step 2: Set Calibration in the GUI
1. Select your **Input Folder** and **Output Folder**.
2. Enter your microscope calibration:
   - **Pixel Size XY (nm):** e.g., `65.0`
   - **Z-Step (nm):** e.g., `200.0`
   - **Expansion Factor:** `1.0` for classical confocal, or your physical gel expansion factor for ExM.
3. Set your **Size Filter** range (minimum and maximum biological volume in $\mu\text{m}^3$) to exclude single-pixel noise and giant clumps.

### Step 3: Define Nuclear ROI (Interactive Napari)
Condensates must be analyzed inside the nucleus to obtain accurate background and partitioning measurements:
- When prompted, draw the nuclear boundary in the integrated **Napari 3D viewer**.
- ExQt automatically extrudes your 2D outline into a 3D cylindrical volume across all optical sections.

### Step 4: Run & Inspect Results
Click **Start Processing**. ExQt will analyze all matched pairs sequentially and generate outputs in your results folder:

- **`Result_Primary_Condensates.csv`** — The clean, publication-ready dataset (passed all QC criteria).
- **`Result_Excluded_Condensates.csv`** — Filtered-out objects with explicit rejection reasons.
- **`_Output_Batch_3d.csv`** — Complete raw dataset with all 30+ physical and morphological parameters.
- **`_Stats.xlsx`** — Formatted Excel workbook.
- **`*_3d_partitioning_analysis.png`** — Automated 4-panel publication dashboard.

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

## License
ExQt is open-source software released under the [MIT License](LICENSE).
