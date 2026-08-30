# ExQt: Expansion Microscopy Quantification of Nuclear Condensates

**ExQt** is an open-source desktop application for automated, reproducible 2D and 3D quantitative analysis of fluorescence microscopy images, designed with an emphasis on **Expansion Microscopy (ExM)** and cellular assemblies. It combines raw multi-dimensional TIFF intensity data with user-provided segmentation masks, applies physical spatial calibration and expansion scaling, and generates auditable object-level statistics, size distributions, intensity measurements, quality-control metrics, and publication-ready reports.

ExQt does **not** perform segmentation itself. Masks can be prepared in any preferred segmentation tool (e.g., [Labkit](https://imagej.net/plugins/labkit/), [ilastik](https://www.ilastik.org/), Fiji/ImageJ, Cellpose, or custom pipelines) capable of exporting standard TIFF masks.

---

## Main Features

- **Multi-Dimensional Processing:** Full 3D stack volume analysis, 2D Maximum Intensity Projections (MIP), or focus-ranked single-slice extraction.
- **Defensive Mask Support:** Accepts binary masks (`0/1`, `0/255` with automatic connected-component labeling) and instance label masks (positive integer IDs).
- **Physical ExM Calibration:** Full physical scaling using lateral XY pixel size, axial Z-step, and sample expansion factor ($\text{ExF}$) to convert pixel voxels to true biological volume ($V_{\text{bio}} = V_{\text{pixel}} / \text{ExF}^3$).
- **Morphometry & Intensity Metrics:** Calibrated volume/area, 3D equivalent spherical diameter, sphericity, total integrated intensity, mean/median intensity, and local background levels.
- **Interactive Size Distribution Preview:** Real-time visual histogram preview with draggable range sliders to inspect population distributions before running full batch analysis.
- **Flexible ROI Workflows:** Whole-image processing (**Auto-ROI**) or interactive 2D/3D polygon drawing via integrated [napari](https://napari.org/) viewers.
- **Reproducible Batch Pipeline:** Generates structured machine-readable CSVs, audit logs, and provenance metadata JSON files alongside multi-sheet formatted Excel reports.
- **Advanced Radial FA Profiling:** Optional 3D Shell–Middle–Core layer erosion using Euclidean distance transforms to quantify geometric fractional anisotropy ($FA$) gradients and topological branching (`z_topology_fail`).
- **Advanced 3D Thermodynamic Partitioning ($K_{\text{part}}$):** Optional single-cell dense-phase concentration ratio extraction ($C_{\text{dense}} / C_{\text{dilute}}$) and 2D biophysical phase diagrams.

---

## Installation

An isolated Conda environment is recommended.

### 1. Clone the repository
```bash
git clone https://github.com/FrancinCZ/ExQt.git
cd ExQt
```

### 2. Create and activate the environment
```bash
conda env create -f environment.yml
conda activate exqt-env
```

Alternatively, install dependencies into an existing environment via pip:
```bash
python -m pip install -r requirements.txt
```

### 3. Launch ExQt
```bash
python App.py
```

---

## Input Data Requirements

For each raw TIFF file, ExQt expects a corresponding segmentation mask in the same input folder with the suffix `_Mask.tif`:

```text
sample_01.tif
sample_01_Mask.tif
```

- **Binary masks:** Foreground is non-zero; distinct 3D connected components are segmented and numbered automatically.
- **Instance masks:** Each object has a unique positive integer label.

---

## Quick Start Guide

1. **Launch ExQt:** Run `python App.py`.
2. **Select Folders:** Choose the input directory with raw TIFFs and `_Mask.tif` files, and set an output directory.
3. **Choose Process Mode:** Select `3d`, `2d`, or `single_slice`.
4. **Set Calibration & Expansion:** Enter acquisition XY pixel size (nm), Z-step (nm), and the physical expansion factor (e.g. `1.0` for unexpanded, `4.0` for standard ExM, `10.0` for TREx).
5. **Adjust Size Range:** Use **Preview size distribution...** to inspect your data and set the lower/upper biological volume bounds.
6. **Configure ROI:** Choose **Auto-ROI** (entire field) or draw a custom manual cell boundary in the Napari viewer.
7. **Select Reports:** Check **Generate selected reports** and click **Configure...** to select desired Excel/CSV sheets and summary plots.
8. **Run Analysis:** Click **Start Processing** and confirm the parameter summary dialog.

---

## Size Settings & Noise Filtering

ExQt cleanly separates three different sizing parameters:
1. **Raw Noise Filter (voxels):** Early noise cutoff to discard single-pixel artifacts before feature extraction.
2. **Analyzed Biological Size Range ($\mu\text{m}^3$ or $\mu\text{m}^2$):** Calibrated biological boundaries used for primary statistics, summary tables, and main plots.
3. **Radial Layer Minimum Voxels:** Minimum layer volume required for valid *Shell–Middle–Core* geometric calculations.

---

## Advanced Analysis Modules

### Radial FA Profiling (Shell–Middle–Core)
Radial FA Profiling subdivides each 3D object into three concentric zones: **Shell** (outer), **Middle** (intermediate), and **Core** (inner). It calculates 3D fractional anisotropy ($FA$) per layer to evaluate geometric elongation and internal structural organization.

### 3D Thermodynamic Partitioning ($K_{\text{part}}$) & Phase Diagrams
ExQt can extract in situ single-cell partitioning coefficients:
$$K_{\text{part}} = \frac{\max(I_{\text{obj}} - I_{\text{offset}}, 0)}{\max(I_{\text{nuc}} - I_{\text{offset}}, 10^{-6})}$$
Combining $K_{\text{part}}$ with $FA$ generates 2D biophysical phase diagrams to assess molecular enrichment relative to surrounding nucleoplasm.  
 *For full mathematical derivations, physical background, and usage details, see [README_PARTITIONING.md](README_PARTITIONING.md).*

---

## Generated Outputs

A completed batch run produces:
- `*_Output_Batch_<mode>.csv` — Comprehensive machine-readable table with all morphological, intensity, and QC metrics.
- `*_Output_Batch_<mode>_metadata.json` — Exact provenance record of applied calibration, thresholds, and software settings.
- `*_Stats.xlsx` — Formatted multi-sheet workbook containing summaries, primary objects, excluded objects, and QC statistics.
- `*_3d_size_intensity_distribution.png` — Standard population overview plots (volume histograms, intensity correlations).
- `*_3d_partitioning_analysis.png` *(optional)* — *Size vs. $K_{\text{part}}$* and 2D Biophysical Phase Diagrams.
- `*_3d_Radial_FA_Profiling_Plots.png` *(optional)* — 4-panel radial layer progression and QC filtering funnel.

---

## Merging Completed Runs

Use **Tools → Merge existing runs...** to scan multiple completed run folders and combine compatible datasets into a unified summary spreadsheet (`Merged_Stats.xlsx`) and comparative multi-run overview figures (`Merged_Stats.png`).

---

## License & Acknowledgments

ExQt is open-source software released under the [MIT License](LICENSE). See [About.md](About.md) for author and institution details.

Built with Python using PySide6, napari, NumPy, pandas, SciPy, scikit-image, tifffile, openpyxl, matplotlib, and seaborn.
