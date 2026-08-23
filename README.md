# ExQt: Expansion Microscopy Quantification of Nuclear Condensates

ExQt is an open-source desktop application for quantitative analysis of 2D and 3D microscopy images, with an emphasis on expansion microscopy (ExM) and nuclear condensates. It combines raw TIFF intensity data with a user-provided segmentation mask, applies physical calibration and the expansion factor, and produces auditable object-level measurements, quality-control results, reports, and plots.

ExQt does **not** perform segmentation itself. Masks can be prepared in software such as [Labkit](https://imagej.net/plugins/labkit/), [ilastik](https://www.ilastik.org/), Fiji, or another segmentation tool capable of exporting TIFF masks.

## Main features

- Analysis of full 3D stacks, maximum-intensity projections, or a selected single slice.
- Binary-mask and instance-label-mask support with defensive validation.
- Physical calibration using XY pixel size, Z-step, and the sample expansion factor.
- Object coordinates, calibrated area or volume, equivalent diameter, and intensity measurements.
- Optional manual ROI selection and napari preview.
- Reproducible batch processing with a machine-readable CSV and metadata JSON.
- Configurable human-readable Excel, CSV, and plot reports.
- Optional **Rezim A** analysis of radial Shell–Middle–Core geometric fractional anisotropy (FA).
- Rezim A quality control for image/ROI edges, incomplete radial layers, and persistent splitting across Z.
- Optional merging of compatible completed runs into one Excel workbook and overview plot.

> [!IMPORTANT]
> Rezim A measures a geometric pattern: how shape anisotropy changes from the outer shell toward the core. It does not directly measure liquidity, viscosity, stiffness, molecular exchange, or phase separation. A positive Shell-to-Core trend may be biologically interesting, but it is not by itself proof of LLPS.

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

The environment file already installs the application dependencies. Alternatively, dependencies can be installed into an existing environment:

```bash
python -m pip install -r requirements.txt
```

### 3. Start ExQt

```bash
python App.py
```

## Input data

For every raw TIFF, ExQt expects a mask in the same input directory with the suffix `_Mask`:

```text
sample.tif
sample_Mask.tif
```

The raw image and mask must describe the same spatial field and have compatible dimensions.

### Supported masks

- **Binary masks:** background is zero and all positive pixels/voxels are foreground. Conventional `0/1` and `0/255` masks are accepted and connected components are labeled automatically.
- **Instance masks:** background is zero and every object has its own positive integer ID.

Non-integer values, negative labels, NaN/Inf values, and reused instance IDs in disconnected regions are rejected instead of being silently interpreted. For automatic connected-component labeling, export a binary mask. For instance segmentation, use one unique integer ID per 3D object.

## Quick start

1. Start ExQt with `python App.py`.
2. Select the folder containing the raw TIFF files and matching `_Mask.tif` files.
3. Select an output folder, or enable **Save to the same folder as input**.
4. Choose the process mode: `3d`, `2d`, or `single_slice`.
5. Set the expansion factor and the analyzed biological-size range.
6. Open **Settings → Advanced...** and verify calibration, channels, the raw noise filter, and optional Rezim A settings.
7. Choose whether to use the whole image (**Auto-ROI**) or draw a manual ROI.
8. Optionally enable napari preview or pause-and-review.
9. Enable **Generate selected reports** and use **Configure...** to select the desired outputs.
10. Start the analysis and confirm the calibration summary shown before processing.

## Calibration and expansion factor

ExQt distinguishes acquisition calibration from effective biological sampling:

- **Pixel Size XY (nm):** lateral sampling of the acquired image.
- **Z-step (nm):** distance between acquired Z planes; used only for 3D analysis.
- **Expansion Factor:** physical sample expansion ratio. Use `1.0` for non-expanded data.

For example, an acquisition pixel size of 58 nm with a 4× expansion factor corresponds to an effective biological XY sampling of 14.5 nm.

ExQt attempts to read usable TIFF metadata and can recommend detected calibration values. Detection is advisory: the values actually used are shown for explicit confirmation before every run. The applied values and detected metadata evidence are recorded in the run metadata JSON.

## Process modes

- **`3d`:** analyzes connected components throughout the complete volume. Z-step and voxel volume are used. Rezim A is available only in this mode.
- **`2d`:** analyzes a maximum-intensity projection and the projected mask.
- **`single_slice`:** extracts one selected/focus-ranked Z plane and analyzes it as a 2D image.

## Size settings: three different purposes

These settings are intentionally separate and should not be interpreted as interchangeable:

1. **Raw noise filter (pixels/voxels)** — an early connected-component noise floor. It removes tiny fragments before the main analysis.
2. **Analyzed biological-size range (µm² or µm³)** — the calibrated lower and upper bounds used by primary statistics, clean report tables, and plots.
3. **Rezim A minimum voxels per FA layer** — the minimum amount of data required for valid Shell, Middle, and Core FA calculations. The same value also sets the minimum valid core size.

The raw machine/audit CSV may therefore contain more components than the final size-eligible or primary QC-valid sets. This is expected and allows exclusions to remain auditable.

## Rezim A: Shell–Middle–Core FA

Rezim A divides each eligible 3D object into three normalized radial regions using a distance-transform-based scheme:

- **Shell:** outer radial layer
- **Middle:** intermediate radial layer
- **Core:** inner radial layer

Geometric fractional anisotropy is calculated independently for each layer using physically scaled coordinates. The main derived values are:

- `Middle − Shell`
- `Core − Middle`
- `Delta A = Core FA − Shell FA`

A positive Delta A means that the object's core is geometrically more anisotropic than its shell. The per-object paired plots show whether this pattern is consistent across condensates rather than being driven only by a pooled average.

### Rezim A quality control

Primary Rezim A plots use only objects that are:

- inside the selected biological-size range;
- complete, with valid Shell, Middle, and Core measurements;
- not touching a disallowed image or ROI edge; and
- accepted by the configured Z-topology policy.

Z-topology assessment ignores insignificant detached fragments and classifies persistent substantial splitting across slices as:

- **PASS:** suitable for primary analysis;
- **REVIEW:** uncertain topology, retained for audit but excluded from the primary set when strict Z-topology is enabled;
- **FAIL:** strong evidence of persistent splitting or clutter, retained for audit but excluded from the primary set.

The option **Require Z-topology PASS for primary comparison** is deliberately conservative: uncertain objects remain visible in the outputs, but do not silently influence the primary Shell–Middle–Core comparison.

## Generated files

The source batch CSV and its metadata JSON are always the reproducibility backbone of a run. Additional human-facing outputs are optional.

### Always generated

- `*_Output_Batch_<mode>.csv` — complete machine/audit table with measurements and QC fields.
- `*_Output_Batch_<mode>_metadata.json` — applied calibration, channels, report settings, Rezim A settings, QC policy, and provenance.

### Optional report outputs

The **Configure...** dialog can enable or disable:

- a readable Excel report;
- a compact primary QC-valid CSV;
- a QC-excluded CSV with exclusion reasons;
- an additional full raw audit CSV;
- standard size/intensity plots;
- Rezim A Shell–Middle–Core and QC plots.

The Excel report separates summary, primary objects, excluded objects, raw measurements, QC policy, and explanatory content into dedicated sheets instead of presenting one unstructured table.

## Merging completed runs

Use **Tools → Merge existing runs...** to recursively scan a folder containing completed ExQt run CSV files.

Each source batch CSV must have its matching metadata file beside it:

```text
sample_Output_Batch_3d.csv
sample_Output_Batch_3d_metadata.json
```

Before merging, ExQt compares the QC-policy fingerprints. Runs with incompatible analysis or QC settings are rejected rather than silently pooled. Acquisition pixel size and Z-step may differ and remain recorded per run, but the settings that define the analyzed population and QC policy must be compatible.

The merge produces:

- `Merged_Stats.xlsx` — per-run summaries, primary objects, QC-excluded objects, QC-policy comparison, and optionally the complete raw audit table;
- `Merged_Stats.png` — per-run median radial FA, median Delta A, QC acceptance, and descriptive object-level Delta A distributions.

Pooled condensates are descriptive observations, not independent biological replicates. For biological inference, the experimental unit should normally remain the independently acquired image, cell, or sample, depending on the study design.

## Interpreting the main plots

- **Per-object Shell → Middle → Core FA:** shows paired radial changes for every accepted object.
- **Layer-to-layer FA changes:** separates the Shell-to-Middle and Middle-to-Core contributions.
- **Delta A distribution:** summarizes `Core FA − Shell FA`; values above zero indicate greater core anisotropy.
- **QC funnel:** shows the transition from all segmented components to size-eligible, complete-FA, and primary QC-valid objects.
- **Standard plots:** describe calibrated size, intensity, density, and the number of condensates per selected ROI/cell.

Interpretation should consider the number of independent runs, run-to-run consistency, segmentation quality, acquisition calibration, and QC acceptance. A consistent trend across several independently acquired images is stronger evidence than many condensates from a single image.

## Reproducibility recommendations

- Keep every batch CSV together with its metadata JSON.
- Confirm XY pixel size, Z-step, and expansion factor for every acquisition.
- Use the same segmentation and QC policy for runs intended for comparison.
- Do not choose size or QC thresholds after inspecting the desired biological result.
- Preserve QC-excluded objects and reasons instead of deleting them.
- Treat Rezim A as a geometric pattern analysis unless validated against an independent reference or biological control.

## License and acknowledgments

ExQt is released under the [MIT License](LICENSE). See [About.md](About.md) for project and author information.

The application relies on open-source projects including PySide6, napari, NumPy, pandas, SciPy, scikit-image, tifffile, openpyxl, matplotlib, and seaborn.


