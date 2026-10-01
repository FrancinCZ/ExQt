# Graphical User Interface (GUI) Overview

The ExQt interface (`app.py` or `App.py`) gives you complete control over your analysis, microscope calibrations, and visual inspection.

The window is divided into two main sections:
- **Left Panel:** All settings, calibration fields, size filters, and execution buttons.
- **Right Panel:** Integrated [Napari](https://napari.org/) 3D interactive viewer for inspecting optical stacks and drawing nuclear boundaries.

---

## The Left Control Panel

### 1. Data Selection
- **Input folder (`Browse`):** Select the folder containing your raw multi-channel `.tif` stacks and their corresponding segmentation masks (`*_Mask.tif`).
- **Output folder (`Browse...`):** Choose where results (`.csv`, `.xlsx`, figures) should be saved.
- **Save next to input:** When checked, results are written into the input folder instead of a separate output folder.

### 2. Processing Configuration
- **Process mode:**
  - **`3d` (Default):** Full volumetric analysis across all optical Z-sections. Required for 3D shape anisotropy and radial core–shell profiling.
  - **`2d`:** Maximum Intensity Projection (MIP) analysis.
  - **`single_slice`:** Automatically selects and analyzes only the single sharpest in-focus plane.
- **Expansion Factor:** Enter the linear physical expansion factor of your hydrogel (e.g., `1.0` for classical unexpanded confocal; ≈ `8` for our Magnify protocol imaged in distilled water – Magnify expands 3–4.5× in 1× PBS, 5–7× in 1:50 PBS and 8–11× in water, so always record the imaging medium). Prefer an ExF measured on nuclei over the macroscopic gel ExF, because nuclei can expand less than the gel. ExQt uses this to calculate true biological volumes ($V_{\text{bio}} = V_{\text{gel}} / \text{ExF}^3$).
- **Show Napari preview:** Keeps the 3D viewer active and synchronized with each processed image.
- **Auto-ROI:** When checked, treats the entire field of view as a single cell/nucleus (convenient for single-nucleus crops or high-throughput screens). When unchecked, ExQt pauses to let you outline the nucleus manually. With Auto-ROI, $K_{\text{part}}$ is not computed (NaN, `auto_roi_no_nucleoplasm`), because there is no nucleus outline to define the nucleoplasm.
- **Review each image:** An optional quality check — pauses after each image to let you visually approve the segmentation mask before moving to the next file.

### 3. Size Filtering & Interactive Preview
- **Min size / Max size** (µm³ in 3D, µm² in 2D; the unit is shown in the field): Sets the physical size boundaries for **Primary Condensates**. Objects smaller than Min (e.g., single-pixel noise) or larger than Max (e.g., giant unresolved clumps) are separated into the excluded dataset.
- **Preview sizes... (Button):** Opens an interactive histogram of your data *before* running the full batch. You can drag vertical sliders on the graph to visually pick the optimal Min and Max size cutoffs for your population.

### 4. Action Buttons
- **Start analysis (Blue button):** Starts the analysis of all paired files in the input folder.
- **Confirm ROI and Continue (Green button):** Appears when manual ROI drawing is active. Click this once you have drawn the nuclear boundary in Napari to proceed with feature extraction. Draw one shape per nucleus on one or more Z-slices. Shapes on different slices that overlap in XY are treated as the same nucleus and interpolated between its own slices (slices outside that range are excluded); a nucleus drawn on a single slice is extruded through the whole stack. Each nucleus gets its own `cell_id`. Shapes that overlap or touch on the same slice cannot be separated and are rejected with a message, so you can fix them. If nothing is drawn, ExQt asks before using the whole field of view (`roi_source = empty_fallback_fov`).
- **Approve & Next Image (Orange button):** Appears during Review Mode to approve the current segmentation and move to the next image.
- **Stop and Discard This Image (Red button):** Discards the currently displayed image and safely stops the batch (all previously approved images remain saved).

---

## Where to Set Pixel Size and Z-Step (Settings)

Microscope calibration parameters are accessed via the top menu bar:

**Navigation:** click `Settings` in the top menu bar; the settings window opens directly.

The **Settings** window contains the optical calibration and all analysis options:

| Setting | Default | Description |
| :--- | :---: | :--- |
| **Pixel Size XY (nm)** | `58.0` | Physical lateral pixel size recorded by the camera (unscaled). |
| **Z-step (nm)** | `250.0` | Axial optical distance between adjacent slices in the stack (3D mode only). |
| **Signal Channel** | `1` | 0-based channel index for the fluorescent condensate signal (e.g., green channel). |
| **DAPI Channel** | `0` | 0-based channel index for the nuclear counterstain (e.g., blue DAPI channel). |
| **Detector offset (ADU)** | `0.0` | Detector/camera offset subtracted from condensate and nucleoplasm intensity for $K_{\text{part}}$. From acquisition settings or a dark frame; `0` = no subtraction. Never estimated from the image. **From .lif…** reads a Leica .lif file: if every active detector is in photon-counting mode the offset is set to 0 and the file is recorded as its source. |
| **Noise filter (voxels)** | `5` | Early connected-component threshold to discard tiny single-pixel camera shot noise before calibration. |
| **Enable Radial FA Profiling (3D)** | Checked | Enables Mode A: concentric layer peeling (Core–Middle–Shell) and the Geometric Null Model. |
| **Min. voxels per layer** | `20` | Minimum number of voxels required in the core shell for valid Fractional Anisotropy calculation. |
| **Require Z-topology PASS** | Checked | Flags condensates that split into multiple fragments across Z-slices, excluding them from primary benchmarks. |
| **Dark mode** | Off | Dark or light appearance; remembered for the next start. |

> [!TIP]
> **Metadata Detection:** When you select an input folder with `Browse`, ExQt reads OME / ImageJ calibration from the TIFF headers (plain TIFF DPI such as the default 72 dpi is ignored; values outside XY 1–2000 nm or Z 1–5000 nm are reported and ignored). ExQt then asks what to apply:
> - **One calibration in all files:** *Yes* copies it into Settings, *No* keeps your settings. Either way, the Settings values are applied to every file.
> - **Different calibrations:** *Yes* applies each file's own TIFF calibration (Settings only for files without one), *No* applies Settings to all files.
>
> If the TIFF calibration ≈ Settings / ExF, ExQt warns that the files seem to be calibrated to the pre-expansion size (applying ExF again would make volumes ExF³× too small). The calibration actually applied to each file and its source are saved in `*_metadata.json` (`calibration_source`, `applied_calibration_by_file`).

---

## Top Menu: Tools

### 1. `Tools` → `Align Z-stacks...`
Opens the **Z-stack aligner** dialog to correct lateral (XY) drift between consecutive Z-slices (e.g. stage or thermal drift) before running condensate analysis. It does not correct chromatic aberration between channels:
- **Reference channel:** Choose the channel with the most continuous signal (or click **Auto-detect...** to let ExQt test inter-slice correlation automatically).
- **Max step per Z-plane (px):** Maximum allowed shift between adjacent slices (default 8.0 px).
- **Expand canvas:** Enlarges the field of view by the cumulative drift so no edge voxels are clipped; the added zero border is recorded in `<stem>_alignment_valid.tif` and excluded from statistics.
- **Align masks:** Shifts matching `*_Mask.tif` files with the exact same translation vectors.

### 2. `Tools` → `Merge existing runs...`
If you ran multiple batches across different days or folders, this tool scans the directory tree and merges all compatible output CSVs into a single unified Excel spreadsheet (`Merged_Stats.xlsx`) and comparison figure (`Merged_Stats.png`).

## Top Menu: Settings
One click opens the Settings window described above (calibration, channels, detector offset, Radial FA Profiling and dark mode).
