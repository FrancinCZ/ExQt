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
- **Save to the same folder as input:** When checked, creates an automated output directory directly inside your input folder.

### 2. Processing Configuration
- **Process mode:**
  - **`3d` (Default):** Full volumetric analysis across all optical Z-sections. Required for 3D shape anisotropy and radial core–shell profiling.
  - **`2d`:** Maximum Intensity Projection (MIP) analysis.
  - **`single_slice`:** Automatically selects and analyzes only the single sharpest in-focus plane.
- **Expansion Factor:** Enter the linear physical expansion factor of your hydrogel (e.g., `1.0` for classical unexpanded confocal, `4.0` for standard ExM, or any custom expansion value). ExQt uses this to calculate true biological volumes ($V_{\text{bio}} = V_{\text{gel}} / \text{ExF}^3$).
- **Show Napari preview:** Keeps the 3D viewer active and synchronized with each processed image.
- **Auto-ROI:** When checked, treats the entire field of view as a single cell/nucleus (convenient for single-nucleus crops or high-throughput screens). When unchecked, ExQt pauses to let you outline the nucleus manually.
- **Pause and review segmentations:** An optional quality check — pauses after each image to let you visually approve the segmentation mask before moving to the next file.

### 3. Size Filtering & Interactive Preview
- **Analyzed Size Min & Max ($\mu\text{m}^3$ or $\mu\text{m}^2$):** Sets the physical volume boundaries for **Primary Condensates**. Objects smaller than Min (e.g., single-pixel noise) or larger than Max (e.g., giant unresolved clumps) are separated into the excluded dataset.
- **Preview size distribution... (Button):** Opens an interactive histogram of your data *before* running the full batch. You can drag vertical sliders on the graph to visually pick the optimal Min and Max size cutoffs for your population.

### 4. Action Buttons
- **Start analysis (Blue button):** Begins automated batch processing across all paired files in the input folder.
- **Confirm ROI and Continue (Green button):** Appears when manual ROI drawing is active. Click this once you have drawn the nuclear boundary in Napari to proceed with feature extraction.
- **Approve & Next Image (Orange button):** Appears during Review Mode to approve the current segmentation and move to the next image.
- **Stop and Discard This Image (Red button):** Discards the currently displayed image and safely stops the batch (all previously approved images remain saved).

---

## Where to Set Pixel Size and Z-Step (Advanced Settings)

Microscope calibration parameters are accessed via the top menu bar:

**Navigation:** `Settings` → `Advanced...` (from the top menu bar)

Opening **Advanced Settings** reveals the core optical calibration dialog:

| Setting | Default | Description |
| :--- | :---: | :--- |
| **Pixel Size XY (nm)** | `58.0` | Physical lateral pixel size recorded by the camera (unscaled). |
| **Z-step (nm)** | `250.0` | Axial optical distance between adjacent slices in the stack (3D mode only). |
| **Signal Channel** | `1` | 0-based channel index for the fluorescent condensate signal (e.g., green channel). |
| **DAPI Channel** | `0` | 0-based channel index for the nuclear counterstain (e.g., blue DAPI channel). |
| **Raw noise filter (voxels)** | `5` | Early connected-component threshold to discard tiny single-pixel camera shot noise before calibration. |
| **Enable Radial FA Profiling (3D)** | Checked | Enables Mode A: concentric layer peeling (Core–Middle–Shell) and the Geometric Null Model. |
| **Min. voxels per layer** | `20` | Minimum number of voxels required in the core shell for valid Fractional Anisotropy calculation. |
| **Require Z-topology PASS** | Checked | Flags condensates that split into multiple fragments across Z-slices, excluding them from primary benchmarks. |

> [!TIP]
> **Automatic Metadata Detection:** When you select an input folder with `Browse`, ExQt automatically reads the TIFF image headers. If valid physical pixel sizes are detected, ExQt auto-fills XY and Z calibration and confirms it in the bottom status bar!

---

## Top Menu: Tools

### 1. `Tools` → `Align Z-stacks...`
Opens the **Multi-Channel Stack Aligner** dialog to correct for thermal stage drift or chromatic aberration before running condensate analysis:
- **Reference channel:** Choose the channel with the most continuous signal (or click **Auto-detect...** to let ExQt test inter-slice correlation automatically).
- **Max step per Z-plane (px):** Maximum allowed shift between adjacent slices (default 8.0 px).
- **Expand canvas:** Automatically enlarges the field of view by the cumulative drift so no edge voxels are clipped.
- **Align masks:** Shifts matching `*_Mask.tif` files with the exact same translation vectors.

### 2. `Tools` → `Merge existing runs...`
If you ran multiple batches across different days or folders, this tool scans the directory tree and merges all compatible output CSVs into a single unified Excel spreadsheet (`Merged_Stats.xlsx`) and comparison figure (`Merged_Stats.png`).

### 3. `Settings` → `Dark Mode`
Toggles between clean dark mode and light mode stylesheets for comfortable viewing in dim microscopy rooms.