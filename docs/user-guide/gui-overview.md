# Graphical User Interface (GUI) Overview

The ExQt window (`App.py`) contains the analysis settings on the left and a Napari viewer on the right.

The window is divided into two main sections:
- **Left panel:** input and output folders, analysis options, size range and buttons.
- **Right panel:** [Napari](https://napari.org/) viewer for inspecting stacks and drawing nuclear ROIs.

---

## The Left Control Panel

### 1. Data Selection
- **Input folder (`Browse`):** folder with the raw `.tif` stacks and their masks (`*_Mask.tif`).
- **Output folder (`Browse...`):** where results are saved.
- **Save next to input:** When checked, results are written into the input folder instead of a separate output folder.

### 2. Processing Configuration
- **Process mode:**
  - **`3d` (default):** analysis of the whole stack; required for FA and layer profiling.
  - **`2d`:** analysis of the maximum intensity projection.
  - **`single_slice`:** analysis of the sharpest slice only.
- **Expansion factor:** linear expansion factor (`1.0` for unexpanded samples). It depends on the protocol and the imaging medium; use a value measured for your samples. ExQt uses it for the estimated pre-expansion sizes ($V_{\text{bio}} = V_{\text{gel}} / \text{ExF}^3$).
- **Show Napari preview:** shows each processed image in the viewer.
- **Auto-ROI:** uses the whole field of view instead of a drawn nucleus ROI. $K_{\text{part}}$ is then not computed, because no nucleoplasm is defined. When unchecked, ExQt pauses for you to outline the nuclei.
- **Review each image:** pauses after each image so you can approve or discard it.

### 3. Size Filtering & Interactive Preview
- **Min size / Max size** (µm³ in 3D, µm² in 2D): size range of the primary objects; objects outside it are not used in the primary statistics.
- **Preview sizes...:** histogram of object sizes before the analysis; the range can be set by dragging.

### 4. Action Buttons
- **Start analysis:** analyses all paired files in the input folder.
- **Preview ROI / Use saved ROI / Confirm ROI and Continue:** shown while drawing the ROI, with step-by-step help above them (see *Nuclear ROI & Napari*).
- **Approve & Next Image:** in review mode, keeps the current image and continues.
- **Stop and Discard This Image:** discards the current image and stops; images approved earlier are saved.

---

## Where to Set Pixel Size and Z-Step (Settings)

Click **Settings** in the top menu bar:

| Setting | Default | Description |
| :--- | :---: | :--- |
| **Pixel Size XY (nm)** | `58.0` | Pixel size of the image, before expansion correction. |
| **Z-step (nm)** | `250.0` | Distance between slices (3D only). |
| **Signal Channel** | `1` | 0-based index of the measured channel. |
| **DAPI Channel** | `0` | 0-based index of the nuclear stain channel. |
| **Detector offset (ADU)** | `0.0` | Subtracted from intensities for $K_{\text{part}}$; `0` = no subtraction. **From .lif…** reads it from a Leica .lif file (photon counting → 0). |
| **Noise filter (voxels)** | `5` | Objects with fewer voxels are discarded before any measurement. |
| **Enable Radial FA Profiling (3D)** | Off | Layer metrics (Core–Middle–Shell) and the geometric null model. |
| **Min. voxels per layer** | `20` | Minimum voxels for a valid layer FA and for the core. |
| **Require Z-topology PASS** | Checked | Objects split along Z are not primary. |
| **QC: Keep objects cut by the first/last Z slice** | Unchecked | For stacks that do not cover the whole nucleus in Z: objects touching only the first/last slice stay primary ($K_{\text{part}}$ usable, size/FA underestimated). Runs with and without it are not merged. |
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
Corrects lateral (XY) drift between consecutive Z-slices and writes aligned copies to a new folder (see *Stack Aligner*):
- **Reference channel:** channel used to estimate the drift (**Auto-detect...** suggests one).
- **Max step per Z-plane (px):** largest allowed shift between neighbouring slices (default 8 px).
- **Expand canvas:** enlarges the image so nothing is cropped; the added border is excluded from statistics.
- **Align masks:** shifts the masks with the same translations.

### 2. `Tools` → `Merge existing runs...`
Merges compatible run CSVs from a folder tree into one Excel file (`Merged_Stats.xlsx`) and figure (`Merged_Stats.png`); runs with different QC settings are not merged.

### 3. `Tools` → `Compare runs...`
Draws groups of runs side by side without merging them, e.g. ExM vs IF. Add one folder per group (`Add folder...`); every run found in that folder and its subfolders belongs to the group, and the group name is editable. Runs inside one group must share the QC settings. Settings that differ between groups (expansion factor, size range, `allow_z_edge` …) are allowed, listed in red under the figure and on the `Settings_Differences` sheet. The **Nuclei** list shows every nucleus with QC-passed objects (e.g. `ExM img001 C1`, with its object count); untick nuclei to leave them out. When two nuclei would get the same name, the full run path is added (e.g. `ExM 2026.08.06/img001/Result/img001 C1`). **Each nucleus as its own column (coloured by group)** draws one box per picked nucleus instead of one per group. Output: `<name>.png` (the six panels of the standard plots on their fixed axes: size, K_part, size vs. K_part, FA, size vs. FA, objects per nucleus) and `<name>.xlsx`: `Per_Group` (runs, nuclei, QC funnel: raw, size-eligible and QC-passed objects with the acceptance rate, nucleus medians), `Per_Nucleus`, `Settings_Differences`, and the combined object tables of all compared runs as in a run's `*_Detailed_Stats.xlsx` (`Primary_Objects`, `QC_Excluded`, `All_Objects_Raw`, each with `group` and `run_id`; with picked nuclei only those nuclei). The K_part axis of the figure is fitted to the plotted values (see [Reports](../reporting.md)). When the runs have particle classes (Radial FA Profiling), Compare runs also saves `<name>_Classes.png`: the % of QC-passed objects in each class, one bar per group (or per nucleus with **Each nucleus as its own column**). A bar pools all objects of its column; white dots are the single nuclei, so a nucleus with many objects does not hide the spread. The numbers are on the `Class_Percentages` and `Class_Per_Nucleus` sheets. Compare groups by nuclei; size is only comparable when the analyzed size ranges match.

## Top Menu: Settings
Opens the Settings window described above.
