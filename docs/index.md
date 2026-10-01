# ExQt: Expansion Microscopy Quantification Tool

**ExQt** (*Expansion Microscopy Quantification Tool*) is an open-source desktop application for measuring objects in 3D fluorescence images of the cell nucleus, from confocal immunofluorescence (IF) or Expansion Microscopy (ExM). It works on existing segmentation masks; it does not segment images itself.

---

## What ExQt Measures

For every object in the mask:

- **Size** in the image (gel) scale and, divided by the expansion factor, in the estimated pre-expansion scale ($V_{\text{bio}} = V_{\text{gel}}/\text{ExF}^3$).
- **Enrichment** $K_{\text{part}}$: mean intensity of the object relative to the surrounding nucleoplasm.
- **Shape**: 3D Fractional Anisotropy (FA) of the mask.
- **Layers** (optional): FA and mean intensity of three concentric layers, and the FA difference compared with a homogeneous ellipsoid of the same shape.
- **Quality-control flags** and a record of the settings, code version and input files of each run.

All lengths use the physical voxel size, so the coarser Z sampling of confocal stacks is taken into account.

---

## Pipeline

`Batch.process_condensates` processes one raw/mask TIFF pair; the GUI worker in `App.py` runs it for every file and writes the CSV and `*_metadata.json`.

1. Load the stack and mask; select the signal channel; alignment QC (`Batch.validate_alignment_qc`).
2. ROI: Auto-ROI (whole field of view) or manual nucleus ROI, interpolated per nucleus (`Batch._interpolate_or_extrude_roi`).
3. Connected components of the mask inside the ROI; raw noise floor `min_voxels`.
4. Per object: volumes, $K_{\text{part}}$, edge flags and, optionally, layer metrics with the ellipsoid reference.
5. Reports: primary set by `postprocessing._prepare_reporting_frames`, Excel, plots.

---

## Installation

ExQt needs **Python 3.10 or newer** (tested with Python 3.14) and runs on Windows, Linux and macOS. Open a terminal (Command Prompt on Windows) and run:

```bash
git clone https://github.com/FrancinCZ/ExQt.git
cd ExQt
python -m pip install -r requirements.txt
python App.py
```

Without Git: on GitHub click **Code → Download ZIP**, unzip it, open a terminal in the unzipped folder and run the last two commands.

With Conda instead of pip: `conda env create -f environment.yml`, then `conda activate exqt-env` and `python App.py`.
