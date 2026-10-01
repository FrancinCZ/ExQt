# ExQt: Expansion Microscopy Quantification Tool

**ExQt** (*Expansion Microscopy Quantification Tool*) is an open-source desktop application for measuring biomolecular condensates in 3D images. It works on existing segmentation masks (it does not segment images itself) and measures size, enrichment, shape and the core–shell profile of every object.

It is designed to work across both **classical confocal immunofluorescence (IF)** and **Expansion Microscopy (ExM)**, providing a unified, auditable pipeline to measure the size, shape, enrichment and internal intensity profile of condensates in the cell nucleus.

---

## Why ExQt?

In diffraction-limited confocal microscopy, structures closer than the Abbe resolution limit can overlap and appear as a single object. Expansion Microscopy addresses this by physically magnifying the specimen in a swellable polymer hydrogel, allowing standard confocal hardware to resolve finer structural detail.

Analyzing ExM data introduces specific computational requirements:

- **True 3D processing** — maintaining connectivity across anisotropic voxel grids ($\Delta z > \Delta x, \Delta y$).
- **Dual-scale reporting** — preserving both the physical hydrogel scale ($V_{\text{gel}}$) and the estimated biological pre-expansion scale ($V_{\text{bio}} = V_{\text{gel}}/\text{ExF}^3$).
- **Internal profiling** — quantifying concentration gradients across concentric layers (Core–Middle–Shell), rather than treating condensates as uniform blobs.
- **Shape analysis** — comparing the FA of inner and outer layers with a homogeneous-ellipsoid null model of the same shape.

---

## Pipeline

`Batch.process_condensates` processes one raw/mask TIFF pair; the GUI worker in `App.py` runs it for every file and writes the CSV and `*_metadata.json`.

1. Load the stack and mask; select the signal channel; alignment QC (`Batch.validate_alignment_qc`).
2. ROI: Auto-ROI (whole field of view) or manual nucleus ROI, interpolated per nucleus (`Batch._interpolate_or_extrude_roi`).
3. Connected components of the mask inside the ROI; raw noise floor `min_voxels`.
4. Per object: gel and biological volume ($V_{\text{bio}} = V_{\text{gel}}/\text{ExF}^3$), $K_{\text{part}}$ with an explicit detector offset, edge/ROI-edge flags, optional Radial FA Profiling with the geometric null model.
5. Reports: primary set by `postprocessing._prepare_reporting_frames`, Excel, plots.


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

---
