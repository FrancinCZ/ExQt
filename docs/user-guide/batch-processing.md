# Batch Processing

The batch analysis (`Batch.py`) processes every raw/mask pair in a folder one after another. With a manual ROI it pauses at each image so you can outline the nuclei; with Auto-ROI it runs without interaction.

---

## Directory Organization & File Pairing Convention

ExQt pairs each image with its mask by file name:

```
Experiment_Folder/
├── cell01.tif              <── image stack
├── cell01_Mask.tif         <── mask (binary or labels)
├── cell02.tif
├── cell02_Mask.tif
└── Result/                 <── Output folder of one run
    ├── cell01_ROI.tif      <── ROI used in this run (one label per nucleus)
    ├── cell02_ROI.tif
    ├── Experiment_Folder_Output_Batch_3d.csv
    └── Experiment_Folder_Output_Batch_3d_metadata.json
```

Each run saves its ROIs in its own output folder, so a new run never overwrites the ROI of an earlier one; the input folder is not modified. Leftover `*_ROI.tif` files in the input folder are ignored.


---

## Batch Lifecycle per Image Pair

For each pair:

```
[1/5] Load stack and mask
       │
[2/5] ROI (Auto-ROI or drawn), saved as *_ROI.tif
       │
[3/5] Objects of the mask inside the ROI
       │
[4/5] Measurements (volumes, K_part, FA, layers)
       │
[5/5] CSV and metadata.json
```

---

## Run Metadata (`*_Output_Batch_<mode>_metadata.json`)

Written by the analysis worker in `App.py` next to the CSV. Before comparing numbers between two runs, first check that `provenance` and `files.input_files` match.

| Key | Content |
| :--- | :--- |
| `provenance.git_commit`, `git_dirty` | Code revision; `git_dirty = true` means the tracked code differed from that commit. |
| `provenance.source_sha256` | SHA-256 of every analysis module – identifies the code even without git. |
| `provenance.python`, `platform`, `packages` | Interpreter, OS and library versions (numpy, scipy, scikit-image, pandas, tifffile, …). |
| `files.input_files` | Name, size and SHA-256 of each analysed TIFF, its mask and alignment sidecars (`_drift.csv`, `_alignment_valid.tif`). |
| `files.analysed_files`, `alignment_review_files`, `excluded_files` | What entered the CSV, which stacks were alignment REVIEW, and every excluded file with its reason. |
| `parameters` | Mode, ExF, `auto_roi`, raw noise floor, size range, GUI calibration, `calibration_source`, `applied_calibration_by_file`. |
| `partitioning` | K offset method and value, Auto-ROI and padding policy. |
| `classification` | Gradient test and significance level. |
| `mode_a` | Radial FA Profiling settings and Z-topology policy. |

The Excel `QC_Policy` sheet derives a fingerprint from the settings that change the meaning of the results (including `auto_roi`, K offset method and gradient test); `Merge existing runs` refuses to pool runs whose fingerprints differ. Offset value, calibration and code revision are shown but not part of the fingerprint, because they may legitimately differ between correct runs.
