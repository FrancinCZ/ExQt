# High-Throughput Batch Processing Engine

The ExQt batch processing engine (`Batch.py`) automates the feature extraction lifecycle across entire experimental directories containing dozens of 3D multi-gigabyte acquisitions without requiring human supervision.

---

## Directory Organization & File Pairing Convention

ExQt matches raw fluorescent intensity images with their corresponding binary segmentation masks using a strict suffix convention:

```
Experiment_Folder/
├── cell01.tif              <── Raw 3D fluorescent intensity stack
├── cell01_Mask.tif         <── 3D binary segmentation mask (labels or binary)
├── cell01_ROI.tif          <── [Auto-generated] Nuclear ROI mask saved by ExQt
├── cell02.tif
├── cell02_Mask.tif
└── cell02_ROI.tif
```

When a batch run is executed, ExQt scans the directory, matches each raw `.tif` with its paired `*_Mask.tif`, validates dimensions, and processes each pair sequentially.

---

## Batch Lifecycle per Image Pair

For each matched pair, ExQt executes the following deterministic sequence:

```
[1/5] Load 3D Stack & TIFF Header Calibration
       │
[2/5] Nuclear ROI Extraction & Audit Save (*_ROI.tif)
       │
[3/5] Mask Intersect (img_mask * roi_mask) & 3D Connected Components
       │
[4/5] True Signal Extraction (Volumes, K_part, Radial EDT, Null Model)
       │
[5/5] Save CSV Table (_Output_Batch_3d.csv) & Metadata (_metadata.json)
```

---

## Production Code: Batch Iteration & Metadata Serialization

Below is the production logic from [`Batch.py`](file:///C:/Users/franc/Desktop/ExQt_Rezim_A_Final/Batch.py) demonstrating automated discovery, progress tracking, and structured metadata persistence:

```python
# --- EXCERPT FROM Batch.py: run_batch_processing() ---

def run_batch_processing(
    image_dir: Path,
    pixel_size_nm: float,
    z_step_nm: float,
    expansion_factor: float,
    mode_a_enabled: bool = True,
    progress_callback: Callable | None = None,
) -> pd.DataFrame:
    all_rows = []
    metadata_log = {}

    # 1. Discover and pair valid raw/mask files
    raw_files = sorted(image_dir.glob("*.tif"))
    raw_files = [f for f in raw_files if not f.name.endswith(("_Mask.tif", "_ROI.tif"))]

    total_files = len(raw_files)
    for idx, tif_path in enumerate(raw_files):
        mask_path = tif_path.with_name(f"{tif_path.stem}_Mask.tif")
        if not mask_path.exists():
            print(f"Skipping {tif_path.name}: Mask not found ({mask_path.name})")
            continue

        if progress_callback:
            progress_callback(int((idx / total_files) * 100), f"Processing {tif_path.name}")

        # 2. Execute full feature extraction on current pair
        df_pair, meta_pair = process_batch_pair(
            tif_path=tif_path,
            mask_path=mask_path,
            pixel_size_nm=pixel_size_nm,
            z_step_nm=z_step_nm,
            expansion_factor=expansion_factor,
            mode_a_enabled=mode_a_enabled,
        )

        all_rows.append(df_pair)
        metadata_log[tif_path.name] = meta_pair

    # 3. Concatenate into master dataset and persist run metadata
    master_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    
    # Save structured audit JSON
    meta_path = image_dir / "Batch_Run_metadata.json"
    meta_path.write_text(json.dumps(metadata_log, indent=2), encoding="utf-8")
    
    return master_df
```

### Key Technical Features:
- **Resilient Matching**: File discovery explicitly ignores `_Mask.tif` and `_ROI.tif` to avoid processing masks as raw images.
- **Auditable Metadata Log**: The runtime metadata (`Batch_Run_metadata.json`) saves the exact timestamp, git hash / script version, optical calibration, and hardware environment under which the data was acquired.