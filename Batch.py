import math
from pathlib import Path
import h5py
import napari
import numpy as np
import pandas as pd
import tifffile
import json
from scipy.ndimage import laplace
from skimage.measure import regionprops, label
from pathlib import Path


def _select_focus_slice(volume):
    #The sharpest slice is often the one with the highest variance in the Laplacian of the image
    scores = [laplace(volume[z].astype(np.float32)).var() for z in range(volume.shape[0])]
    return int(np.argmax(scores))


def load_ilastik_h5(h5_path, prob_channel=0):

    with h5py.File(h5_path, "r") as f:
        #Find the dataset
        if "exported_data" in f:
            dataset = f["exported_data"]
        else:
            key = list(f.keys())[0]
            dataset = f[key]
            
        data = dataset[:]
        
        #Axis tags analysis to determine the correct channel and time axes
        if "axistags" in dataset.attrs:
            try:
                axistags_str = dataset.attrs["axistags"]
                if isinstance(axistags_str, bytes):
                    axistags_str = axistags_str.decode('utf-8')
                
                axistags = json.loads(axistags_str)
                axes = [ax["key"] for ax in axistags.get("axes", [])]
                
                if 'c' in axes:
                    c_idx = axes.index('c')
                    data = np.take(data, prob_channel, axis=c_idx)
                    axes.pop(c_idx) 
                
                if 't' in axes:
                    t_idx = axes.index('t')
                    data = np.take(data, 0, axis=t_idx)
                
                return data.astype(np.float32)
                
            except Exception as e:
                print(f"      [Varování] Nelze analyzovat 'axistags', fallback na shape: {e}")

        #FALLBACK: if axistags are not present or cannot be parsed, try to infer the channel axis based on shape heuristics
        data = np.squeeze(data)
        if data.ndim == 4:
            #(Z, Y, X, C)
            if data.shape[-1] <= 10:       
                data = data[..., prob_channel]
            #(Z, C, Y, X)
            elif data.shape[1] <= 10:      
                data = data[:, prob_channel, :, :]
            #(C, Z, Y, X)
            elif data.shape[0] <= 10:      
                data = data[prob_channel, ...]
            else:
                data = data[..., prob_channel]
        elif data.ndim == 3:
            #(Y, X, C)
            if data.shape[-1] <= 10:       
                data = data[..., prob_channel]
            #(C, Y, X)
            elif data.shape[0] <= 10:      
                data = data[prob_channel, ...]

    return data.astype(np.float32)


def process_condensates_hybrid_h5(
    tif_path,
    h5_path,
    mode="3d",  # "3d" | "2d" | "single_slice"
    target_z_slice=None,
    expansion_factor=1.0,
    log_threshold=0.25,  #Threshold for Ilastik probability map to consider a pixel as condensate
    prob_channel=0,     #Which Ilastik channel corresponds to condensates (0 or 1)
    show_napari=True,
    pixel_size_nm=58.0,
    z_step_nm=250.0,
    signal_channel=0,
):
    tif_path = Path(tif_path)
    h5_path = Path(h5_path)
    print(f"\n[1/4] Loading Pair: {tif_path.name} & {h5_path.name}")

    #Loading raw TIF and Ilastik probability map
    img_raw = tifffile.imread(tif_path)
    img_prob = load_ilastik_h5(h5_path, prob_channel=prob_channel)

    #Control of channels
    if img_raw.ndim == 4:
        img_raw = np.take(img_raw, signal_channel, axis=1)

    #Control of dimensions
    if img_prob.shape != img_raw.shape:
        print(f"      WARNING: Inconsistent data shapes! TIF: {img_raw.shape}, H5: {img_prob.shape}")
        if img_prob.shape == img_raw.shape[::-1]: # If axes X and Z are transposed
            img_prob = np.transpose(img_prob)
            print("      -> H5 map was automatically transposed.")

    is_stack = img_raw.ndim == 3

    print(f"      Mode: {mode}")
    #MIP or single slice selection based on mode
    if mode == "single_slice":
        if not is_stack:
            raise ValueError("single_slice mode needs a 3D (Z,Y,X) stack.")
        z_idx = (
            target_z_slice
            if target_z_slice is not None
            else _select_focus_slice(img_raw)
        )
        print(f"      Using Z-slice {z_idx} (auto-focus)")

        img_intensity = img_raw[z_idx]
        img_prob_process = img_prob[z_idx]
        is_3d = False

    elif mode == "2d":
        img_intensity = img_raw if not is_stack else img_raw.max(axis=0)
        img_prob_process = img_prob if not is_stack else img_prob.max(axis=0)
        is_3d = False

    elif mode == "3d":
        if not is_stack:
            raise ValueError("3d mode needs a 3D (Z,Y,X) stack.")
        img_intensity = img_raw
        img_prob_process = img_prob
        is_3d = True

    #Normalization Probability map
    if img_prob_process.max() > 1.5:
        img_prob_process = img_prob_process / img_prob_process.max()

    print("[2/4] Segmenting condensates directly from Ilastik Probability Map...")
    binary_mask = img_prob_process > log_threshold
    labeled_mask = label(binary_mask)

    print("[3/4] Extracting true signal metrics...")
    eff_pixel_size_nm = pixel_size_nm / expansion_factor
    eff_z_step_nm = z_step_nm / expansion_factor

    #Extracting region properties and calculating metrics
    props = regionprops(labeled_mask, intensity_image=img_intensity)
    objects_data = []

    for region in props:
        mean_int = region.intensity_mean
        if region.area < 3:  
            continue

        row = {
            "filename": tif_path.name,
            "mode": mode,
            "is_3d": is_3d,
            "object_id": region.label,
            "mean_intensity": round(mean_int, 2),
            "max_intensity": round(region.intensity_max, 2),
            "integrated_density": round(region.area * mean_int, 2),
        }

        if is_3d:
            voxel_volume_bio_um3 = ((eff_pixel_size_nm**2) * eff_z_step_nm) / 1e9
            row["volume_px"] = region.area
            row["volume_bio_um3"] = round(region.area * voxel_volume_bio_um3, 5)
            row["shape_metric_bio"] = row["volume_bio_um3"]
        else:
            pixel_area_bio_um2 = (eff_pixel_size_nm**2) / 1e6
            row["area_px"] = region.area
            row["area_bio_um2"] = round(region.area * pixel_area_bio_um2, 5)
            row["shape_metric_bio"] = row["area_bio_um2"]

        objects_data.append(row)

    print(f"      Found {len(objects_data)} valid condensates.")

    print("[4/4] Generating Preview...")
    if show_napari:
        viewer = napari.Viewer(title=f"Hybrid H5 Check: {tif_path.name} [{mode}]")

        #Raw Signal
        viewer.add_image(
            img_intensity, name="Raw Signal", colormap="gray", blending="additive"
        )
        
        #Ilastik Probability Map
        prob_max = img_prob_process.max()
        viewer.add_image(
            img_prob_process,
            name="Ilastik Probability",
            colormap="magenta",
            opacity=0.6,
            blending="additive",
            contrast_limits=(0, prob_max if prob_max > 0 else 1)
        )

        #Detected Condensates
        if objects_data:
            coords = [region.centroid for region in props]
            sizes = [max(math.sqrt(region.area) * 2.0, 5.0) for region in props]

            if coords:
                viewer.add_points(
                    coords,
                    size=sizes,
                    name="Detected Condensates",
                    symbol="disc",            
                    face_color=[0, 0, 0, 0],  #Native transparent face color
                    border_color="yellow" if is_3d else "cyan",
                    out_of_slice_display=False #Prevent points from being shown in slices where they dont exist
                )

        print("\n[Napari] Viewer opened. Press Ctrl+C in terminal to abort completely.")
        napari.run()

    return pd.DataFrame(objects_data)


if __name__ == "__main__":
    folder_path = Path(r"C:\Users\franc\Desktop\BRD4")
    all_dataframes = []

    MODE = "3d"  # "3d" | "2d" | "single_slice"
    enable_night_preview = True
    expansion_factor = 1.0

    #Ilastik channel index for condensate probability map (0 or 1)
    ILASTIK_PROB_CHANNEL = 0

    print(f"Looking for TIF + H5 pairs in folder: {folder_path}")

    #Going through all TIF files in the folder and finding corresponding H5 files
    raw_files = [
        f
        for f in sorted(folder_path.glob("*.tif"))
        if "Probabilities" not in f.name and "Final" not in f.name
    ]

    for raw_tif in raw_files:
        #Auto detect corresponding H5 file based on naming convention
        h5_file = raw_tif.with_name(f"{raw_tif.stem}_Probabilities.h5")

        if not h5_file.exists():
            print(
                f"  [SKIPPED] No corresponding H5 file found for {raw_tif.name} (looked for {h5_file.name})"
            )
            continue

        df_file = process_condensates_hybrid_h5(
            tif_path=raw_tif,
            h5_path=h5_file,
            mode=MODE,
            expansion_factor=expansion_factor,
            log_threshold=0.5,  #Threshold for Ilastik probability map 
            prob_channel=ILASTIK_PROB_CHANNEL,
            show_napari=enable_night_preview,
            pixel_size_nm=58.0,
            z_step_nm=250.0,
            signal_channel=0,
        )
        if not df_file.empty:
            all_dataframes.append(df_file)

    if all_dataframes:
        final_df = pd.concat(all_dataframes, ignore_index=True)
        
        folder_name = Path(folder_path).name 
        
        output_csv_filename = f"{folder_name}_Output_Batch_3d.csv"
        final_df.to_csv(output_csv_filename, index=False)
        print("\n[Success] Batch processing with H5 completed.")
        print(f"Results saved to: {output_csv_filename}")
    else:
        print(
            "\n[Error or No Data] No files were processed (check h5 file names)."
        )
