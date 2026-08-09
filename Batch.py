import math
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import tifffile
import json
from scipy.ndimage import laplace
from skimage.measure import regionprops, label
from skimage.filters import gaussian

def _select_focus_slice(volume):
    scores = [laplace(volume[z].astype(np.float32)).var() for z in range(volume.shape[0])]
    return int(np.argmax(scores))


def _resolve_channel_axis(img_raw, expected_axis=1, max_channels=6):
    if img_raw.ndim != 4:
        return expected_axis

    sizes = img_raw.shape
    if sizes[expected_axis] <= max_channels:
        return expected_axis

    candidates = [ax for ax, s in enumerate(sizes) if s <= max_channels]
    if len(candidates) == 1:
        found_axis = candidates[0]
        print(f"      WARNING: axis {expected_axis} has size {sizes[expected_axis]}, which doesn't "
              f"look like a channel axis. Found a plausible channel axis at position {found_axis} "
              f"(size {sizes[found_axis]}) instead - auto-adjusting. Full shape was {sizes}. "
              f"Please verify this is correct for your data; pass channel_axis= explicitly to override.")
        return found_axis

    raise ValueError(
        f"Could not confidently identify the channel axis for TIF shape {sizes}. "
        f"Expected axis {expected_axis} to hold <= {max_channels} channels, but it has "
        f"{sizes[expected_axis]}, and {'no' if not candidates else 'more than one'} other axis "
        f"looks like a plausible channel axis either. Pass channel_axis= explicitly to "
        f"process_condensates_h5() to resolve this manually."
    )


def load_ilastik_h5(h5_path, prob_channel=0):
    with h5py.File(h5_path, "r") as f:
        #Find the dataset
        if "exported_data" in f:
            dataset = f["exported_data"]
        else:
            key = list(f.keys())[0]
            dataset = f[key]
            
        data = dataset[:]
        
        #Axis tags 
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
                print(f"WARNING: Ilastik axistags not found or invalid for {h5_path.name}. Using fallback heuristic based on shape: {data.shape}")

        #FALLBACK
        data = np.squeeze(data)
        if data.ndim == 4:
            if data.shape[-1] <= 10:       
                data = data[..., prob_channel]
            elif data.shape[1] <= 10:      
                data = data[:, prob_channel, :, :]
            elif data.shape[0] <= 10:      
                data = data[prob_channel, ...]
            else:
                data = data[..., prob_channel]
        elif data.ndim == 3:
            if data.shape[-1] <= 10:       
                data = data[..., prob_channel]
            elif data.shape[0] <= 10:      
                data = data[prob_channel, ...]

    return data.astype(np.float32)


def process_condensates_h5(
    tif_path,
    h5_path,
    mode="3d",  #3d  2d single_slice
    target_z_slice=None,
    expansion_factor=1.0,
    prob_threshold=0.3,  #Threshold for Ilastik probability map
    sigma=1.0,
    min_voxels=5,
    prob_channel=0,      
    show_napari=True,
    pixel_size_nm=58.0,
    z_step_nm=250.0,
    signal_channel=1,    #Ch01 = BRD4/MED1 
    dapi_channel=0,      #Ch00 = DAPI 
    channel_axis=1,      #Which axis of a 4D TIF holds the channel dimension
    auto_roi=True,
    send_layer_func=None,
    request_roi_func=None,
):
    tif_path = Path(tif_path)
    h5_path = Path(h5_path)
    print(f"\n[1/5] Loading Pair: {tif_path.name} & {h5_path.name}")

    if send_layer_func:
        send_layer_func({"type": "clear_layers"})

    #TIF and PRobability map 
    img_raw = tifffile.imread(tif_path)
    img_prob = load_ilastik_h5(h5_path, prob_channel=prob_channel)

    if img_raw.ndim == 4:
        ch_axis = _resolve_channel_axis(img_raw, expected_axis=channel_axis)
        img_dapi = np.take(img_raw, dapi_channel, axis=ch_axis)
        img_intensity = np.take(img_raw, signal_channel, axis=ch_axis)
    else:
        img_dapi = img_raw
        img_intensity = img_raw

    if img_prob.shape != img_intensity.shape:
        print(f"WARNING: Inconsistent data shapes! TIF: {img_intensity.shape}, H5: {img_prob.shape}")
        if img_prob.shape == img_intensity.shape[::-1]: 
            img_prob = np.transpose(img_prob)
            print("H5 map was automatically transposed.")

    if img_prob.shape != img_intensity.shape:
        raise ValueError(f"Shape mismatch. TIF intensity shape {img_intensity.shape} does not match H5 probability shape {img_prob.shape} even after transposing.")

    is_stack = img_intensity.ndim == 3

    print(f"      Mode: {mode}")
    #MIP or single slice selection
    if mode == "single_slice":
        if not is_stack:
            raise ValueError("single_slice mode needs a 3D (Z,Y,X) stack.")
        z_idx = (target_z_slice if target_z_slice is not None else _select_focus_slice(img_intensity))
        print(f"      Using Z-slice {z_idx} (auto-focus)")

        img_intensity = img_intensity[z_idx]
        img_dapi_process = img_dapi[z_idx]
        img_prob_process = img_prob[z_idx]
        is_3d = False

    elif mode == "2d":
        img_intensity = img_intensity if not is_stack else img_intensity.max(axis=0)
        img_dapi_process = img_dapi if not is_stack else img_dapi.max(axis=0)
        img_prob_process = img_prob if not is_stack else img_prob.max(axis=0)
        is_3d = False

    elif mode == "3d":
        if not is_stack:
            raise ValueError("3d mode needs a 3D (Z,Y,X) stack.")
        img_intensity = img_intensity
        img_dapi_process = img_dapi
        img_prob_process = img_prob
        is_3d = True

    if img_prob_process.max() > 1.5:
        img_prob_process = img_prob_process / img_prob_process.max()

    # ROI
    print(f"\n[2/5] ROI extraction (Auto-ROI: {auto_roi})")
    if auto_roi or request_roi_func is None:
        roi_mask = np.ones_like(img_prob_process, dtype=bool)
        if is_3d:
            extruded_mask = np.ones_like(img_prob_process, dtype=int)
        else:
            extruded_mask = np.ones_like(img_prob_process, dtype=int)
    else:
        if send_layer_func:
            send_layer_func({
                "type": "image",
                "name": f"Signal ({tif_path.name})",
                "data": img_intensity,
                "kwargs": {"colormap": "gray", "blending": "additive"}
            })
            
            prob_max = img_prob_process.max()
            send_layer_func({
                "type": "image",
                "name": "Ilastik Tip Probability",
                "data": img_prob_process,
                "kwargs": {
                    "colormap": "magenta", 
                    "opacity": 0.6,
                    "blending": "additive",
                    "contrast_limits": (0, prob_max if prob_max > 0 else 1)
                }
            })
        
        print("Waiting for user to draw ROI in Napari...")
        extruded_mask = request_roi_func(img_intensity.shape, is_3d)
        roi_mask = extruded_mask > 0

    #Segmentation
    print("\n[3/5] Segmenting condensates directly from Ilastik Probability Map...")
    
    img_prob_process = img_prob_process * roi_mask

    img_prob_smoothed = gaussian(img_prob_process, sigma=sigma)

    binary_mask = img_prob_smoothed > prob_threshold
    labeled_mask = label(binary_mask)

    print("[4/5] Extracting true signal metrics...")
    eff_pixel_size_nm = pixel_size_nm / expansion_factor
    eff_z_step_nm = z_step_nm / expansion_factor

    props = regionprops(labeled_mask, intensity_image=img_intensity)
    objects_data = []

    for region in props:
        mean_int = region.intensity_mean
        if region.area < min_voxels:
            continue
            
        centroid_coords = tuple(int(round(c)) for c in region.centroid)
        centroid_values = region.centroid
        if len(centroid_values) >= 3:
            z_px, y_px, x_px = [round(float(v), 3) for v in centroid_values[:3]]
        elif len(centroid_values) == 2:
            z_px, y_px, x_px = 0.0, round(float(centroid_values[0]), 3), round(float(centroid_values[1]), 3)
        else:
            z_px = y_px = x_px = np.nan
        
        try:
            cell_id = extruded_mask[centroid_coords] 
        except IndexError:
            continue
            
        if cell_id == 0:
            continue

        row = {
            "filename": tif_path.name,
            "mode": mode,
            "is_3d": is_3d,
            "cell_id": cell_id,
            "object_id": region.label,
            "Z_px": z_px,
            "Y_px": y_px,
            "X_px": x_px,
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

    print(f"      Found {len(objects_data)} valid condensates inside ROI.")

    print("[5/5] Generating Preview")
    if send_layer_func:
        
        send_layer_func({
            "type": "image",
            "name": f"Raw Signal ({tif_path.name})",
            "data": img_intensity,
            "kwargs": {"colormap": "gray", "blending": "additive"}
        })
        
        send_layer_func({
            "type": "labels",
            "name": "ROI Boundaries",
            "data": roi_mask.astype(int),
            "kwargs": {"opacity": 0.2}
        })
        
        prob_max = img_prob_process.max()
        send_layer_func({
            "type": "image",
            "name": "Ilastik Probability",
            "data": img_prob_process,
            "kwargs": {
                "colormap": "magenta", 
                "opacity": 0.6, 
                "blending": "additive",
                "contrast_limits": (0, prob_max if prob_max > 0 else 1)
            }
        })


        coords = [region.centroid for region in props]
    
        if objects_data and len(coords) > 0:
            if is_3d:
                sizes = [max((3 * region.area / (4 * math.pi))**(1/3) * 2.0, 3.0) for region in props]
            else:
                sizes = [max(math.sqrt(region.area / math.pi) * 2.0, 3.0) for region in props]
                
            send_layer_func({
                "type": "points",
                "name": "Detected Condensates",
                "data": coords,
                "kwargs": {
                    "size": sizes,
                    "symbol": "disc",
                    "face_color": "yellow" if is_3d else "cyan",
                    "out_of_slice_display": False
                }
            })

    return pd.DataFrame(objects_data)


if __name__ == "__main__":
    folder_path = Path(r"")
    all_dataframes = []

    MODE = "3d"  #"3d"  "2d"  "single_slice"
    enable_night_preview = True
    expansion_factor = 1.0

    ILASTIK_PROB_CHANNEL = 0
    DAPI_CHANNEL = 0
    SIGNAL_CHANNEL = 1

    print(f"Looking for TIF + H5 pairs in folder: {folder_path}")

    raw_files = [
        f
        for f in sorted(folder_path.glob("*.tif"))
        if "Probabilities" not in f.name and "Final" not in f.name
    ]

    for raw_tif in raw_files:
        h5_file = raw_tif.with_name(f"{raw_tif.stem}_Probabilities.h5")

        if not h5_file.exists():
            print(f"  No corresponding H5 file found for {raw_tif.name} (looked for {h5_file.name})")
            continue

        df_file = process_condensates_h5(
            tif_path=raw_tif,
            h5_path=h5_file,
            mode=MODE,
            expansion_factor=expansion_factor,
            prob_threshold=0.3, 
            sigma=1.0,
            prob_channel=ILASTIK_PROB_CHANNEL,
            show_napari=enable_night_preview,
            pixel_size_nm=58.0,
            z_step_nm=250.0,
            signal_channel=SIGNAL_CHANNEL,
            dapi_channel=DAPI_CHANNEL,
        )
        if not df_file.empty:
            all_dataframes.append(df_file)

    if all_dataframes:
        final_df = pd.concat(all_dataframes, ignore_index=True)
        folder_name = Path(folder_path).name 
        output_csv_filename = f"{folder_name}_Output_Batch_{MODE}.csv"
        final_df.to_csv(output_csv_filename, index=False)
        print("\nSuccess - Batch processing done.")
        print(f"Results saved to: {output_csv_filename}")
    else:
        print("\n Error or No Data - No files were processed.")
