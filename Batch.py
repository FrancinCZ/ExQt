from pathlib import Path
import h5py
import napari
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage.io import imread
from skimage.measure import regionprops
from skimage.segmentation import watershed


def process_pair(
    tif_path,
    h5_path,
    expansion_factor= 4.0,
    prob_thresh=0.6,
    min_size_px=20,
    min_area_filter=100,
    show_napari=False,
    pixel_size_nm=58.0,
    z_step_nm=250.0,
    signal_channel=0, 
):
    tif_path = Path(tif_path)
    h5_path = Path(h5_path)

    #Load TIF image
    print("\n[1/6] Loading TIF image...")
    img = tifffile.imread(tif_path)

    #Load Ilastik H5 data
    print("[2/6] Loading Ilastik H5 data...")
    with h5py.File(h5_path, "r") as f:
        raw_data = f["exported_data"]
        print(f"      H5 Shape: {raw_data.shape}, ndim: {raw_data.ndim}")

        if raw_data.ndim == 4:
            chan_idx = signal_channel if raw_data.shape[1] > signal_channel else 0
            prob_map = raw_data[:, chan_idx, :, :]
            is_3d = True

        elif raw_data.ndim == 3:
            if raw_data.shape[0] == 2:
                prob_map = raw_data[signal_channel, :, :]
            else:
                prob_map = raw_data[:, :, signal_channel]
            is_3d = False
        else:
            raise ValueError(f"Unsupported data dimensions: {raw_data.ndim}D")

    #Normalize probability map
    if prob_map.max() > 1.0:
        prob_map_norm = prob_map / 255.0
    else:
        prob_map_norm = prob_map

    #Clean binary mask
    print(f"[3/6] Creating mask and cleaning (min_size={min_size_px})...")
    binary_mask = prob_map_norm > prob_thresh

    labeled_binary, _ = ndi.label(binary_mask)
    component_sizes = np.bincount(labeled_binary.ravel())
    too_small = component_sizes < min_size_px
    too_small_mask = too_small[labeled_binary]

    cleaned_mask = binary_mask.copy()
    cleaned_mask[too_small_mask] = False

    #Detect local maxima for Watershed
    print("[4/6] Starting Watershed segmentation...")
    prob_map_masked = prob_map_norm.copy()
    prob_map_masked[~cleaned_mask] = 0.0

    local_max = ndi.maximum_filter(prob_map_masked, size=3) == prob_map_masked
    local_max[prob_map_masked <= prob_thresh] = False

    markers, num_foci = ndi.label(local_max)
    print(f"      Found {num_foci} condensates.")

    #Perform Watershed segmentation
    print("[5/6] Performing Watershed segmentation...")
    labeled_mask = watershed(-prob_map_norm, markers, mask=cleaned_mask)
    print("      Segmentation completed.")

    if show_napari:
        viewer = napari.Viewer(title=f"Check: {tif_path.name}")
        viewer.add_image(img, name="Former TIF", opacity=0.8)
        viewer.add_image(prob_map, name="Ilastik ProbMap", visible=False)
        viewer.add_labels(
            labeled_mask, name="Segmented Objects", opacity=0.5
        )

        print("\n[Napari] Viewer opened.")
        napari.run()

    #Extract metrics and calculate biological scaling
    print("[6/6] Extracting data")
    exf = expansion_factor if expansion_factor else 1.0

    if is_3d:
        #Voxel dimensions in nanometers divided by expansion factor
        bio_x_nm = pixel_size_nm / exf
        bio_y_nm = pixel_size_nm / exf
        bio_z_nm = z_step_nm / exf

        voxel_volume_bio_nm3 = bio_x_nm * bio_y_nm * bio_z_nm
    else:
        # For 2D images (fallback)
        effective_pixel_size = pixel_size_nm / exf
        pixel_area_bio_nm2 = effective_pixel_size**2

    #Select correct intensity channel from TIF image
    img_intensity = img
    if img.ndim == 4:
        if img.shape[1] in [2, 3, 4]:
            img_intensity = img[:, min(signal_channel, img.shape[1] - 1), :, :]
        elif img.shape[-1] in [2, 3, 4]:
            img_intensity = img[:, :, :, min(signal_channel, img.shape[-1] - 1)]
        elif img.shape[0] in [2, 3, 4]:
            img_intensity = img[min(signal_channel, img.shape[0] - 1), :, :, :]

    img_intensity = np.squeeze(img_intensity)

    props = regionprops(labeled_mask, intensity_image=img_intensity)

    objects_data = []
    for region in props:
        if region.area >= min_area_filter:
            mean_int = region.mean_intensity
            int_density = region.area * mean_int

            row = {
                "filename": tif_path.name,
                "object_id": region.label,
                "mean_intensity": round(mean_int, 2),
                "integrated_density": round(int_density, 2),
                "volume_voxels": region.area,
            }

            if is_3d:
                row["volume_bio_nm3"] = region.area * voxel_volume_bio_nm3
            else:
                row["area_bio_nm2"] = region.area * pixel_area_bio_nm2

            objects_data.append(row)

    return pd.DataFrame(objects_data)



folder_path = Path(r"/home/martinfranc/Stažené/Magnify/2D")
all_dataframes = []

enable_night_preview = False  #Set to True to enable Napari preview for each pair
expansion_factor = 4.0        #Set according to your ExM expansion factor 

for tif_file in folder_path.glob("*.tif"):
    matching_h5 = list(folder_path.glob(f"{tif_file.stem}*.h5"))

    if matching_h5:
        h5_file = matching_h5[0]
        print(f"\n--- Processing pair: {tif_file.name} + {h5_file.name} ---")
        df_pair = process_pair(
            tif_file,
            h5_file,
            expansion_factor=expansion_factor,
            show_napari=enable_night_preview,
            pixel_size_nm=60.0,    #Microscope resolution (0.06 µm)
            z_step_nm=250.0,      #Microscope Z-step (0.25 µm)
            signal_channel=0   
        )
        all_dataframes.append(df_pair)
    else:
        print(f"WARNING: No corresponding H5 file found for {tif_file.name}")

#Save combined results to CSV
if all_dataframes:
    final_df = pd.concat(all_dataframes, ignore_index=True)
    output_csv = folder_path / "Final_Output_Batch.csv"
    final_df.to_csv(output_csv, index=False)

    print("\nBatch completed.")
    print(f"Detected objects: {len(final_df)}")
    print(f"Saved to: {output_csv}")
    print(final_df.head())