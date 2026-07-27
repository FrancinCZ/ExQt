import json
from pathlib import Path
import h5py
import napari
import numpy as np
import pandas as pd
import tifffile
from csbdeep.utils import normalize
from skimage.measure import regionprops
from stardist.models import StarDist2D


def process_pair(
    tif_path,
    h5_path,
    model,
    expansion_factor=1.0,
    min_area_filter=100,
    show_napari=True,
    pixel_size_nm=58.0,
    z_step_nm=250.0,
    signal_channel=0,
):
    tif_path = Path(tif_path)
    h5_path = Path(h5_path)

    #Load TIF image
    print("\n[1/5] Loading TIF image...")
    img = tifffile.imread(tif_path)

    #Load Ilastik H5 data
    print("[2/5] Loading Ilastik H5 data...")
    with h5py.File(h5_path, "r") as f:
        raw_data = f["exported_data"]
        print(f"      H5 Shape: {raw_data.shape}, ndim: {raw_data.ndim}")

        #Read axistags if present, to know real axis order
        axis_keys = None
        if "axistags" in raw_data.attrs:
            try:
                tags = json.loads(raw_data.attrs["axistags"])["axes"]
                axis_keys = "".join(t["key"].lower() for t in tags)
            except Exception:
                axis_keys = None

        raw_data = raw_data[...]

        if axis_keys is not None:
            #Use axistags for a reliable channel/axis selection
            if "c" in axis_keys:
                c_idx = axis_keys.index("c")
                prob_map = np.take(
                    raw_data, indices=signal_channel, axis=c_idx
                )
            else:
                prob_map = raw_data
            is_3d = "z" in axis_keys

        elif raw_data.ndim == 4:
            chan_idx = (
                signal_channel if raw_data.shape[1] > signal_channel else 0
            )
            prob_map = raw_data[:, chan_idx, :, :]
            is_3d = True

        elif raw_data.ndim == 3:
            # Fallback (no axistags)
            if raw_data.shape[0] == 2:
                prob_map = raw_data[signal_channel, :, :]
                is_3d = False
            elif raw_data.shape[-1] in (2, 3, 4):
                prob_map = raw_data[:, :, signal_channel]
                is_3d = False
            else:
                prob_map = raw_data
                is_3d = True
        else:
            raise ValueError(f"Unsupported data dimensions: {raw_data.ndim}D")

        prob_map = np.squeeze(prob_map)
        if prob_map.ndim == 3:
            prob_map = np.max(prob_map, axis=0)

    #Normalize probability map
    if prob_map.dtype == np.uint16:
        prob_map_norm = prob_map.astype(np.float32) / 65535.0
    elif prob_map.max() > 1.0:
        prob_map_norm = prob_map / 255.0
    else:
        prob_map_norm = prob_map

    #StarDist segmentation
    print("[3/5] Starting StarDist segmentation...")
    img_norm = normalize(prob_map_norm, 1, 99.8)
    labeled_mask, details = model.predict_instances(img_norm)
    num_foci = len(details["coord"])
    print(f"      Found {num_foci} condensates.")

    #Optional Napari Preview
    if show_napari:
        viewer = napari.Viewer(title=f"Check: {tif_path.name}")
        viewer.add_image(img, name="Former TIF", opacity=0.8)
        viewer.add_image(prob_map, name="Ilastik ProbMap", visible=False)
        viewer.add_labels(labeled_mask, name="Segmented Objects", opacity=0.5)
        print("\n[Napari] Viewer opened.")
        napari.run()

    #Extract metrics and calculate biological scaling
    print("[4/5] Extracting data...")
    exf = expansion_factor if expansion_factor else 1.0

    if is_3d:
        #Voxel dimensions in nanometers divided by expansion factor
        bio_x_nm = pixel_size_nm / exf
        bio_y_nm = pixel_size_nm / exf
        bio_z_nm = z_step_nm / exf
        voxel_volume_bio_nm3 = bio_x_nm * bio_y_nm * bio_z_nm
    else:
        #For 2D images (fallback)
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

    if img_intensity.ndim == 3 and labeled_mask.ndim == 2:
        img_intensity = np.max(img_intensity, axis=0)

    #Make sure TIF and H5-derived arrays actually line up before regionprops
    assert img_intensity.shape == labeled_mask.shape, (
        f"Shape mismatch: TIF {img_intensity.shape} vs mask {labeled_mask.shape}"
    )

    props = regionprops(labeled_mask, intensity_image=img_intensity)
    objects_data = []


    max_area_filter = 800   
    min_intensity_filter = 20  

    for region in props:

        if (min_area_filter <= region.area <= max_area_filter) and (region.mean_intensity > min_intensity_filter):
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



if __name__ == "__main__":
    folder_path = Path(r"C:\Users\franc\Desktop\MED1")
    all_dataframes = []

    enable_night_preview = True  #Napari preview
    expansion_factor = 1.0  #ExM factor

    print("Loading StarDist model...")
    model = StarDist2D.from_pretrained("2D_versatile_fluo")

    for tif_file in sorted(folder_path.glob("*.tif")):
        matching_h5 = sorted(folder_path.glob(f"{tif_file.stem}*.h5"))

        if matching_h5:
            h5_file = matching_h5[0]
            if len(matching_h5) > 1:
                print(
                    f"WARNING: Multiple H5 matches for {tif_file.name}, using {h5_file.name}"
                )

            print(f"\nProcessing pair: {tif_file.name} + {h5_file.name}")

            df_pair = process_pair(
                tif_file,
                h5_file,
                model,
                expansion_factor=expansion_factor,
                show_napari=enable_night_preview,
                pixel_size_nm=58.0,
                z_step_nm=250.0,
                signal_channel=0,
            )
            all_dataframes.append(df_pair)
        else:
            print(
                f"WARNING: No corresponding H5 file found for {tif_file.name}"
            )

    #Save into CSV
    if all_dataframes:
        final_df = pd.concat(all_dataframes, ignore_index=True)
        output_csv = folder_path / "Final_Output_Batch.csv"
        final_df.to_csv(output_csv, index=False)

        print("\nBatch completed.")
        print(f"Detected objects: {len(final_df)}")
        print(f"Saved to: {output_csv}")
        print(final_df.head())