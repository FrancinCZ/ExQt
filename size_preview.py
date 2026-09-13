from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from skimage.measure import regionprops

from Batch import (
    _positive_finite,
    _prepare_labeled_mask,
    _read_tiff_with_axes,
    _resolve_channel_axis,
    _select_focus_slice,
    get_metadata_from_tif,
)



PREVIEW_COLUMNS = [
    "filename",
    "mask_filename",
    "mode",
    "object_id",
    "raw_size",
    "biological_size",
    "unit",
]


def _empty_preview_table():
    """Return an empty result with the same schema as a populated preview."""
    return pd.DataFrame(columns=PREVIEW_COLUMNS)


def _load_processing_mask(source_path, mask_path, mode, signal_channel):
    """Load a source/mask pair and reproduce Batch.process_condensates geometry."""
    image, image_axes = _read_tiff_with_axes(source_path)
    mask, _ = _read_tiff_with_axes(mask_path)

    if image.ndim == 4:
        channel_axis = image_axes.index("C") if image_axes.count("C") == 1 else _resolve_channel_axis(image, expected_axis=1)
        channel_count = image.shape[channel_axis]
        if not 0 <= signal_channel < channel_count:
            raise ValueError(
                f"Signal channel {signal_channel} is outside channel axis "
                f"{channel_axis} with {channel_count} channels for {source_path.name}."
            )
        image = np.take(image, signal_channel, axis=channel_axis)

    if image.ndim not in (2, 3):
        raise ValueError(
            f"Unsupported source shape {image.shape} for {source_path.name}; "
            "the selected signal image must be 2D or 3D."
        )

    if mask.shape != image.shape:
        if mask.shape == image.shape[::-1]:
            mask = np.transpose(mask)
        else:
            raise ValueError(
                f"TIF shape {image.shape} does not match mask shape {mask.shape} "
                f"for {source_path.name}."
            )

    is_stack = image.ndim == 3
    if mode == "single_slice":
        if not is_stack:
            raise ValueError(
                f"single_slice mode needs a 3D (Z,Y,X) stack; "
                f"{source_path.name} has shape {image.shape}."
            )
        focus_slice = _select_focus_slice(image)
        return mask[focus_slice]

    if mode == "2d":
        return mask if not is_stack else mask.max(axis=0)

    if mode == "3d":
        if not is_stack:
            raise ValueError(
                f"3d mode needs a 3D (Z,Y,X) stack; "
                f"{source_path.name} has shape {image.shape}."
            )
        return mask

    raise ValueError(
        f"Unsupported processing mode {mode!r}; expected '3d', '2d', or 'single_slice'."
    )


def collect_size_preview(
    input_folder,
    mode="3d",
    expansion_factor=1.0,
    pixel_size_nm=None,
    z_step_nm=None,
    min_voxels=5,
    signal_channel=1,
):

    folder = Path(input_folder)
    if not folder.is_dir():
        raise ValueError(f"Input folder does not exist or is not a directory: {folder}")

    if mode not in {"3d", "2d", "single_slice"}:
        raise ValueError(
            f"Unsupported processing mode {mode!r}; expected '3d', '2d', or 'single_slice'."
        )

    pixel_size_nm = _positive_finite(pixel_size_nm, "Pixel size XY")
    expansion_factor = _positive_finite(expansion_factor, "Expansion factor")
    if mode == "3d":
        z_step_nm = _positive_finite(z_step_nm, "Z-step")

    try:
        min_voxels = int(min_voxels)
    except (TypeError, ValueError) as error:
        raise ValueError("Raw minimum pixels/voxels must be a positive integer.") from error
    if min_voxels < 1:
        raise ValueError("Raw minimum pixels/voxels must be at least 1.")

    try:
        signal_channel = int(signal_channel)
    except (TypeError, ValueError) as error:
        raise ValueError("Signal channel must be a non-negative integer.") from error
    if signal_channel < 0:
        raise ValueError("Signal channel must be a non-negative integer.")

    source_files = sorted(
        path
        for path in folder.iterdir()
        if path.is_file()
        and path.suffix.lower() in {".tif", ".tiff"}
        and not path.stem.lower().endswith("_mask")
        and "final" not in path.stem.lower()
        and "_alignment_" not in path.stem.lower()
    )
    if not source_files:
        return _empty_preview_table()

    rows = []
    unit = "µm³" if mode == "3d" else "µm²"

    for source_path in source_files:
        # Determine effective calibration for this file, falling back to batch parameters
        file_meta = get_metadata_from_tif(source_path)
        file_xy_nm = pixel_size_nm
        file_z_nm = z_step_nm
        if file_meta:
            if "pixel_size" in file_meta and file_meta["pixel_size"] is not None:
                file_xy_nm = float(file_meta["pixel_size"])
            if "z_step" in file_meta and file_meta["z_step"] is not None:
                file_z_nm = float(file_meta["z_step"])

        file_eff_xy_nm = file_xy_nm / expansion_factor
        if mode == "3d":
            file_eff_z_nm = (file_z_nm if file_z_nm is not None else 350.0) / expansion_factor
            file_size_per_element = (file_eff_xy_nm**2 * file_eff_z_nm) / 1e9
        else:
            file_size_per_element = file_eff_xy_nm**2 / 1e6

        mask_candidates = [
            folder / f"{source_path.stem}_Mask.tif",
            folder / f"{source_path.stem}_Mask.tiff",
        ]
        mask_path = next((path for path in mask_candidates if path.is_file()), mask_candidates[0])
        if not mask_path.is_file():
            raise ValueError(f"Missing matching mask for {source_path.name}: {mask_path.name}")

        processing_mask = _load_processing_mask(
            source_path,
            mask_path,
            mode,
            signal_channel,
        )
        labeled_mask = _prepare_labeled_mask(processing_mask)

        for region in regionprops(labeled_mask):
            if region.area < min_voxels:
                continue
            rows.append(
                {
                    "filename": source_path.name,
                    "mask_filename": mask_path.name,
                    "mode": mode,
                    "object_id": int(region.label),
                    "raw_size": int(region.area),
                    "biological_size": float(region.area * file_size_per_element),
                    "unit": unit,
                }
            )


    if not rows:
        return _empty_preview_table()
    return pd.DataFrame.from_records(rows, columns=PREVIEW_COLUMNS)
