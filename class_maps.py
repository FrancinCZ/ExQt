"""Images that show which object got which class.

Batch writes <image>_Objects.tif (voxel value = object_id). After QC this module colours it by
the reported class: <image>_Classes.tif (one code per voxel, opens over the raw image in
Fiji/napari) and <image>_Classes.png (XY and XZ projection with object ids).
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from Batch import object_map_path

#Voxel codes of Classes.tif.
CLASS_CODES = {
    "background": 0,
    "Globular (Low FA)": 1,
    "Elongated (High FA)": 2,
    "Primary, no shape class": 3,
    "Not primary (QC / size range)": 4,
}
CLASS_COLOURS = ["#000000", "#2ca02c", "#ff7f0e", "#1f77b4", "#7f7f7f"]
#In a projection a primary object is drawn over an excluded one.
_PROJECTION_RANK = np.array([0, 3, 4, 2, 1], dtype=np.uint8)
_RANK_TO_CODE = np.argsort(_PROJECTION_RANK).astype(np.uint8)
MAX_LABELLED_OBJECTS = 300


def classes_paths(tif_name, output_dir):
    stem = Path(tif_name).stem
    return Path(output_dir) / f"{stem}_Classes.tif", Path(output_dir) / f"{stem}_Classes.png"


def _truthy(value):
    if isinstance(value, str):
        return value.strip().lower() in ("true", "1", "yes")
    return bool(value) if pd.notna(value) else False


def object_class_codes(report_df):
    """Class code per object_id of one image (report_df rows of that file)."""
    codes = {}
    for _, row in report_df.iterrows():
        if not _truthy(row.get("primary_qc_valid", False)):
            code = CLASS_CODES["Not primary (QC / size range)"]
        else:
            name = str(row.get("condensate_class", ""))
            if name.startswith("Globular"):
                code = CLASS_CODES["Globular (Low FA)"]
            elif name.startswith("Elongated"):
                code = CLASS_CODES["Elongated (High FA)"]
            else:
                code = CLASS_CODES["Primary, no shape class"]
        codes[int(row["object_id"])] = code
    return codes


def _project(code_map, axis):
    return _RANK_TO_CODE[_PROJECTION_RANK[code_map].max(axis=axis)]


def _draw_panel(ax, image, title, labels, xlabel, ylabel):
    ax.imshow(image, cmap=ListedColormap(CLASS_COLOURS), vmin=0, vmax=len(CLASS_COLOURS) - 1,
              interpolation="nearest", aspect="auto")
    for x, y, text in labels:
        ax.text(x, y, text, color="white", fontsize=5, ha="center", va="center")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def export_class_map(object_map, file_report, tif_name, output_dir, z_step_nm=None, pixel_size_nm=None):
    """Write Classes.tif and Classes.png for one image; returns the two paths."""
    object_map = np.asarray(object_map)
    codes = object_class_codes(file_report)
    lut = np.zeros(int(object_map.max()) + 1, dtype=np.uint8)
    for object_id, code in codes.items():
        if object_id < lut.size:
            lut[object_id] = code
    code_map = lut[object_map]
    tif_path, png_path = classes_paths(tif_name, output_dir)
    tifffile.imwrite(str(tif_path), code_map, compression="zlib", photometric="minisblack")

    primary = file_report[file_report["primary_qc_valid"].map(_truthy)] if "primary_qc_valid" in file_report else file_report.iloc[0:0]
    show_ids = len(primary) <= MAX_LABELLED_OBJECTS
    is_3d = code_map.ndim == 3
    if is_3d:
        #XZ height in proportion to the real Z step, so a Z elongation looks like one.
        z_scale = (z_step_nm / pixel_size_nm) if z_step_nm and pixel_size_nm else 1.0
        xy_height = code_map.shape[1]
        xz_height = code_map.shape[0] * z_scale
        fig, (ax_xy, ax_xz) = plt.subplots(
            2, 1, figsize=(8, 8 * (xy_height + xz_height) / max(code_map.shape[2], 1) + 1.2),
            gridspec_kw={"height_ratios": [xy_height, max(xz_height, 1)]},
        )
        xy_labels = [(r.X_px, r.Y_px, int(r.object_id)) for r in primary.itertuples()] if show_ids else []
        xz_labels = [(r.X_px, r.Z_px, int(r.object_id)) for r in primary.itertuples()] if show_ids else []
        _draw_panel(ax_xy, _project(code_map, 0), "Top view (max projection over Z)", xy_labels, "X (px)", "Y (px)")
        _draw_panel(ax_xz, _project(code_map, 1), "Side view (max projection over Y)", xz_labels, "X (px)", "Z (slice)")
    else:
        fig, ax_xy = plt.subplots(figsize=(8, 8 * code_map.shape[0] / max(code_map.shape[1], 1) + 1.2))
        labels = [(r.X_px, r.Y_px, int(r.object_id)) for r in primary.itertuples()] if show_ids else []
        _draw_panel(ax_xy, code_map, "Objects by class", labels, "X (px)", "Y (px)")

    counts = {code: sum(1 for c in codes.values() if c == code) for code in range(1, len(CLASS_COLOURS))}
    handles = [
        Patch(color=CLASS_COLOURS[code], label=f"{name} ({counts[code]})")
        for name, code in CLASS_CODES.items() if code > 0
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8, frameon=False)
    note = "Numbers = object_id in the tables." if show_ids else f"Object ids not drawn (> {MAX_LABELLED_OBJECTS} primary objects)."
    fig.suptitle(f"{tif_name}: objects by class\n{note} A Z elongation looks round in the top view.", fontsize=10)
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    fig.savefig(png_path, dpi=200)
    plt.close(fig)
    return tif_path, png_path


def export_class_maps(report_df, output_dir, calibrations=None):
    """Classes.tif/png for every image of a run whose Objects.tif exists; returns written paths."""
    written = []
    calibrations = calibrations or {}
    for tif_name, file_report in report_df.groupby("filename"):
        map_path = object_map_path(tif_name, output_dir)
        if not map_path.is_file():
            continue
        calibration = calibrations.get(tif_name, {})
        written.extend(export_class_map(
            tifffile.imread(str(map_path)), file_report, tif_name, output_dir,
            z_step_nm=calibration.get("z_step_nm"), pixel_size_nm=calibration.get("pixel_size_nm"),
        ))
    return written
