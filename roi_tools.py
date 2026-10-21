"""Helpers for the manual ROI step: instructions, shape numbering, preview checks, saved ROI reuse."""
from pathlib import Path

import numpy as np
import tifffile

#An ROI covering more than this fraction of slices without any mask object is reported.
EMPTY_SLICE_WARNING_FRACTION = 0.2


def roi_instructions(is_3d, n_shapes, has_saved_roi=False):
    """HTML help shown in the left panel while the user draws nuclei."""
    steps = [
        "Draw <b>one shape around each nucleus</b>: ellipse <b>E</b>, polygon <b>P</b>, "
        "select/move <b>S</b>, delete selected <b>Backspace</b>.",
        "Shapes of different nuclei must not touch.",
    ]
    if is_3d:
        steps.append(
            "Z: a nucleus drawn on <b>one slice</b> is used through the whole stack. Draw it on its "
            "<b>first and last visible slice</b> (overlapping shapes = same nucleus) to leave out "
            "slices without the nucleus."
        )
    steps.append(
        "The number on each shape is the nucleus it will be in the results (<b>!</b> = shapes touch). "
        "<b>Preview ROI</b> shows the final nuclei in all slices; then <b>Confirm ROI</b>."
    )
    if has_saved_roi:
        steps.append("An ROI saved by a previous run exists: <b>Use saved ROI</b> skips drawing.")
    items = "".join(f"<li>{step}</li>" for step in steps)
    return f"<b>Shapes drawn: {n_shapes}</b><ol style='margin-left:-20px'>{items}</ol>"


def shape_nucleus_ids(drawn_labels, cell_labels, n_shapes):
    """cell_id each drawn shape ends up in (shape i is label i+1 of napari's to_labels), or None."""
    drawn = np.asarray(drawn_labels)
    cells = np.asarray(cell_labels)
    ids = []
    for shape_label in range(1, n_shapes + 1):
        values = cells[(drawn == shape_label) & (cells > 0)]
        ids.append(int(np.bincount(values).argmax()) if values.size else None)
    return ids


def roi_coverage_warnings(cell_labels, condensate_mask, max_empty_fraction=EMPTY_SLICE_WARNING_FRACTION):
    """Warn per nucleus when its ROI spans many Z slices that contain no mask object inside it."""
    cells = np.asarray(cell_labels)
    if cells.ndim != 3:
        return []
    mask = np.asarray(condensate_mask) > 0
    warnings = []
    for cell_id in [int(v) for v in np.unique(cells) if v > 0]:
        nucleus = cells == cell_id
        covered = nucleus.any(axis=(1, 2))
        with_objects = (nucleus & mask).any(axis=(1, 2))
        empty = int(np.sum(covered & ~with_objects))
        total = int(np.sum(covered))
        if total and empty / total > max_empty_fraction:
            warnings.append(
                f"Nucleus {cell_id}: {empty} of {total} ROI slices contain no mask objects. "
                "The stack may not cover this nucleus in Z: objects at the first/last slice are then "
                "excluded by QC (see Settings: Keep objects cut by the first/last Z slice)."
            )
    return warnings


def load_saved_roi(path, shape):
    """A drawn ROI saved by an earlier run, or None if missing, of another shape, or whole-field."""
    if not path or not Path(path).is_file():
        return None
    try:
        roi = tifffile.imread(str(path))
    except Exception:
        return None
    if roi.shape != tuple(shape) or roi.max() == 0 or np.all(roi > 0):
        return None
    return roi.astype(np.int32)


def live_shape_numbers(shapes_data, shape_types, image_shape):
    """Nucleus number of each drawn shape, computed from the drawn slices only (fast, no full stack).

    Uses the same grouping and ordering as Batch._interpolate_or_extrude_roi, so the numbers equal
    the final cell_ids. Returns (numbers, error); shapes that cannot be assigned are marked "!".
    """
    import napari
    from Batch import _group_roi_shapes

    n_shapes = len(shapes_data)
    if n_shapes == 0:
        return [], None
    plane_shape = tuple(image_shape[-2:])
    if len(image_shape) == 2:
        return [str(i + 1) for i in range(n_shapes)], None

    z_of_shape = [int(round(float(np.asarray(vertices)[0, 0]))) for vertices in shapes_data]
    drawn_z = sorted(set(z_of_shape))
    #Work in the bounding box of all shapes: row-major pixel order (used for numbering) is unchanged.
    all_yx = np.concatenate([np.asarray(vertices)[:, -2:] for vertices in shapes_data])
    low = np.maximum(np.floor(all_yx.min(axis=0)).astype(int) - 2, 0)
    high = np.minimum(np.ceil(all_yx.max(axis=0)).astype(int) + 3, plane_shape)
    crop_shape = tuple(int(v) for v in np.maximum(high - low, 1))
    stack = np.zeros((len(drawn_z),) + crop_shape, dtype=np.int32)
    for index, (vertices, shape_type) in enumerate(zip(shapes_data, shape_types)):
        #Rasterize each shape only in its own box, then paste it into the crop.
        yx = np.asarray(vertices)[:, -2:] - low
        box_low = np.maximum(np.floor(yx.min(axis=0)).astype(int) - 1, 0)
        box_high = np.minimum(np.ceil(yx.max(axis=0)).astype(int) + 2, crop_shape)
        box_shape = tuple(int(v) for v in np.maximum(box_high - box_low, 1))
        plane = napari.layers.Shapes([yx - box_low], shape_type=shape_type)
        mask = plane.to_masks(mask_shape=box_shape)[0]
        target = stack[drawn_z.index(z_of_shape[index]),
                       box_low[0]:box_low[0] + box_shape[0], box_low[1]:box_low[1] + box_shape[1]]
        target[mask[:target.shape[0], :target.shape[1]]] = index + 1
    try:
        groups = _group_roi_shapes(stack)
    except ValueError as error:
        return ["!"] * n_shapes, str(error)

    numbers = ["?"] * n_shapes
    for cell_id, group in enumerate(groups, start=1):
        for z_index, mask in group:
            for value in np.unique(stack[z_index][mask]):
                if value > 0:
                    numbers[int(value) - 1] = str(cell_id)
    return numbers, None
