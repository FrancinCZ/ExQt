"""Describe saved size calibration without changing measurements or GUI state."""
import json
import math
from pathlib import Path


def load_scale_metadata(csv_path):
    path = Path(csv_path)
    sidecar = path.with_name(path.stem + "_metadata.json")
    try:
        with sidecar.open(encoding="utf-8") as stream:
            value = json.load(stream)
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def scale_correction_label(metadata, columns):
    if "volume_bio_um3" in columns:
        dimension, quantity = 3, "Volume"
    elif "area_bio_um2" in columns:
        dimension, quantity = 2, "Area"
    else:
        return "Size scale unavailable: no explicit biological volume/area column."
    parameters = metadata.get("parameters", {}) if isinstance(metadata, dict) else {}
    value = parameters.get("expansion_factor") if isinstance(parameters, dict) else None
    try:
        factor = float(value)
        divisor = factor ** dimension
        valid = not isinstance(value, bool) and factor > 0 and math.isfinite(factor) and math.isfinite(divisor) and divisor > 0
    except (TypeError, ValueError, OverflowError):
        valid = False
    if not valid:
        return f"{quantity}: biological-scale export; ExF unavailable in saved metadata. Correction not verified."
    return (
        f"{quantity}: estimated pre-expansion scale. Saved linear ExF = {factor:g}x. "
        f"Acquired {quantity.lower()} / {factor:g}^{dimension} (= / {divisor:g})."
    )
