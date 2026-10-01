"""Record which code and which input files produced a run (VALIDATION_PROTOCOL §2)."""
from __future__ import annotations

import hashlib
import platform
import subprocess
import sys
from importlib import metadata as importlib_metadata
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent

# Modules whose content determines the numbers in the output CSV and reports.
SOURCE_FILES = (
    "App.py", "Batch.py", "calibration_policy.py", "defaults.py", "lif_metadata.py", "null_model.py", "partitioning_plots.py",
    "postprocessing.py", "provenance.py", "reference_values.py", "rezim_a_anisotropy.py", "rezim_a_core_shell.py",
    "rezim_a_metrics.py", "scale_reporting.py", "size_preview.py", "stack_aligner.py",
)
PACKAGES = (
    "numpy", "scipy", "scikit-image", "pandas", "tifffile", "matplotlib", "seaborn",
    "openpyxl", "napari", "PySide6",
)


def file_sha256(path, chunk_size=1 << 20):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git(args, cwd):
    try:
        result = subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout.strip() if result.returncode == 0 else None


def code_provenance(project_dir=PROJECT_DIR):
    """Git revision (if available), SHA-256 of every analysis module, Python and package versions.
    The file hashes identify the code even without git (e.g. a copied or zipped folder)."""
    project_dir = Path(project_dir)
    commit = _git(["rev-parse", "HEAD"], project_dir)
    status = _git(["status", "--porcelain", "--untracked-files=no"], project_dir) if commit else None
    packages = {}
    for name in PACKAGES:
        try:
            packages[name] = importlib_metadata.version(name)
        except importlib_metadata.PackageNotFoundError:
            packages[name] = None
    return {
        "git_commit": commit,
        # True = the tracked source differs from git_commit; the file hashes are then the reference.
        "git_dirty": bool(status) if commit else None,
        "source_sha256": {
            name: file_sha256(project_dir / name)
            for name in SOURCE_FILES if (project_dir / name).is_file()
        },
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": packages,
    }


def input_file_record(raw_path, mask_path=None):
    """Name, size and SHA-256 of a source TIFF, its mask and any alignment sidecars."""
    raw_path = Path(raw_path)
    record = {
        "name": raw_path.name,
        "size_bytes": raw_path.stat().st_size,
        "sha256": file_sha256(raw_path),
    }
    if mask_path is not None and Path(mask_path).is_file():
        record["mask_name"] = Path(mask_path).name
        record["mask_sha256"] = file_sha256(mask_path)
    stem = raw_path.name
    for suffix in (".ome.tiff", ".ome.tif", ".tiff", ".tif"):
        if stem.lower().endswith(suffix):
            stem = stem[:-len(suffix)]
            break
    sidecars = {}
    for name in (f"{stem}_drift.csv", f"{stem}_alignment_valid.tif"):
        candidate = raw_path.with_name(name)
        if candidate.is_file():
            sidecars[name] = file_sha256(candidate)
    record["sidecars"] = sidecars
    return record
