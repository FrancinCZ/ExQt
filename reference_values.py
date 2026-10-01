"""Reference values used by classification and plots, each with its meaning and status.

Change a value here only; it is recorded in metadata.json ("reference_values") and the
FA boundary is part of the QC fingerprint because it defines the Globular/Elongated classes.
FA equivalents are for a homogeneous prolate ellipsoid a:1:1 (FA from covariance eigenvalues).
"""

REFERENCE_VALUES = {
    "fa_globular_boundary": {
        "value": 0.65,
        "meaning": "A_object below this value is classified 'Globular (Low FA)', otherwise 'Elongated (High FA)'.",
        "geometry": "FA 0.65 ≈ axis ratio 1.84 of a prolate ellipsoid.",
        "status": "working value [UNVERIFIED] - not derived from data or literature.",
    },
    "fa_psf_elongation_limit": {
        "value": 0.82,
        "meaning": "Reference line: elongation comparable to an optical PSF, so shape alone is not interpretable.",
        "geometry": "FA 0.82 = axis ratio 2.5 of a prolate ellipsoid.",
        "status": (
            "working value [UNVERIFIED] - depends on the set-up; a confocal PSF near the "
            "coverslip has ratio ~2.05 (FA ~0.72) and deep in mounting medium up to ~5 (FA ~0.96) "
            "(PLoS One 2015, doi:10.1371/journal.pone.0121096; values not checked in the primary "
            "source); not measured for the ExM gel. Measure with beads before using it as a limit."
        ),
    },
    "k_low_contrast": {
        "value": 1.5,
        "meaning": "Reference band: K_part below this value is shown as low contrast to the nucleoplasm.",
        "geometry": "",
        "status": "working value [UNVERIFIED] - visual guide only, not a validity criterion.",
    },
}


def reference_value(name):
    """Current numeric value of a reference entry (looked up at call time)."""
    return float(REFERENCE_VALUES[name]["value"])
