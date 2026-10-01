"""Read detector settings from the XML header of a Leica .lif file (image data is not read).

Only the unambiguous case yields an offset: every active detector in photon-counting mode
counts photons without an electronic offset, so the offset is 0 ADU. In analog mode Leica
stores "Offset" as an instrument setting, not in ADU, so no value is derived from it.
"""
from __future__ import annotations

import struct
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

from provenance import file_sha256

LIF_MAGIC = 0x70
LIF_XML_BLOCK = 0x2A
PHOTON_COUNTING = "PhotonCounting"


def read_lif_xml(path):
    """Return the XML header of a .lif file as text; ValueError if it is not a LIF file."""
    path = Path(path)
    with open(path, "rb") as stream:
        head = stream.read(13)
        if len(head) < 13:
            raise ValueError(f"{path.name} is not a Leica LIF file (too short).")
        magic, _size = struct.unpack("<ii", head[:8])
        block, n_chars = struct.unpack("<bi", head[8:13])
        if magic != LIF_MAGIC or block != LIF_XML_BLOCK or n_chars <= 0:
            raise ValueError(f"{path.name} is not a Leica LIF file (unexpected header).")
        raw = stream.read(n_chars * 2)
    if len(raw) != n_chars * 2:
        raise ValueError(f"{path.name}: LIF XML header is truncated.")
    return raw.decode("utf-16-le")


def detector_summary(path):
    """Active detectors (name, type, channel, mode, offset, gain) and how often each setting occurs in the header."""
    path = Path(path)
    try:
        root = ET.fromstring(read_lif_xml(path))
    except ET.ParseError as error:
        raise ValueError(f"{path.name}: LIF XML header cannot be parsed ({error}).") from error

    n_images = sum(1 for element in root.iter("Element") if element.find("Data/Image") is not None)
    counts = Counter()
    for detector in root.iter("Detector"):
        attrs = detector.attrib
        if attrs.get("IsActive") != "1":
            continue
        counts[(
            attrs.get("Name", ""), attrs.get("Type", ""), attrs.get("Channel", ""),
            attrs.get("AcquisitionModeName", ""), attrs.get("Offset", ""), attrs.get("Gain", ""),
        )] += 1
    active = [
        {"name": name, "type": kind, "channel": channel, "mode": mode,
        "offset": offset, "gain": gain, "occurrences": count}
        for (name, kind, channel, mode, offset, gain), count in sorted(counts.items())
    ]
    return {"file": path.name, "n_images": n_images, "active_detectors": active}


def offset_from_lif(path):
    """Return (offset_adu or None, source dict). A value is returned only when every active
    detector is in photon-counting mode; otherwise the message says why not."""
    summary = detector_summary(path)
    source = {
        "method": "lif_photon_counting",
        "file": summary["file"],
        "sha256": file_sha256(path),
        "n_images": summary["n_images"],
        "detectors": summary["active_detectors"],
    }
    detectors = summary["active_detectors"]
    if not detectors:
        source.update(method="lif_unresolved", message="No active detector found in the LIF header.")
        return None, source
    analog = sorted({d["name"] for d in detectors if d["mode"] != PHOTON_COUNTING})
    if analog:
        source.update(
            method="lif_unresolved",
            message=(
                f"Detectors {', '.join(analog)} are not photon counting; their 'Offset' is an instrument "
                "setting, not ADU. Measure a dark frame (laser off, same settings) and enter its median."
            ),
        )
        return None, source
    return 0.0, source
