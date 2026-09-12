# Radial Profiling (Mode A: Core–Shell Architecture)

Are condensates uniform droplets, or do they have layers like an onion?

Mode A answers this by dividing each 3D condensate into three concentric layers: **Core**, **Middle**, and **Shell**. This allows you to check whether a protein is concentrated deep inside or sits mostly on the surface.

---

## The Three Layers (Peeling the Onion)

Using the 3D distance from the condensate surface to its center (Euclidean Distance Transform), every voxel inside the object is assigned to one of three layers:

| Layer | Region | Plain-language meaning |
| :--- | :--- | :--- |
| **Core** | Deepest center ($d > 2/3$) | The innermost heart of the condensate. |
| **Middle** | Intermediate zone ($1/3 < d \le 2/3$) | The transition layer between core and surface. |
| **Shell** | Outer surface ($0 < d \le 1/3$) | The outer boundary facing the surrounding nucleoplasm. |

---

## Measuring Internal Gradients

For each layer, ExQt measures the average fluorescence brightness:
- $I_{\text{core}}$: Average brightness in the center
- $I_{\text{middle}}$: Average brightness in the middle layer
- $I_{\text{shell}}$: Average brightness at the surface

The difference between core and shell brightness gives the gradient:

$$\Delta I_{\text{core-shell}} = I_{\text{core}} - I_{\text{shell}}$$

- **$\Delta I > 0$ (Positive):** The center is brighter than the surface — the protein concentrates in the core.
- **$\Delta I \approx 0$ (Zero):** The droplet is roughly uniform throughout.
- **$\Delta I < 0$ (Negative):** The protein is enriched at the rim or shell.

---


## Production Code: Radial Partitioning (`rezim_a_core_shell.py`)

Below is the production implementation from [`rezim_a_core_shell.py`](file:///C:/Users/franc/Desktop/ExQt_Rezim_A_Final/rezim_a_core_shell.py):

```python
#EXCERPT FROM rezim_a_core_shell.py: split_core_middle_shell()

def split_core_middle_shell(
    object_mask: np.ndarray,
    min_core_voxels: int = 20,
    sampling: tuple | None = None,
    shell_end: float = 1.0 / 3.0,
    core_start: float = 2.0 / 3.0,
) -> dict:
    mask = np.asarray(object_mask, dtype=bool)

    #1. Pad by 1 voxel so bounding-box tangent faces are treated as external edges
    padded_mask = np.pad(mask, 1, mode="constant", constant_values=False)
    padded_distance = distance_transform_edt(padded_mask, sampling=sampling)
    slices = tuple(slice(1, -1) for _ in range(mask.ndim))
    distance_from_edge = padded_distance[slices]
    max_distance = float(distance_from_edge.max())

    if max_distance == 0.0:
        raise ValueError("Object is too small to form radial layers.")

    normalized_distance = distance_from_edge / max_distance

    #2. Fixed thirds partition scheme
    shell = mask & (normalized_distance <= shell_end)
    middle = mask & (normalized_distance > shell_end) & (normalized_distance <= core_start)
    core = mask & (normalized_distance > core_start)

    core_voxels = int(core.sum())
    complete_coverage = int(shell.sum() + middle.sum() + core.sum()) == int(mask.sum())

    return {
        "core": core,
        "middle": middle,
        "shell": shell,
        "normalized_distance": normalized_distance,
        "qc": {
            "object_voxels": int(mask.sum()),
            "core_voxels": core_voxels,
            "max_distance": max_distance,
            "core_valid": core_voxels >= min_core_voxels,
            "complete_coverage": complete_coverage,
        },
    }
```