---Work in progress---


# ExQt: Quantitative Analysis of Nuclear Condensates via Expansion Microscopy

ExQt is an open-source, PySide6-based desktop application designed for the flexible analysis of 3D/2D biological image data, specifically optimized for Expansion Microscopy (ExM). It bridges the gap between raw imaging data, probability maps (e.g., from Ilastik), and interactive visualization using `napari`.

## Features
*   **Interactive 3D/2D Visualization:** Integration with `napari` for visual inspection of raw data, probability maps, and segmented condensates.
*   **Segmentation:** Employs customizable probability thresholds and minimum voxel constraints to strictly filter out sub-resolution noise and artifacts.
*   **Batch Processing:** Automated pipeline for processing multiple `.tif` files, outputting comprehensive spatial coordinates (Z, Y, X), biological volumes, and signal intensities.
*   **Auditability & Reproducibility:** Automatically generates statistics and graphs with metadata logs (JSON) for every batch, detailing all applied parameters (expansion factor, pixel size, thresholds).

---

## Installation

Installing ExQt is recommended inside an isolated virtual environment using [Conda](https://docs.conda.io/en/latest/miniconda.html) (or Mamba) to prevent dependency conflicts with other Python packages on your system.

### 1. Clone the repository
First, download the source code to your local machine:

```Bash
git clone [https://github.com/FrancinCZ/ExQt.git](https://github.com/FrancinCZ/ExQt.git)
cd ExQt
```
### 2. Create a Conda environment

Create a new virtual environment named exqt-env with Python 3.10 (or your preferred compatible version) and activate it:

```Bash
conda create -n exqt-env python=3.10 -y
conda activate exqt-env
```

### 3. Install dependencies
Install the required packages. Ensure your activated environment is running, then install the dependencies via pip:

```Bash
pip install -r requirements.txt
```
(Note: If you encounter issues with napari or PySide6 on specific operating systems, refer to their official documentation for OS-specific binaries).

### 4. Run the application

Once everything is installed, you can launch the graphical user interface by running:

```Bash
python App.py
```

---


## Acknowledgments and Third-Party Licenses

ExQt is open-source software released under the [MIT License](LICENSE). However, it relies on several excellent open-source projects, which are distributed under their respective licenses:

*   **PySide6** (Qt for Python) is used for the graphical user interface. It is licensed under the [LGPLv3](https://www.gnu.org/licenses/lgpl-3.0.html).
*   **napari** is used for multidimensional image visualization and is licensed under the [BSD 3-Clause License](https://github.com/napari/napari/blob/main/LICENSE).
*   **SciPy**, **NumPy**, and **pandas** are used for data processing and spatial evaluation. They are licensed under the BSD 3-Clause License.

The development of this tool was tailored for workflows involving [Ilastik](https://www.ilastik.org/), an interactive learning and segmentation toolkit (GPLv2).
