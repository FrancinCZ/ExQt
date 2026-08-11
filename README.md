# ExQt: Expansion microscopy Quantification analysis of Nuclear Condensates

ExQt is an open-source, PySide6-based desktop application designed for the flexible analysis of 3D/2D biological image data, specifically optimized for Expansion Microscopy (ExM). It bridges the gap between raw imaging data, segmentation masks (e.g., from Labkit or Fiji), and interactive visualization using `napari`.

## Features
*   **Interactive 3D/2D Visualization:** Integration with `napari` for visual inspection of raw data, segmentation masks, and analyzed condensates.
*   **Segmentation & Filtering:** Employs user-provided TIF masks combined with minimum size constraints to strictly filter out sub-resolution noise and artifacts.
*   **Batch Processing:** Automated pipeline for processing multiple `.tif` files, outputting comprehensive spatial coordinates (Z, Y, X), biological volumes, equivalent diameters, and signal intensities.
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

environment.yml is provided so that it automatically handles the installation of Python and all required dependencies (including napari and PySide6). Create and activate the environment by running:

```Bash
conda env create -f environment.yml
conda activate exqt-env

```

(Note: If you prefer using pure Python without Conda, you can alternatively install the dependencies via pip install -r requirements.txt)

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

## Quick Start

To help you get familiar with ExQt, a sample dataset is provided in the `example_data/` directory.

1. **Launch the app:** Run `python App.py`.
2. **Choose the folder:** Click **Browse** and select your data directory. ExQt expects each raw image `sample.tif` to have a matching segmentation mask named `sample_Mask.tif` in the same folder.
3. **Set Parameters:** Select a **Process mode** (`3d`, `2d`, or `single_slice`) and adjust settings (e.g., Expansion Factor, Z-step, Pixel size, Min. size) according to your sample metadata.
4. **Run Analysis:** Click **Start analysis** to begin processing the data.
5. **Review:** If you enable the **"Pause and review segmentations"** option, ExQt will pause after each image to let you visually inspect the intermediate steps in `napari` before continuing. Otherwise, processing proceeds automatically.
6. **Output:** The software outputs a CSV file containing spatial coordinates (Z, Y, X), biological volumes, and morphological data of detected condensates. A JSON metadata file is created alongside it for reproducibility. If enabled, Excel statistics and distribution plots are generated automatically in the output folder.

## Configuration Parameters

### Main Settings

The main interface of ExQt provides parameters required to correctly process and scale your data.

* **Input / Output Directories:** Define the source folder containing your raw `.tif` images and corresponding `.tif` masks (named `sample_Mask.tif`). Select the destination folder for results, or check **Save to the same folder as input**.
* **Process mode:** ExQt handles full 3D volumes (`3d`), 2D images via Maximum Intensity Projection (`2d`), or single-slice extraction (`single_slice`). In `single_slice` mode, ExQt automatically selects the sharpest Z-slice using a focus-quality metric (variance of the Laplacian).
* **Expansion Factor:** The physical expansion ratio of your biological sample (e.g., 4.0). Essential for calculating true biological volumes and spatial coordinates. Set to 1.0 for non-expanded data.
* **Min. size (voxels / pixels):** The minimum number of connected voxels (in `3d` mode) or pixels (in `2d`/`single_slice` mode) required to classify a segmented object as a valid condensate, filtering out sub-resolution noise.
* **Process Controls:**
* *Show Napari preview:* Toggles visual rendering of image data and detected condensates in `napari`.
* *Auto-ROI:* When enabled, treats the entire field of view as a single region of interest. If disabled, prompts the user to draw a custom ROI for per-cell (`cell_id`) statistics.
* *Pause and review segmentations:* Halts the batch processing pipeline after each image for visual verification.
* *Generate Excel stats & Generate plots:* Toggles automatic creation of descriptive statistics and distribution plots.



### Advanced Settings

Accessible via the top menu bar (`Settings -> Advanced...`), this menu defines hardware and channel configurations as well as physical dimensions.

* **Pixel Size XY (nm):** The physical lateral dimension of a single pixel as acquired by the microscope (e.g., 205.2 nm). Essential for calculating true areas and volumes.
* **Z-step (nm):** The physical axial distance between consecutive slices in a 3D stack (e.g., 1000.0 nm). Used in 3D mode for accurate biological volume calculations.
* **Signal Channel:** The channel index in your raw multi-channel `.tif` image containing the specific condensate signal to quantify (default is channel 1).
* **DAPI Channel:** The channel index corresponding to the nuclear stain (default is channel 0).

## Acknowledgments and Third-Party Licenses

ExQt is open-source software released under the [MIT License](https://www.google.com/search?q=LICENSE). It relies on several awesome open-source libraries:

* **PySide6** (Qt for Python) for the graphical user interface (LGPLv3).
* **napari** for multidimensional image visualization (BSD 3-Clause License).
* **SciPy**, **NumPy**, and **pandas** for data processing and spatial evaluation (BSD 3-Clause License).

The workflow is tailored for segmentation masks generated via tools like [Labkit](https://imagej.net/plugins/labkit/) or [Ilastik](https://www.ilastik.org/).

```

