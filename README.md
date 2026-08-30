# ExQt: Expansion Microscopy Quantification of Nuclear Condensates

**ExQt** is an open-source desktop application for quantitative 3D analysis of biomolecular condensates and sub-cellular assemblies in **Expansion Microscopy (ExM)** and standard confocal microscopy. It combines raw 3D TIFF intensity stacks with user-provided segmentation masks, applies physical calibration and the expansion factor ($\text{ExF}$), and extracts single-cell thermodynamic partitioning coefficients ($K_{\text{part}}$), 2D biophysical phase diagrams, geometric fractional anisotropy ($FA$), and radial layer profiles.

ExQt does **not** perform segmentation itself. Masks can be prepared in software such as [Labkit](https://imagej.net/plugins/labkit/), [ilastik](https://www.ilastik.org/), Fiji, or any tool capable of exporting binary or instance TIFF masks.

---

## Main Features

- **3D Thermodynamic Partitioning ($K_{\text{part}}$):** Automated extraction of in situ dense-to-dilute partitioning coefficients ($C_{\text{dense}} / C_{\text{dilute}}$) with camera dark noise subtraction and active Z-slice background isolation.
- **2D Biophysical Phase Diagrams:** Automated classification of assemblies into distinct physical regimes (*Spherical LLPS Droplets*, *Chromatin-Wetted Condensates*, and *Hollow Membrane Structures*) based on $FA$ vs. $K_{\text{part}}$.
- **Radial FA Profiling (Shell–Middle–Core):** Sub-nanoscale geometric layer erosion using Euclidean distance transforms to identify flat lamellar sheets vs. isotropic spheres.
- **Topological QC Funnel:** Automated filtering for image/ROI edges, incomplete radial layers, and persistent Z-branching (`z_topology_fail`) to detect complex interconnected organelle networks (e.g. Golgi cisternae).
- **Physical ExM Calibration:** True biological volume scaling ($V_{\text{bio}} = V_{\text{pixel}} / \text{ExF}^3$) supporting $4\times$, $10\times$ (TREx), and custom expansion protocols.
- **Reproducible Batch Processing & Auditing:** Generates machine-readable CSVs, metadata JSON provenance, and publication-ready multi-sheet Excel reports.
- **Interactive Napari & ROI Previews:** Full 2D/3D Napari integration for manual ROI drawing and calibrated size distribution previews.

---

## Installation

An isolated Conda environment is recommended.

### 1. Clone the repository
```bash
git clone https://github.com/FrancinCZ/ExQt.git
cd ExQt
```

### 2. Create and activate the environment
```bash
conda env create -f environment.yml
conda activate exqt-env
```

Alternatively, install dependencies via pip:
```bash
python -m pip install -r requirements.txt
```

### 3. Start ExQt
```bash
python App.py
```

---

## Input Data Structure

For every raw TIFF stack, ExQt expects a segmentation mask in the same folder with the suffix `_Mask.tif`:

```text
sample_01.tif
sample_01_Mask.tif
```

- **Binary masks (`0/1` or `0/255`):** Connected components are labeled automatically.
- **Instance masks:** Unique positive integer per 3D object.

---

## Core Quantitative Modules

### 1. 3D Thermodynamic Partitioning ($K_{\text{part}}$)
Quantifies the thermodynamic propensity of molecules to concentrate into the condensed phase:
$$K_{\text{part}} = \frac{\max(I_{\text{obj}} - I_{\text{offset}}, 0)}{\max(I_{\text{nuc}} - I_{\text{offset}}, 10^{-6})}$$
*See [README_PARTITIONING.md](README_PARTITIONING.md) for full biophysical details, equations, and benchmarks.*

### 2. 2D Biophysical Phase Diagrams ($FA$ vs. $K_{\text{part}}$)
Separates unconstrained spherical droplets (**SON**, $FA = 0.533$) from chromatin-wetted condensates (**POL II**, $FA = 0.810$) and hollow membrane structures (**Golgi**, $K = 0.89\times$).

### 3. Radial FA Profiling (Shell–Middle–Core)
Measures how geometric fractional anisotropy changes from the outer surface to the innermost core:
$$\Delta A = \text{Core } FA - \text{Shell } FA$$
Positive $\Delta A$ with $\Delta I < 0$ uniquely identifies hollow 2D membrane sheets (e.g. Golgi cisternae).

---

## Experimental Benchmark Table

| Metric (3D Analysis) | 🔵 **POL II (Transcriptional Hubs)** | 🔴 **GOLGI / GM130 (Negative Control)** | 🟢 **SON (Nuclear Speckles)** |
| :--- | :---: | :---: | :---: |
| **Biological State** | Chromatin-tethered condensate | Hollow membrane labyrinth | Unconstrained spherical droplet |
| **Segmented Objects** | **$11\,900$** | **$207$** | **$96$** |
| **Median Intensity ($I_{\text{obj}}$)** | **$4\,724.72\text{ ADU}$** | **$679.46\text{ ADU}$** | **$5\,120.30\text{ ADU}$** |
| **Background ($I_{\text{nuc}}$)** | **$249.58\text{ ADU}$** | **$761.23\text{ ADU}$** | **$275.40\text{ ADU}$** |
| **Partitioning ($K_{\text{part}}$)** | **$\mathbf{18.93\times}$** | **$\mathbf{0.89\times}$** | **$\mathbf{18.57\times}$** |
| **Fractional Anisotropy ($FA$)** | **$0.810$** *(prolate ellipsoid)* | **$0.698\text{ to }0.99$** *(planar sheets)* | **$0.533$** *(globular spheroid)* |
| **Hollow Structures ($\Delta I < 0$)** | **$0\,\%$** | **$31.1\,\%$ in large cisternae** | **$0\,\%$** |

---

## Generated Output Files

Every run automatically outputs:
- `*_Output_Batch_3d.csv` — Full machine-readable table with all 3D morphological and partitioning parameters.
- `*_Output_Batch_3d_metadata.json` — Complete audit trail of applied calibration, thresholds, and software versions.
- `*_3d_partitioning_analysis.png` — *Size vs. $K_{\text{part}}$* and *2D Biophysical Phase Diagram*.
- `*_3d_Radial_FA_Profiling_Plots.png` — 4-panel *Shell–Middle–Core* progression and QC Funnel.
- `exm_partitioning_report_EN.docx` / `exm_partitioning_report_CZ.docx` — Ready-to-publish validation reports.

---

## License & Citation

ExQt is released under the [MIT License](LICENSE). See [About.md](About.md) for author information.

If you use ExQt in your research, please cite:
- **Cho, Spille, Cisse et al. (*Science* 2018):** *Mediator and RNA polymerase II clusters associate in transcription-dependent condensates.*
- **Feric, Pappu, Brangwynne et al. (*Cell* 2016):** *Coexisting Liquid Phases Underlie Nucleolar Subcompartments.*
