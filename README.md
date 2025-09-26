## nymphalid-phenomics

</div><br>
<div align="center">
    <p><img src="assets/composite_dorsal1.png"></p>
</div>

### About
This repository contains the full pipeline and code for the analysis in the paper on global color pattern variation and chemical defense in Nymphalidae butterflies. The project combines large-scale image acquisition, segmentation, feature extraction, and statistical analysis to quantify aposematic signals and similarity metrics across thousands of species. All code, data, and instructions are provided for full reproducibility.

### Quickstart
1. **Clone the repository:**
	```bash
	git clone https://github.com/mluerig/nymphalid-phenomics.git
	cd nymphalid-phenomics
	```
2. **Download the data:**
	- Download the full image and metadata dataset from Zenodo (see Data section below).
	- Place raw images in `data_raw/images/` and metadata tables in `data_raw/tables/`.
3. **Run the analysis scripts:**
	- See Analysis section for details.

### Installation
1. **Create the environment:**
	```bash
	mamba env create -f environment.yml -n nymphalidae1
	conda activate nymphalidae1
	```
2. **Install additional packages:**
	- UNICOM (image encoder):
	  ```bash
	  pip install --no-cache-dir "unicom @ git+https://github.com/deepglint/unicom.git@4d84a3b496a47bcad68467d71c5ca787b0366042"
	  ```
	- PyTorch (choose wheel matching your CUDA):
	  ```bash
	  pip install --index-url https://download.pytorch.org/whl/cu126 torch torchvision
	  ```
	- Or use mamba for PyTorch:
	  ```bash
	  mamba install pytorch torchvision -c pytorch -c nvidia
	  ```

### Data

- **Download:**
  - All raw images, segmentation masks, and metadata tables are available from Zenodo:
	 - [Zenodo DOI link here — add actual link]
  - Place images in `data_raw/images/` and tables in `data_raw/tables/`.
- **Structure:**
  - `data_raw/images/` — raw specimen images
  - `data_raw/segmentation_masks/` — segmented masks (produced by pipeline)
  - `data_raw/tables/` — metadata, results, and analysis tables

### Analysis

</div><br>
<div align="center">
    <p><img src="assets/figureS1.png" width="500"></p>
</div>

Run the following scripts in order to reproduce the results:
1. **Segmentation:**
	- `scripts/01_segmentation_demo.ipynb` — segments butterflies from images using GroundedSAM.
2. **Feature extraction & cleaning:**
	- `scripts/01_feature_extraction.ipynb` — extracts UNICOM embeddings, cleans masks, and prepares features.
3. **Data assembly:**
	- `scripts/02_make_data.R` — assembles specimen-level and species-level tables for analysis.
4. **Statistical analysis:**
	- `scripts/03_analysis.R` — runs all statistical models and generates figures and tables for the paper.

### Citation
If you use this code or data, please cite:

> Luerig, M. et al. (2025). Global color pattern variation and chemical defense in Nymphalidae butterflies. [Journal Name, Volume, Pages]. DOI: [add DOI]

Zenodo data DOI: [add Zenodo DOI]