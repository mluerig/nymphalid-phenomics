## nymphalid-phenomics

</div><br>
<div align="center">
	<p><img src="assets/composite_dorsal1.png" width="90%"></p>
</div>

### About
This repository contains the analytical pipeline and code to reproduce the results of the paper "Aposematic color patterns are the dominant axis of phenotypic diversification in Nymphalid butterflies" (Lürig et al. 2025). The pipeline combines image acquisition (mostly from online sources), segmentation, feature extraction, and phylogenetically informed statistical analysis to quantify aposematic color patterns in over one third of all species of the butterfly family Nymphalidae. 

### Quickstart
1. **Clone repo and install environment:**
	```bash
	git clone https://github.com/mluerig/nymphalid-phenomics.git
	cd nymphalid-phenomics
	mamba env create -f environment.yml -n nymphalidae1
	conda activate nymphalidae1
	```
2. **Download the data:**
	- Download the archived repository that contains all images and metadata from [Zenodo](https://doi.org/10.5281/zenodo.17204330). 
	- Either place the segmented images in `data_raw/images/` and meta-data tables in `data/data_primary/`, or simply `cd` into the downloaded Zenodo repo.
3. **Run the analysis scripts:**
	- See Analysis section for details.

### Installation
1. **Create the environment:**
	```bash
	mamba env create -f environment.yml -n nymphalidae1
	conda activate nymphalidae1
	```
2. **Install additional packages:**

	UNICOM (image encoder):
    ```bash
    pip install --no-cache-dir "unicom @ git+https://github.com/deepglint/unicom.git@4d84a3b496a47bcad68467d71c5ca787b0366042"
    ```

	PyTorch (choose wheel matching your CUDA - may need the `--force-reinstall` flag):
    ```bash
    pip install --index-url https://download.pytorch.org/whl/cu126 torch torchvision
    ```

### Data

**Download:**
  - A sample of raw images, all segmentation masks, as well as primary tabular- and meta-data are available through Zenodo: [https://doi.org/10.5281/zenodo.17204330](https://doi.org/10.5281/zenodo.17204330)
  - The Zenodo-repo is simply a fully populated version of this GitHub-repo. So, you can either place the segmented images in `data_raw/images/` and meta-data tables in `data/data_primary/`, or just `cd` into the downloaded Zenodo repo and then run the scripts.

**Structure:**
- `data_raw/`
    - `images_sample/` — Sample of raw images (from GBIF)
    - `segmentation_masks_clean/` — Segmented masks (produced by pipeline)
    - `segmentation_masks_moths/` — Segmented masks (moth dataset)
    - `tables/` — Embeddings and features (produced by pipeline)
- `data/`
    - `data_primary/` — Primary tabular and meta-data (labels, feature-key, etc.)
    - `data_secondary/` — Derived from primary with make_data script (LD-scores, similarity, etc.)
    - `analyses_secondary/` — Regressions, phylogenetic modelling, etc.
- `scripts/` - scripts to reproduce all results.
- `figures/` - figures from the manuscript.
- `tables/` - tables from the manuscript.


### Analysis

</div><br>
<div align="center">
    <p><img src="assets/figureS1.png" width="400"></p>
</div>

For the paper, we ran the full pipeline, which included 1) image downloads from GBIF and idigbio, and own imaging, 2) segmentation (see script 01)), 3) image encoding (see script 02), 4) cleaning, 5) wing surface classification, 6) literature review, and 7) statistical analysis (see scripts 03 and 04). We do not provide the code for cleaning and classifier training, as it was an iterative process with several manual steps (e.g., selection from interactive plots). The full procedure is described in detail in the methods-section of the manuscript. 

The following scripts will allow you to reproduce the results:

**In Python:**
- [`01_segmentation_demo.ipynb`](scripts/01_segmentation_demo.ipynb) — segments butterflies from images using GroundedSAM (demo, without cleaning-steps).
- [`02_feature_extraction.ipynb`](scripts/02_feature_extraction.ipynb) — extracts embeddings using UNICOM.

**In R:**
- [`03_make_data.R`](scripts/03_make_data.R) — assembles specimen-level and species-level tables for analysis.
- [`04_analysis.R`](scripts/04_analysis.R) — runs all statistical models and generates figures and tables for the paper.


### Citation
