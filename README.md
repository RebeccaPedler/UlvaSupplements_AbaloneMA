# UlvaSupplements_AbaloneMA

Code, data, and outputs supporting the manuscript:

**"Feeding behaviour, growth performance and nutrient utilisation of abalone with dietary *Ulva* sp. supplementation: A systematic review and multi-level meta-analysis"**

Rebecca L. Pedler, Matthew S. Bansemer, James O. Harris, Ondi L. Crino

---

## Overview

This repository contains the full analytical pipeline, raw data, bibliometric search files, and screening records used to conduct a systematic review and meta-analysis of *Ulva* sp. supplementation in abalone (*Haliotis* spp.) diets. 

## Repository structure

```
UlvaSupplements_AbaloneMA/
├── code/
│   ├── metaDigitiser.R
│   ├── all_outcomes.R
│   ├── deduplicate.R
│   ├── feed_intake.R
│   ├── growth_performance.R
│   ├── nutrient_utilisation.R
│   ├── orchard_plots.R
│   └── study_characteristics.R
├── data/
│   ├── Ulva inclusion in Haliotis sp. diets_ A Meta-analysis.csv
│   ├── Ulva inclusion in Haliotis sp. diets_ A Meta-analysis - Metadata.csv
│   └── figures_for_extraction/
│       ├── Falade/
│       └── Mwangudza/
├── searches/
│   ├── primary_literature/
│   │   ├── WOS_11092025.csv
│   │   ├── SCOPUS_11092025.csv
│   │   ├── WOS_and_SCOPUS_combined_11092025.csv
│   │   └── WOS_and_SCOPUS_duplicate_removed_11092025.csv
│   └── grey_literature/
│       ├── BASE_ALL_10102025.csv
│       └── BASE_duplicate_removed_10102025.csv
└── screening/
    ├── RP_abstracts.csv
    └── RP_fulltext.csv
```

## Code

### `code/deduplicate.R`

R script used to identify and remove duplicate records from the combined Web of Science and SCOPUS search export (`WOS_and_SCOPUS_combined_11092025.csv`). Outputs the deduplicated file (`WOS_and_SCOPUS_duplicate_removed_11092025.csv`).

### `code/metaDigitiser.R`

R script used to extract numerical data from figures in the `data/figures_for_extraction/` subfolder using the `metaDigitiser` package. Produces extracted data suitable for inclusion in the primary analysis dataset.

### `code/meta_analysis.R`

R script implementing the full MLMA, subsequent MLMR, publication bias testing, and figure generation.

### `code/study_characteristics.R`

R script used to summarise and visualise study characteristics across the full primary dataset, including species, *Ulva* inclusion levels, experimental duration, and geographic distribution of included studies.

## Data

### `data/Ulva inclusion in Haliotis sp. diets_ A Meta-analysis.csv`

The primary analysis dataset. Contains extracted data from all eligible articles.

### `data/Ulva inclusion in Haliotis sp. diets_ A Meta-analysis - Metadata.csv`

Metadata description for the primary dataset.

### `data/figures_for_extraction/`

Figures used for numerical data extraction via MetaDigitiser. Organised by study, with the following sources:

**Falade** — Figures from:
Falade AE. *Optimising Integrated Multitrophic Aquaculture (IMTA) on a South African Abalone Farm*. Rhodes University; 2023.

**Mwangudza** — Figures from:
Mwangudza PM. *Assessment and Mitigation of Biosecurity Risks Associated with Macroalgae Inclusion in Farmed Abalone Diets in South Africa*. Rhodes University; 2024.

## Searches

### `searches/primary_literature/`

Bibliometric search exports from Web of Science Core Collection and SCOPUS (conducted 11 September 2025), along with the combined file and the deduplicated output used for screening.

### `searches/grey_literature/`

Search exports from BASE (Bielefeld Academic Search Engine; conducted 10 October 2025), along with the deduplicated output.

## Screening

### `screening/RP_abstracts.csv`

Bibliometric information and eligibility outcomes for all records proceeding through abstract screening.

### `screening/RP_fulltext.csv`

Bibliometric information and eligibility outcomes for all records proceeding through full-text screening.

## Citation

Pedler RL, Bansemer MS, Harris JO, Crino OL. Feeding behaviour, growth performance and nutrient utilisation of abalone with dietary *Ulva* sp. supplementation: A systematic review and multi-level meta-analysis. *In preparation*.

## Contact

Rebecca Pedler
rebecca.pedler@yumbah.com | Rebecca.pedler@flinders.edu.au

