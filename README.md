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
├── Code/
│   ├── deduplicate.R
│   ├── metaDigitiser.R
│   ├── meta_analysis.R
│   └── study_characteristics.R
├── Data/
│   ├── cleaned_data_for_meta_analysis.csv
│   ├── ulva_meta_analysis_raw_data.csv
│   ├── metadata.csv
│   └── figures_for_extraction/
│       ├── Boarder/
│       │   ├── Figure_9a.jpg
│       │   └── Figure_9b.jpg
│       ├── Falade/
│       │   ├── S006.1_CF.jpg
│       │   ├── S006.1_FinalWeight.jpg
│       │   ├── S006.1_SGR.jpg
│       │   ├── S006.2_condition_factor.jpg
│       │   ├── S006.2_FCR.jpg
│       │   ├── S006.2_FeedIntake.jpg
│       │   ├── S006.2_final_shell_length.jpg
│       │   ├── S006.2_PER.jpg
│       │   └── S006.2_SGR.jpg
│       └── Mwangudza/
│           ├── S005.1_Final_length.jpg
│           └── S005.1_Final_weight.jpg
├── Figures/
│   ├── all_data_mlmr_plot.png          # meta_analysis.R
│   ├── all_data_orchard_plot.png       # meta_analysis.R
│   ├── bubble_plot_species.png         # study_characteristics.R
│   ├── combined_corr.png               # meta_analysis.R
│   ├── combined_fb.png                 # meta_analysis.R
│   ├── combined_gp.png                 # meta_analysis.R
│   ├── combined_nutr.png               # meta_analysis.R
│   ├── combined_plot.png               # meta_analysis.R
│   ├── combined_plot_sens.png          # meta_analysis.R
│   ├── decision_tree.png               # [source script not identified]
│   ├── funnel_plot_all.png             # meta_analysis.R
│   ├── funnel_plot_sens.png            # meta_analysis.R
│   ├── journal_year.png                # study_characteristics.R
│   ├── PRISMA_MA.png                   # [source script not identified]
│   ├── Sankey_plot.png                 # study_characteristics.R
│   ├── sens_data_orchard_plot.png      # meta_analysis.R
│   ├── sens_mlmr_plot.png              # meta_analysis.R
│   ├── ulva_inclusion_relationship.png # meta_analysis.R
│   └── world_map.png                   # study_characteristics.R
├── Screening/
│   [screening records to be added]
├── Searches/
│   ├── grey_literature/
│   │   ├── BASE_ALL_10102025.csv
│   │   └── BASE_duplicate_removed_10102025.csv
│   └── primary_literature/
│       ├── SCOPUS_11092025.csv
│       ├── WOS_11092025.xls
│       ├── WOS_and_SCOPUS_combined_11092025.csv
│       └── WOS_and_SCOPUS_duplicates_removed_11092025.csv
└── README.md
```

## Code

### `code/deduplicate.R`

R script used to identify and remove duplicate records from the combined Web of Science and SCOPUS search export (`WOS_and_SCOPUS_combined_11092025.csv`). Outputs the deduplicated file (`WOS_and_SCOPUS_duplicate_removed_11092025.csv`).

### `code/metaDigitiser.R`

R script used to extract numerical data from figures in the `data/figures_for_extraction/` subfolder using the `metaDigitiser` package. Produces extracted data suitable for inclusion in the primary analysis dataset.

### `code/study_characteristics.R`

R script used to clean data, summarise and visualise study characteristics across the full primary dataset, including species, *Ulva* inclusion levels, experimental duration, and geographic distribution of included studies.

### `code/meta_analysis.R`

R script implementing the full MLMA, subsequent MLMR, publication bias testing, and figure generation.

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

