Feeding behaviour, growth performance and nutrient utilisation of abalone with dietary Ulva sp. supplementation: A meta-analysis 

This repository stores the data, bibliometric files, and code used for this study. Please find the description for each folder below. Kindly contact Rebecca Pedler (Rebecca.pedler@yumbah.com) for any queries, or to reuse any data or analysis for future studies.

Code

This folder contains all code used to obtain results for the meta-analysis:

•	metaDigitiser. This R script was used to extract data from figures in the figures_for_extraction subfolder found in Data.

•	all_outcomes. This R script was used to conduct the MLMA on the entire dataset.

•	deduplicate. This R script was used to remove duplicate hits from the WOS_and_SCOPUS_combined_11092025.csv file.

•	feed_intake. This R script was used to conduct the MLMA, subsequent MLMR, publication bias testing and figure generation for the dataset where outcome_category in Ulva inclusion in Haliotis sp. diets_ A Meta-analysis.csv is feed_behaviour.

•	growth_performance. This R script was used to conduct the MLMA, subsequent MLMR, publication bias testing and figure generation for the dataset where outcome_category in Ulva inclusion in Haliotis sp. diets_ A Meta-analysis.csv is growth_performance.

•	nutrient_utilisation. This R script was used to conduct the MLMA, subsequent MLMR, publication bias testing and figure generation for the dataset where outcome_category in Ulva inclusion in Haliotis sp. diets_ A Meta-analysis.csv is nutrient_utilisation.

•	orchard_plots. This R script was used to generate orchard plots for the MLMA on all data as well as feed behaviour, growth performance and nutrient utilisation datasets.

•	study_characteristics. This R script was used to assess study characteristics of the entire Ulva inclusion in Haliotis sp. diets_ A Meta-analysis.csv dataset.

Data

This folder contains data and metadata for the meta-analysis:

•	Ulva inclusion in Haliotis sp. diets_ A Meta-analysis.csv: This csv contains extracted data from all eligible articles.

•	Ulva inclusion in Haliotis sp. diets_ A Meta-analysis - Metadata.csv: This csv contains the metadata description of information extracted from eligible articles.

•	sub-folder figures_for_extraction which houses figures used for data extraction. Data was extracted using the MetaDigitiser.R code. 

Within this sub-folder, the following is stored:   

•	Falade. This contains all figures from the thesis:

Falade AE. Optimising Integrated Multitrophic Aquaculture (IMTA) on a South African Abalone Farm. Rhodes University; 2023.

-	S006.1_CF.jpg
-	S006.1_FinalWeight.jpg
-	S006.1_SGR.jpg
-	S006.2_FCR.jpg
-	S006.2_FeedIntake.jpg
-	S006.2_PER.jpg
-	S006.2_SGR.jpg
-	S006.2_condition_factor.jpg
-	S006.2_final_shell_length.jpg

•	Mwangudza. This contains all figures from the thesis:
Mwangudza PM. Assessment and Mitigation of Biosecurity Risks Associated with Macroalgae Inclusion in Farmed Abalone Diets in South Africa. Rhodes University; 2024.
-	S005.1_Final_length.jpg
-	S005.1_Final_weight.jpg

Searches

This folder contains bibliometric files created during primary and grey literature searches. These files are stored in the following subfolders:

primary_literature: within this folder, the following files are stored:
•	WOS_11092025.csv: This file contains all hits returned from the search string on Web of Science Core Collection (11092025).

•	SCOPUS_11092025.csv: This file contains all hits returned from the search string on SCOPUS (11092025).

•	WOS_and_SCOPUS_combined_11092025.csv: This file contains the combined hits returned from Web of Science Core Collection and SCOPUS (11092025).

•	WOS_and_SCOPUS_duplicate_removed_11092025.csv: This file contains the combined hits returned from Web of Science Core Collection and SCOPUS (11092025) after deduplication using Deduplicate.

Grey_literature: within this folder, the following files are stored:
•	BASE_ALL_10102025.csv: This file contains all hits returned from BASES (10102025) using a modified string “Haliot” AND “Diet”.

•	BASE_duplicate_removed_10102025.csv: This file contains all hits returned from BASE (10102025) using a modified string “Haliot” AND “Diet” and after deduplication using Deduplicate.R.

Screening This folder contains abstract and full-text screening files:
•	RP_abstracts.csv: This csv contains bibliometric information and outcomes for articles proceeding through abstract screening.

•	RP_fulltext.csv: This csv contains bibliometric information and outcomes for articles proceeding through full-text screening.


