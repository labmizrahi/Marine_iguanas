Marine Iguana Microbiota Analysis
This repository contains the R scripts used for the ecological and evolutionary data analysis in the paper investigating the fecal microbiota of the Galapagos marine iguanas.

Citation
If you use this code or data, please cite the original publication:

Grinshpan, I., Lavy, O., Zorea, A., Amit, I., Levin, L., Furman, O., Somekh, D., Guevara, N., Moraïs, S., Cordero, O. X., & Mizrahi, I. (2026). Endemic within endemics: the microbiota of the Galapagos marine iguanas. ISME Communications, 6(1), ycag040. https://doi.org/10.1093/ismeco/ycag040

Repository Structure
The analysis is divided into three primary R scripts.

Reading and filtering data.R: Contains the code for importing the QIIME 2 pipeline outputs, subsampling reads to an even depth, and filtering erroneous ASVs using the Phyloseq package.

Neutral models.R: Includes the modified neutral community model (NCM) used to assess dispersal patterns across the metacommunity and partition taxa into ecologically selected or neutrally assembled groups.

Tree models and comparissons.R: Contains the code for phylogenetic tree construction, clustering analysis, comparative analysis with Fijian terrestrial iguanas, and calculations for the ASV Diversification Index.

Data Availability
Raw sequence reads, sample metadata, and environmental variables are deposited in the ENA under project PRJEB93862.
