# Symbioses between desert plants, arbuscular mycorrhizal fungi, and endofungal bacteria promote drought tolerance in arid ecosystems
This repository contains the code used to perform the upstream and downstream statistical analyses in Chavez-Gonzalez, et al. 2023 which was modified from the pipeline in Flores-Nuñes et al., 2023

0.-AMFData.bash : Pipeline that processes 16SrRNAV4 and ITS2 amplicon sequencing paired reads to generate OTUs, OTU table, and Taxonomic classification. see a citation on the manuscript. Also, contains the functions to make the reverse complement for 16S sequences.

1-Parse.output-AMF.R: ./data : includes 16S and ITS filtered OTUs, OTU table, taxonomic classification, and metadata

2.-getOTUdata-ITS2-16S.R: R script that generates filtered OTU, taxa and metadata tables for future analyses. 

3.-AMFdata-root.colonization.R: R script and test analyses with the data generated in Chavez-Gonzalez, et al. 2023 

4.-Barplots-VennDiagramsT12-AMF-16Spore.R and Barplots-VennDiagramsT12-ITS-16S.R: R script that generates all barplots of taxonomic classification, VennDiagrams and all packages used. 


