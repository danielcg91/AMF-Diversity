# AMF-Diversity
This repository contains the code used to perform the upstream and downstream statistical analyses in Chavez-Gonzalez, et al. 2023 which was modified from pipeline in Flores-Nuñes et al., 2023
AMFData.bash : Pipeline that processes 16SrRNAV4 and ITS2 amplicon sequencing paired reads to generate OTUs, OTU table and Taxonomic classification. see citation on manuscript. Also, contains the functions to make the reverse complement for 16S sequences.

1-Parse.output-AMF.R: ./data : includes 16S and ITS filtered OTUs, OTU table, taxonomic classification, and metadata

2.-getOTUdata-ITS2-16S.R: R script that generates filtered OTU, taxa and metadata tables for future analyses. 

3.-AMFdata-root.colonization.R: R script and test analyses with the data generated in Chavez-Gonzalez, et al. 2023 
