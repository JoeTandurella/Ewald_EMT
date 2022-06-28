# EMT_Ewald
Curated set of analysis scripts for single cell data from the Ewald Lab.

## Publication Title
Triple negative breast cancer metastasis involves complex epithelial-mesenchymal transition dynamics and requires vimentin

Authors: Eloïse M. Grasset1*, Matthew Dunworth1, Gaurav Sharma2, Melanie Loth3, Joseph Tandurella3, Ashley Cimino-Mathews3,4, Melissa Gentz1, Sydney Bracht1, Meagan Haynes1, Elana J. Fertig2,3,5,6, Andrew J. Ewald1,2,3,6*

Affiliations:
1 Department of Cell Biology, School of Medicine, Johns Hopkins University, Baltimore, MD 21205, USA
2 Department of Biomedical Engineering, School of Medicine, Johns Hopkins University, Baltimore, MD 21205, USA
3 Department of Oncology, Sidney Kimmel Comprehensive Cancer Center, Johns Hopkins University School of Medicine, Baltimore, MD 21205, USA
4 Department of Pathology, School of Medicine, Johns Hopkins University, Baltimore, MD 21205, USA
5 Department of Applied Mathematics and Statistics, Johns Hopkins University, Baltimore, MD 21205, USA
6 Convergence Institute, Johns Hopkins University, Baltimore, MD 21205, USA

### 1_CombineData.ipynb
#### Description - Joins the counts data together
##### Input 
* D1-D3-D5-byHaiping/D1/filtered_feature_bc_matrix/matrix.mtx
* D1-D3-D5-byHaiping/D1/filtered_feature_bc_matrix/barcodes.tsv
* D1-D3-D5-byHaiping/D1/filtered_feature_bc_matrix/features.tsv
* D1-D3-D5-byHaiping/D3/filtered_feature_bc_matrix/matrix.mtx
* D1-D3-D5-byHaiping/D3/filtered_feature_bc_matrix/barcodes.tsv
* D1-D3-D5-byHaiping/D3/filtered_feature_bc_matrix/features.tsv
* D1-D3-D5-byHaiping/D5/filtered_feature_bc_matrix/matrix.mtx
* D1-D3-D5-byHaiping/D5/filtered_feature_bc_matrix/barcodes.tsv
* D1-D3-D5-byHaiping/D5/filtered_feature_bc_matrix/features.tsv
##### Output
* mergedAdata.h5ad

### 2_cleanData.ipynb
#### Description - Preprocesses the combined data
##### Input 
* mergedAdata.h5ad
##### Output
* labelledAdata.h5ad
* data/dataMatrix.mtx'
* data/obs.csv
* data/var.csv
* rawCleanadata.h5ad
* data/dataMatrixClean.mtx
* data/obsClean.csv
* data/varClean.csv

### 3_ClusteringAfterBatchCorrection.Rmd
#### Description - Create CDS with the “clean” data. CDS is preprocessed.
##### Input 
* data/dataMatrixClean.mtx
* data/obsClean.csv
* data/varClean.csv
##### Output
* intermediateData/cdsBasic.RDS
* intermediateData/cdsBasic.RDS
* data/normalizedMatrix.mtx
* intermediateData/cdsNormalizedDimReduced.rds
* intermediateData/cdsDayCorrected.RDS

### 4_EJFPreprocessingAnalysis.Rmd
#### Description - Run pseudotime for matrigel and collagen
##### Input 
* intermediateData/cdsBasic.RDS
* geneList.csv
##### No Output

### 5_EJFPreprocessingAnalysis.R
#### Description - Run pseudotime for matrigel and collagen
##### Input 
* cdsBasic.RDS  
* geneList.csv
##### Output
* EJFPsuedotimeResults2Mar2020.Rda

### 6_PseudotimeCorrectedData.Rmd
#### Description - Process and plot the corrected data
##### Input 
* intermediateData/cdsDayCorrected.RDS
* geneList.csv
##### No Output

### 7_ScanpyEMT.ipynb
#### Description - Calculate and export cell cycle 
##### Input 
* rawCleanadata.h5ad
* regev_lab_cellCycle_smouse.txt
* regev_lab_cellCycle_mmouse.txt
##### Output
* obs.csv

### 8_addCellCycle.Rmd
#### Description - Add cell cycle to cds objects
##### Input 
* EJFPsuedotimeResults2Mar2020.Rda
* obs.csv
##### Output
* cellCycle_EJFPsuedotimeResults2Mar2020.Rda

### 9_CoGAPSData_collagen.Rmd
#### Description - Calculate top collagen genes and subset the data to those genes for CoGAPS
##### Input 
* EJFPsuedotimeResults2Mar2020.Rda
##### Output
* sbstCds.rds
* countsSbstCds.rds

### 10_cogapsObjects_collagen.Rmd
#### Description - Create CoGAPS objects for the collagen gene subset
##### Input 
* countsSbstCds.rds
##### Output
* emtMat.mtx
* emtParams.rds

### 11_ewald_gwCogapsAnalysis_collagen.Rmd
#### Description - Analyze CoGAPS results
##### Input 
* sbstCds.rds
* CoGAPS results
##### Output
* UMAPS
* Boxplots
* patternMarkers_emtGwCogaps_20.csv

### 12_cogapsVisuals_collagen.Rmd
#### Description - Create heatmap and gene rank file for CoGAPS result
##### Input 
* gwCogapsRes/emtGwCogaps_20.rds
* geneList.csv
##### Output
* emtGeneRanks.csv

### Ran gene set analysis using MSigDb
#### Input
* Pattern marker genes
#### Output
* Pathway results

### 13_Figure.Rmd
#### Description - Create the first set of figures
##### Input 
* cellCycle_EJFPsuedotimeResults2Mar2020.Rda
* gwCogapsRes/emtGwCogaps_20.rds
* msigdb/emtPatterns.csv
* patternMarkers_emtGwCogaps_20.csv
* geneList.csv
* msigdb/emtPatterns.csv
* emtResults.csv
##### Output
* Figures

