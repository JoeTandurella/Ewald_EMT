# Psuedotime analysis of EMT genes
# Gaurav Sharma and Elana Fertig
# 2 Mar 2020

# R setup
library(monocle3)
library(dplyr)
library(Matrix)
library(pheatmap)
sessionInfo()

# R version 3.6.2 (2019-12-12)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Sierra 10.12.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] shiny_1.4.0                 pheatmap_1.0.12             Matrix_1.2-18               dplyr_0.8.3                
# [5] monocle3_0.2.0              SingleCellExperiment_1.6.0  SummarizedExperiment_1.14.1 DelayedArray_0.10.0        
# [9] BiocParallel_1.18.1         matrixStats_0.55.0          GenomicRanges_1.36.1        GenomeInfoDb_1.20.0        
# [13] IRanges_2.18.3              S4Vectors_0.22.1            Biobase_2.44.0              BiocGenerics_0.30.0        
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-6             RcppAnnoy_0.0.13         RColorBrewer_1.1-2       tools_3.6.2              DT_0.9                  
# [6] R6_2.4.0                 irlba_2.3.3              vipor_0.4.5              uwot_0.1.5               lazyeval_0.2.2          
# [11] colorspace_1.4-1         withr_2.1.2              tidyselect_0.2.5         gridExtra_2.3            compiler_3.6.2          
# [16] BiocNeighbors_1.2.0      labeling_0.3             scales_1.1.0             proxy_0.4-23             speedglm_0.3-2          
# [21] stringr_1.4.0            digest_0.6.21            rmarkdown_1.16           XVector_0.24.0           scater_1.12.2           
# [26] pkgconfig_2.0.3          htmltools_0.4.0          fastmap_1.0.1            htmlwidgets_1.5.1        rlang_0.4.4             
# [31] rstudioapi_0.10          DelayedMatrixStats_1.6.1 farver_2.0.3             jsonlite_1.6             RCurl_1.95-4.12         
# [36] magrittr_1.5             BiocSingular_1.0.0       GenomeInfoDbData_1.2.1   Rcpp_1.0.2               ggbeeswarm_0.6.0        
# [41] munsell_0.5.0            viridis_0.5.1            lifecycle_0.1.0          stringi_1.4.3            yaml_2.2.0              
# [46] MASS_7.3-51.4            zlibbioc_1.30.0          plyr_1.8.4               grid_3.6.2               promises_1.1.0          
# [51] ggrepel_0.8.1            crayon_1.3.4             lattice_0.20-38          splines_3.6.2            knitr_1.25              
# [56] pillar_1.4.2             igraph_1.2.4.1           reshape2_1.4.3           codetools_0.2-16         glue_1.3.1              
# [61] evaluate_0.14            leidenbase_0.1.0         RcppParallel_4.4.4       vctrs_0.2.3              batchelor_1.0.1         
# [66] httpuv_1.5.2             gtable_0.3.0             RANN_2.6.1               purrr_0.3.3              tidyr_1.0.0             
# [71] assertthat_0.2.1         ggplot2_3.2.1            xfun_0.10                rsvd_1.0.2               mime_0.7                
# [76] xtable_1.8-4             RSpectra_0.15-0          later_1.0.0              viridisLite_0.3.0        tibble_2.1.3            
# [81] beeswarm_0.2.3   

# Read data
cds <- readRDS('cdsBasic.RDS')
pData(cds)$media <- as.character(pData(cds)$media)
pData(cds)$media[is.na(pData(cds)$media)] <- 'CellLine'
pData(cds)$mediaday <- paste(as.character(pData(cds)$media),as.character(pData(cds)$day),sep="-")

# Normalize data
set.seed('1234')
cdsNorm <- cds
cdsNorm <- preprocess_cds(cdsNorm, num_dim = 50, method = 'PCA')
cdsNorm <- reduce_dimension(cdsNorm, reduction_method="UMAP")
cdsNorm = align_cds(cdsNorm, num_dim = 50, alignment_group = "day")
cdsNorm <- reduce_dimension(cdsNorm, reduction_method="UMAP")

pdf('mediadayumap_2Mar2010.pdf')
plot_cells(cdsNorm, color_cells_by="mediaday")
dev.off()

# Psuedotime analysis for Collagen

cdsNormCollagen <- cdsNorm[,pData(cdsNorm)$media=='Collagen']

cdsNormCollagen <- cluster_cells(cdsNormCollagen)
cdsNormCollagen <- learn_graph(cdsNormCollagen)

cdsNormCollagen <- order_cells(cdsNormCollagen)

pdf('collagenumaps_2mar2020.pdf')
plot_cells(cdsNormCollagen, color_cells_by = 'day')
plot_cells(cdsNormCollagen, color_cells_by = 'pseudotime')
dev.off()

EMTGenes <- read.table('geneList.csv', stringsAsFactors = F)[,1]

pdf('collagenEMTGenes_2mar2020.pdf')
plot_genes_in_pseudotime(cdsNormCollagen[c(EMTGenes[1:6],'Vim'),], color_cells_by = 'day')
dev.off()



# Psuedotime analysis for Matrigel
cdsNormMatrigel <- cdsNorm[,pData(cdsNorm)$media=='Matrigel']

cdsNormMatrigel <- cluster_cells(cdsNormMatrigel)
cdsNormMatrigel <- learn_graph(cdsNormMatrigel)

cdsNormMatrigel <- order_cells(cdsNormMatrigel)

pdf('matrigelumaps_2mar2020.pdf')
plot_cells(cdsNormMatrigel, color_cells_by = 'day')
plot_cells(cdsNormMatrigel, color_cells_by = 'pseudotime')
dev.off()

pdf('matrigelEMTGenes_2mar2020.pdf')
plot_genes_in_pseudotime(cdsNormMatrigel[c(EMTGenes[1:6],'Vim'),], color_cells_by = 'day')
dev.off()

# Ouptut CDS for subsequent analysis
save(list=c('cdsNorm','cdsNormMatrigel','cdsNormCollagen','EMTGenes'), file='EJFPsuedotimeResults2Mar2020.Rda')

