## Project: Genome-Wide Epigenetic Gene Silencing at the Wound-Edge of Chronic Wound Patients Impairs Epithelial to Mesenchymal Transition Opposing Closure
### Running title: " Spatial transcriptomics profile of signature genes for Kera1 and Kera2 clusters in human skin"
#### Author: "Rajneesh Srivastava"
#### Date: "06/01/2021"

### Set directory for data analysis
```setwd ("/path/SEN_visium_dataset/Result/")```

### Load libraries in R
```
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(hdf5r)
library(patchwork)
library(future)
library(sctransform)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 15000 * 1024^2)
```

### Upload Data
Note - require h5 file along with image files containing histology pics (most likely named "spatial.tar.gz" or "spatial" folder)
#### Data uploaded in GEO with GSE ID
```
KS.dir = "/path/"                                             
KS_data <- Load10X_Spatial(data.dir = paste0(KS.dir, NS009),
						filename = "filtered_feature_bc_matrix.h5", 
						assay = "Spatial", 
						slice = NS009, 
						filter.matrix = TRUE)
            assign (NS009, KS_data)
```
## Compute mitochondrial content and filter
```
NS009[["percent.mt"]] <- PercentageFeatureSet(NS009, pattern = "^MT-")
NS009 <- subset(x = NS009, subset = percent.mt < 15)
```
## Data transformation and scaling using 'SCT' module in Seurat
```
NS009 <- SCTransform(NS009, 
						assay = "Spatial", 
						vars.to.regress = "percent.mt",
						verbose = FALSE)
```
## Check if the sample meets the analysis criteria (Please see methods)
```
VlnPlot(NS009, 
		features = c("nCount_Spatial",
					 "nFeature_Spatial",     
					  "percent.mt"), 
					  ncol = 3, pt.size=0)
```
## Spatial transcriptomic profile of signature genes for Kera1 and Kera 2 clusters.
```
Kera1=c("KRT14","KRT1")
Kera2=c("KRT19","KRT7")
SpatialFeaturePlot(NS009, features = Kera1)
SpatialFeaturePlot(NS009, features = Kera2)
```
## Spatial transcriptome profile (and localization) of top transcription factors,enzymes and metabolically active genes with increased expression levels in Kera 2
```
Genes=c("GAPDH", "CITED4", "COX7C", "COX8A", "COX5B", "NDUFA4", "COX7B", "PRDX5", "COX6A1")
SpatialFeaturePlot(NS009, features = Genes)
```
## Save RDS file
```
setwd ("/path/")
saveRDS (NS009, file = "NS009.rds")
```
## Thank you
