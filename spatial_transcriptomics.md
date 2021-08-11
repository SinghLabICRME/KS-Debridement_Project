## Project: Genome-Wide Epigenetic Gene Silencing at the Wound-Edge of Chronic Wound Patients Impairs Epithelial to Mesenchymal Transition Opposing Closure
#### Running title: " Spatial transcriptomics profile of signature genes for Kera1 and Kera2 clusters in human skin"
#### Author: "Rajneesh Srivastava"
#### Date: "06/01/2021"

### 1. Setup the directory for data analysis in R (version 4.0.4)
```setwd ("/path/SEN_visium_dataset/Result/")```

### 2. Software installation
```
install.packages("Seurat")   # version 4.0.2
install.packages("hdf5r")
install.packages("ggplot2")  # version 3.3.4
install.packages("future")
```
### 3. Load libraries in R
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

### 4. Upload Data
Data were submitted in Gene Expression Omnibus (GEO) repository with GSE ID: **GSE176417**.
Spatial transcriptomics data can be accessed/ downloaded via below link:

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5364340

Alternatively, processed .rds file can be accessed via below link and ready to use from **step 8**
https://drive.google.com/file/d/10304LEvPBRgDnhBew8U6LzYJ0rNUTepS/view?usp=sharing

Note - This analysis require h5 (or split files; GSM5364340) file along with image files containing histology pics (most likely named "spatial.tar.gz" or "spatial" folder)

```
KS.dir = "/path/"                                             
KS_data <- Load10X_Spatial(data.dir = paste0(KS.dir, NS9),
						filename = "filtered_feature_bc_matrix.h5", 
						assay = "Spatial", 
						slice = NS9, 
						filter.matrix = TRUE)
            assign (NS009, KS_data)
```
### 5. Compute mitochondrial content and filter
```
NS9[["percent.mt"]] <- PercentageFeatureSet(NS9, pattern = "^MT-")
NS9 <- subset(x = NS9, subset = percent.mt < 15)
```
### 6. Data transformation and scaling using 'SCT' module in Seurat
```
NS9 <- SCTransform(NS9, 
		assay = "Spatial", 
		vars.to.regress = "percent.mt",
		verbose = FALSE)
```
### 7. Check if the sample meets the analysis criteria (Please see methods)
```
VlnPlot(NS9, 
	features = c("nCount_Spatial",
			"nFeature_Spatial",     
			"percent.mt"), 
			 ncol = 3, 
			 pt.size=0)
```
### 8. Save RDS file
```
setwd ("/path/")
saveRDS (NS9, file = "NS9.rds")
NS9=readRDS(file="./path/NS9.rds")
```
### 9. Spatial transcriptomic profile of signature genes for Kera1 and Kera 2 clusters.
```
Kera1=c("KRT14","KRT1")
Kera2=c("KRT19","KRT7")
SpatialFeaturePlot(NS9, features = Kera1)
SpatialFeaturePlot(NS9, features = Kera2)
```
### 10. Spatial transcriptome profile (and localization) of top transcription factors,enzymes and metabolically active genes with increased expression levels in Kera 2
```
Genes=c("GAPDH", "CITED4", "COX7C", "COX8A", "COX5B", "NDUFA4", "COX7B", "PRDX5", "COX6A1")
SpatialFeaturePlot(NS9, features = Genes)
```
## Thank you
