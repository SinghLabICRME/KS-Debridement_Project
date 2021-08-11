# KS-Debridement_Project
## Project: Genome-Wide Epigenetic Gene Silencing at the Wound-Edge of Chronic Wound Patients Impairs Epithelial to Mesenchymal Transition Opposing Closure
### Running title: " scRNA-seq transcriptomics profile of normal skin and chronic wound"
#### Author: "Ahmed S Abouhashem"
#### Date: "06/01/2020"

### Set directory for data analysis
```setwd ("/path/SEN_scRNA-seq_dataset/Result/")```

### Load libraries in R
```
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
```

### Upload Data
Note - require 3 files for each sample (matrix, barcodes and features) or directly upload the rds file: 
https://drive.google.com/file/d/18NM8KxyTua7dbk6CvZMz9jKhd4wFo0ak/view?usp=sharing)
# If proceeding with .rds file
## reading the Seurat object and set the object identity to seurat_clusters
```
object = readRDS('Insert_Folder_Path_here/object.rds')
Idents(object)='seurat_clusters'
```

## plot total number of counts, genes and mitochondrial gene expression in cells
```
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
```

## tSNE plot
```
DimPlot(object, label = T, split.by = 'condition')
```

## Figure 3B: subsetting normal skin from the main Seurat object and label keratinocyte clusters/ subsetting keratinocyte clusters and plotting cells
```
normal_skin = subset(object, subset=condition=='Normal skin')
kera = subset(normal_skin, subset = seurat_clusters == 5|seurat_clusters == 6)
kera <- FindVariableFeatures(object = kera, selection.method = "vst", nfeatures = 2000)
kera <- ScaleData(kera, verbose = TRUE)
kera <- RunPCA(kera, npcs = 50, verbose = FALSE)
kera = RunTSNE(kera, dims = 1:10)
DimPlot(normal_skin, cols = c('grey','grey','grey','grey','grey','black','black','grey','grey','grey','grey','grey','grey'),split.by = 'condition')
DimPlot(kera)
```

## Figure 3C & D: plotting genes (ELOB, PHB, CITED4, GAPDH, COX7A1 and AKR1C1)
```
VlnPlot(kera,'ELOB')
VlnPlot(kera,'PHB')
VlnPlot(kera,'CITED4')
VlnPlot(kera,'GAPDH')
VlnPlot(kera,'COX7A1')
VlnPlot(kera,'AKR1C1')
```


## Figure 3E & F: for heatmaps, scale genes values and define lists of genes of interest
```
metabolism_genes = c('COX8A','COX7B','SURF1','NDUFA4','COX4I1','COX6A1','COX5B','COX6C','COX7C','COX5A','COX6B1','GLS','PRDX2','PRDX5','COX7A2L','PRDX1','CYCS','COX14','LAMTOR1','COX20','LAMTOR5')
glycolysis_genes = c('PKM','TPI1','PPA1','PGAM1','PGK1','ENO1','ALDOA','GAPDH')
kera = ScaleData(kera, features = rownames(kera))
DoHeatmap(kera, metabolism_genes)
DoHeatmap(kera, glycolysis_genes)
```

## Figure 4A
```
VlnPlot(object,'KRT14')
```

## Figure 4B
```
DimPlot(object, cols = c('grey','grey','grey','grey','grey','purple','grey','grey','grey','grey','grey'))
cluster5 = subset(object,subset = seurat_clusters==5)
cluster5 = RunTSNE(cluster5)
FeaturePlot(cluster5,cols = c('grey','blue','black'), 'KRT14')
```

## Figure 4C & D
```
FeaturePlot(cluster5,cols = c('grey','blue','black'), split.by = 'condition', 'TP53')
FeaturePlot(cluster5,cols = c('grey','blue','black'), split.by = 'condition', 'ADAM17')
FeaturePlot(cluster5,cols = c('grey','blue','black'), split.by = 'condition', 'NOTCH1')
FeaturePlot(cluster5,cols = c('grey','blue','black'), split.by = 'condition', 'SMURF1')
FeaturePlot(cluster5,cols = c('grey','blue','black'), split.by = 'condition', 'TWIST1')
```

## Supplementary Figure 3A: dotplot
```
c('MLANA','PMEL','DCT','JCHAIN','IGHM','IGHG1','LYVE1','TFF3','CCL21','TPSB2','TPSAB1','CTSG','KRT7','KRT19','AQP5','KRT1','KRT5','KRT14','IL1B','LYZ','CXCL8','TPM2','ACTA2','TAGLN','CXCR4','GNLY','NKG7','CDH5','VWF','SELE','COL1A2','DCN','COL1A1')
DotPlot(object, features = markers)+ theme(axis.text.x = element_text(angle = 45, hjust=1))
```
## Supplementary Figure 3B: signature genes
```
FeaturePlot(object, cols = c('grey','blue','black'), 'COL1A1')
FeaturePlot(object, cols = c('grey','blue','black'), 'CDH5')
FeaturePlot(object, cols = c('grey','blue','black'), 'NKG7')
FeaturePlot(object, cols = c('grey','blue','black'), 'ACTA2')
FeaturePlot(object, cols = c('grey','blue','black'), 'LYZ')
FeaturePlot(object, cols = c('grey','blue','black'), 'KRT14')
FeaturePlot(object, cols = c('grey','blue','black'), 'KRT7')
FeaturePlot(object, cols = c('grey','blue','black'), 'TPSAB1')
FeaturePlot(object, cols = c('grey','blue','black'), 'TFF3')
FeaturePlot(object, cols = c('grey','blue','black'), 'JCHAIN')
FeaturePlot(object, cols = c('grey','blue','black'), 'MLANA')
```
## Supplementary Figure 3C: Extract cell numbers and percentage from each cluster
```
cell_numbers = table(object@meta.data$condition,object@meta.data$seurat_clusters)
row_sums = rowSums(cell_numbers)
cell_percent = cell_numbers/row_sums*100
write.csv(cell_number,'/cell_number.csv')
write.csv(cell_percent,'/cell_percent.csv')
```

## Supplementary Figure 5B: Violin plots
```
VlnPlot(object, split.by = 'condition', 'TP53')
VlnPlot(object, split.by = 'condition', 'ADAM17')
VlnPlot(object, split.by = 'condition', 'NOTCH1')
VlnPlot(object, split.by = 'condition', 'SMURF1')
VlnPlot(object, split.by = 'condition', 'TWIST1')
```
# If proceeding with the matrix file:
### Files can be found under accession number GSM5364333, GSM5364334, GSM5364335 and GSM5364336 for normal skin samples and GSM5364337, GSM5364338 and GSM5364339 for skin with chronic wound

## Crate a variable for each sample folder including the 3 files (barcodes, features and the matrix)
```
normal1 = "/GSM5364333"
normal2 = "/GSM5364334"
normal3 = "/GSM5364335"
normal4 = "/GSM5364336"
wound1 = "/GSM5364337"
wound2 = "/GSM5364338"
wound3 = "/GSM5364339"
```

## Read the data from the folders in R
```
normal1 = Read10X(normal1)
normal2 = Read10X(normal2)
normal3 = Read10X(normal3)
normal4 = Read10X(normal4)
wound1 = Read10X(wound1)
wound2 = Read10X(wound2)
wound3 = Read10X(wound3)
```

## Create a Seurat object for each sample
#### Genes expressed (>0) in less than 3 cells were filtered out
#### Create a Seurat object for each sample

```
normal1 = CreateSeuratObject(counts = normal1, min.cells = 3, min.features = 200, project = "NS2" )
normal2 = CreateSeuratObject(counts = normal2, min.cells = 3, min.features = 200, project = "NS8" )
normal3 = CreateSeuratObject(counts = normal3, min.cells = 3, min.features = 200, project = "NS9" )
normal4 = CreateSeuratObject(counts = normal4, min.cells = 3, min.features = 200, project = "NS10" )
wound1 = CreateSeuratObject(counts = wound1, min.cells = 3, min.features = 200, project = "HU2" )
wound2 = CreateSeuratObject(counts = wound2, min.cells = 3, min.features = 200, project = "HU4" )
wound3 = CreateSeuratObject(counts = wound3, min.cells = 3, min.features = 200, project = "HU9" )
```

## Calculate the percentage of mitochondrial gene expression for each cell within each sample
```
normal1[["percent.mt"]] <- PercentageFeatureSet(normal1, pattern = "^MT-")
normal2[["percent.mt"]] <- PercentageFeatureSet(normal2, pattern = "^MT-")
normal3[["percent.mt"]] <- PercentageFeatureSet(normal3, pattern = "^MT-")
normal4[["percent.mt"]] <- PercentageFeatureSet(normal4, pattern = "^MT-")
wound1[["percent.mt"]] <- PercentageFeatureSet(wound1, pattern = "^MT-")
wound2[["percent.mt"]] <- PercentageFeatureSet(wound2, pattern = "^MT-")
wound3[["percent.mt"]] <- PercentageFeatureSet(wound3, pattern = "^MT-")
```

## Perform log normalization with 10,000 as a scaling factor
```
normal1 <- NormalizeData(normal1, normalization.method = "LogNormalize", scale.factor = 10000)
normal2 <- NormalizeData(normal2, normalization.method = "LogNormalize", scale.factor = 10000)
normal3 <- NormalizeData(normal3, normalization.method = "LogNormalize", scale.factor = 10000)
normal4 <- NormalizeData(normal4, normalization.method = "LogNormalize", scale.factor = 10000)
wound1 <- NormalizeData(wound1, normalization.method = "LogNormalize", scale.factor = 10000)
wound2 <- NormalizeData(wound2, normalization.method = "LogNormalize", scale.factor = 10000)
wound3 <- NormalizeData(wound3, normalization.method = "LogNormalize", scale.factor = 10000)
```

## Identify the top 2,000 variable genes within each sample
```
normal1 <- FindVariableFeatures(normal1, selection.method = "vst", nfeatures = 2000)
normal2 <- FindVariableFeatures(normal2, selection.method = "vst", nfeatures = 2000)
normal3 <- FindVariableFeatures(normal3, selection.method = "vst", nfeatures = 2000)
normal4 <- FindVariableFeatures(normal4, selection.method = "vst", nfeatures = 2000)
wound1 <- FindVariableFeatures(wound1, selection.method = "vst", nfeatures = 2000)
wound2 <- FindVariableFeatures(wound2, selection.method = "vst", nfeatures = 2000)
wound3 <- FindVariableFeatures(wound3, selection.method = "vst", nfeatures = 2000)
```

## Filter cells
```
normal1 <- subset(normal1, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15 & nCount_RNA < 40000 & nCount_RNA > 1000)
normal2 <- subset(normal2, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15 & nCount_RNA < 40000 & nCount_RNA > 1000)
normal3 <- subset(normal3, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15 & nCount_RNA < 40000 & nCount_RNA > 1000)
normal4 <- subset(normal4, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15 & nCount_RNA < 40000 & nCount_RNA > 1000)
wound1 <- subset(wound1, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15 & nCount_RNA < 40000 & nCount_RNA > 1000)
wound2 <- subset(wound2, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15 & nCount_RNA < 40000 & nCount_RNA > 1000)
wound3 <- subset(wound3, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15 & nCount_RNA < 40000 & nCount_RNA > 1000)
```

## Integrating individual samples into one Seurat object
```
anchors <- FindIntegrationAnchors(object.list = list(normal1, normal2, normal3, normal4, wound1, wound2, wound3), dims = 1:20)
object <- IntegrateData(anchorset = anchors, dims = 1:20)
```

## Run principal component analysis
```
object <- ScaleData(object, verbose = TRUE)
object <- RunPCA(object = object, npcs = 50, verbose = FALSE)
```

## Run tSNE
```
object = RunTSNE(object, dims = 1:21)
DimPlot(object)
```

## Assign conditions to samples (Normal skin and chronic wound)
```
object$condition = 1
temp = object@meta.data

x = dim(temp)[1]
for (i in 1:x)
{
    if(temp[i,1] == 'wound1')
    {
        temp[i,5] = "Chronic wound"
    }
    else if(temp[i,1] == 'wound2')
    {
        temp[i,5] = "Chronic wound"
    }
    else if(temp[i,1] == 'wound3')
    {
        temp[i,5] = "Chronic wound"
    }
    else
    {
        temp[i,5] = "Normal skin"
    }
}
object@meta.data = temp
Idents(object) = "condition"
```

## Clustering
```
object <- FindNeighbors(object = object, dims = 1:21)
object <- FindClusters(object, resolution = 0.15)
```
