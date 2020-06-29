Cory Williams
Programmer Script
R-verion 3.6.0

#loading necessary downstream packages
library(dplyr)
library(Seurat)
library(patchwork)


#all packages and scripts were obtained from the publishing website and ran in an interactive R-studio envrironment 
provided by BU clustering computing network

#assigns the location of the alevin output file to a variable to be called on
files <- "/projectnb/bf528/users/group_8/project_4/programmer/GSM2230760__salmon_quant/alevin/quants_mat.gz"


#loads the package needed to import alevin data
library("tximport", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.6")

#uses txiport to import alevin file ad txi variable
txi <- tximport(files, type="alevin")

#crearting Seurat object using txi$counts which is the cell count and assinging it to vairable "pbmc"
pbmc <- CreateSeuratObject(counts = txi$counts , min.cells = 3, min.features = 200, project = "10X_PBMC")

#plot initial numbers of cell count and gene expression count per cell pre filtering and noramlization
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3

#filtering out genes higher than 200 but lower than 2500 and applying it to
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#account for variation
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#set a standard for gene variation
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


#cluster cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)