# Final analysis script of the COVID Tcell scRNAseq project.
# This script produces the scRNAseq plots of the paper.


# Loading packages
library(openxlsx)
library(presto)
library(msigdbr)
library(fgsea)
library(Seurat)
library(scRepertoire)
library(umap)
library(tidyverse)
library(SeuratData)
library(cowplot)
library(dplyr)
library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(VennDiagram)
library(patchwork)
library(ggpubr)
library(pals)
library(ggvenn)
library(eulerr)
options(ggrepel.max.overlaps = Inf)

# Clearing the environment
rm(list = ls())

###########################################################################################################################################
# Part 1 - Loading and preparation of the main data
###########################################################################################################################################
# Setting the working directory
setwd("~/NAS/Jan/Experiments and Data/210426 - Exp 033 - CD8 Tcell scRNAseq (Boyman)/cellranger_multi_outputs")

# Loading the data into R
Set1_sorted <- "./multi_Set_1/outs/count/filtered_feature_bc_matrix"
Set2_sorted <- "./multi_Set_2_Sorted/outs/count/filtered_feature_bc_matrix"
Set3_sorted <- "./multi_Set_3/outs/count/filtered_feature_bc_matrix"
Set4_sorted <- "./multi_Set_4_Sorted/outs/count/filtered_feature_bc_matrix"

Set1_sorted.data <- Read10X(data.dir = Set1_sorted)
Set2_sorted.data <- Read10X(data.dir = Set2_sorted)
Set3_sorted.data <- Read10X(data.dir = Set3_sorted)
Set4_sorted.data <- Read10X(data.dir = Set4_sorted)


# Removing not needed file links
remove(Set1_sorted,Set2_sorted,Set3_sorted,Set4_sorted)


# Seurat Objects are created using the Gene expression data.
Set1_sorted <- CreateSeuratObject(counts = Set1_sorted.data$`Gene Expression`)
Set2_sorted <- CreateSeuratObject(counts = Set2_sorted.data$`Gene Expression`)
Set3_sorted <- CreateSeuratObject(counts = Set3_sorted.data$`Gene Expression`)
Set4_sorted <- CreateSeuratObject(counts = Set4_sorted.data$`Gene Expression`)


# One assay is created for the surface protein quantification, one assay is created for the dextramers
# Protein capture assays (Only for the sorted datasets):
Set1_sorted_Protein_assay <- CreateAssayObject(counts = Set1_sorted.data$`Antibody Capture`[c(8:9),])
Set2_sorted_Protein_assay <- CreateAssayObject(counts = Set2_sorted.data$`Antibody Capture`[c(8:9),])
Set3_sorted_Protein_assay <- CreateAssayObject(counts = Set3_sorted.data$`Antibody Capture`[c(8:9),])
Set4_sorted_Protein_assay <- CreateAssayObject(counts = Set4_sorted.data$`Antibody Capture`[c(8:9),])

# Dextramer capture assays (Only for the sorted datasets):
Set1_sorted_Dextramer_assay <- CreateAssayObject(counts = Set1_sorted.data$`Antibody Capture`[c(1:7),])
Set2_sorted_Dextramer_assay <- CreateAssayObject(counts = Set2_sorted.data$`Antibody Capture`[c(1:7),])
Set3_sorted_Dextramer_assay <- CreateAssayObject(counts = Set3_sorted.data$`Antibody Capture`[c(1:7),])
Set4_sorted_Dextramer_assay <- CreateAssayObject(counts = Set4_sorted.data$`Antibody Capture`[c(1:7),])

# Now the assays are added to the previously created Seurat objects
Set1_sorted[["Protein"]] <- Set1_sorted_Protein_assay
Set2_sorted[["Protein"]] <- Set2_sorted_Protein_assay
Set3_sorted[["Protein"]] <- Set3_sorted_Protein_assay
Set4_sorted[["Protein"]] <- Set4_sorted_Protein_assay

Set1_sorted[["Dextramer"]] <- Set1_sorted_Dextramer_assay
Set2_sorted[["Dextramer"]] <- Set2_sorted_Dextramer_assay
Set3_sorted[["Dextramer"]] <- Set3_sorted_Dextramer_assay
Set4_sorted[["Dextramer"]] <- Set4_sorted_Dextramer_assay

remove(Set1_sorted_Dextramer_assay,Set1_sorted_Protein_assay,Set2_sorted_Dextramer_assay,
       Set2_sorted_Protein_assay,Set3_sorted_Dextramer_assay,Set3_sorted_Protein_assay,Set4_sorted_Dextramer_assay,Set4_sorted_Protein_assay)

remove(Set1_sorted.data,Set2_sorted.data,Set3_sorted.data,Set4_sorted.data)

# QC and selecting cells for further analysis
Set1_sorted[["percent.mt"]] <- PercentageFeatureSet(Set1_sorted, pattern = "^MT-")
Set2_sorted[["percent.mt"]] <- PercentageFeatureSet(Set2_sorted, pattern = "^MT-")
Set3_sorted[["percent.mt"]] <- PercentageFeatureSet(Set3_sorted, pattern = "^MT-")
Set4_sorted[["percent.mt"]] <- PercentageFeatureSet(Set4_sorted, pattern = "^MT-")

Set1_sorted <- subset(Set1_sorted, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
Set2_sorted <- subset(Set2_sorted, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
Set3_sorted <- subset(Set3_sorted, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
Set4_sorted <- subset(Set4_sorted, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# How many cells are left in each of the datasets?
nrow(Set1_sorted@meta.data)
nrow(Set2_sorted@meta.data)
nrow(Set3_sorted@meta.data)
nrow(Set4_sorted@meta.data)

# Now, normalization of the Gene expression data and of the assay data is done.
Set1_sorted <- NormalizeData(Set1_sorted)
Set2_sorted <- NormalizeData(Set2_sorted)
Set3_sorted <- NormalizeData(Set3_sorted)
Set4_sorted <- NormalizeData(Set4_sorted)

Set1_sorted <- NormalizeData(Set1_sorted, assay= "Protein", normalization.method = "CLR", margin = 2)
Set2_sorted <- NormalizeData(Set2_sorted, assay= "Protein", normalization.method = "CLR", margin = 2)
Set3_sorted <- NormalizeData(Set3_sorted, assay= "Protein", normalization.method = "CLR", margin = 2)
Set4_sorted <- NormalizeData(Set4_sorted, assay= "Protein", normalization.method = "CLR", margin = 2)


# Adding SNP cluster identities to cells using the souporcell output
# 1) Load data into 
SNPs.1_2 <- read.csv("snp_demux_set_1_2_merged/clusters.tsv", sep = "\t")
SNPs.1_2 <- SNPs.1_2[,c(1:3)]
table(SNPs.1_2$status)
SNPs.1_2 <- SNPs.1_2[SNPs.1_2$status=="singlet",]

SNPs.3_4 <- read.csv("snp_demux_set_3_4_merged/clusters.tsv", sep = "\t")
SNPs.3_4 <- SNPs.3_4[,c(1:3)]
table(SNPs.3_4$status)
SNPs.3_4 <- SNPs.3_4[SNPs.3_4$status=="singlet",]

# Making vectors containing rownames (part of preparation for the merging by rowname, which follows after)
rownames.Set1_sorted <- rownames(Set1_sorted@meta.data)
rownames.Set2_sorted <- rownames(Set2_sorted@meta.data)
rownames.Set3_sorted <- rownames(Set3_sorted@meta.data)
rownames.Set4_sorted <- rownames(Set4_sorted@meta.data)

rownames(Set1_sorted@meta.data) <- paste0(rownames(Set1_sorted@meta.data),"_Set_1")
rownames(Set2_sorted@meta.data) <- paste0(rownames(Set2_sorted@meta.data),"_Set_2_Sorted")
rownames(Set3_sorted@meta.data) <- paste0(rownames(Set3_sorted@meta.data),"_Set_3")
rownames(Set4_sorted@meta.data) <- paste0(rownames(Set4_sorted@meta.data),"_Set_4_Sorted")

rownames(SNPs.1_2) <- SNPs.1_2$barcode
rownames(SNPs.3_4) <- SNPs.3_4$barcode

# Do the merge and replace NA values with 'unknown'
Set1_sorted@meta.data <- merge(Set1_sorted@meta.data,SNPs.1_2,by=0,all.x = TRUE)
Set2_sorted@meta.data <- merge(Set2_sorted@meta.data,SNPs.1_2,by=0,all.x = TRUE)
Set3_sorted@meta.data <- merge(Set3_sorted@meta.data,SNPs.3_4,by=0,all.x = TRUE)
Set4_sorted@meta.data <- merge(Set4_sorted@meta.data,SNPs.3_4,by=0,all.x = TRUE)

nrow(Set1_sorted@meta.data)
table(Set1_sorted@meta.data$status,useNA = "always")
table(Set2_sorted@meta.data$status,useNA = "always")
table(Set3_sorted@meta.data$status,useNA = "always")
table(Set4_sorted@meta.data$status,useNA = "always")

colnames(Set1_sorted@meta.data)
table(Set1_sorted@meta.data$status, useNA = "always")

Set1_sorted@meta.data$status[is.na(Set1_sorted@meta.data$status)] <- "unknown"
Set2_sorted@meta.data$status[is.na(Set2_sorted@meta.data$status)] <- "unknown"
Set3_sorted@meta.data$status[is.na(Set3_sorted@meta.data$status)] <- "unknown"
Set4_sorted@meta.data$status[is.na(Set4_sorted@meta.data$status)] <- "unknown"

Set1_sorted@meta.data$assignment[is.na(Set1_sorted@meta.data$assignment)] <- "unknown"
Set2_sorted@meta.data$assignment[is.na(Set2_sorted@meta.data$assignment)] <- "unknown"
Set3_sorted@meta.data$assignment[is.na(Set3_sorted@meta.data$assignment)] <- "unknown"
Set4_sorted@meta.data$assignment[is.na(Set4_sorted@meta.data$assignment)] <- "unknown"

rownames(Set1_sorted@meta.data) <- rownames.Set1_sorted 
rownames(Set2_sorted@meta.data) <- rownames.Set2_sorted 
rownames(Set3_sorted@meta.data) <- rownames.Set3_sorted 
rownames(Set4_sorted@meta.data) <- rownames.Set4_sorted 

remove(rownames.Set1_sorted, rownames.Set2_sorted, rownames.Set3_sorted, rownames.Set4_sorted)
remove(SNPs.1_2,SNPs.3_4)


# Adding suffix to SNP assignment to keep info from which Souporcell run it comes
Set1_sorted@meta.data$assignment <- paste0(Set1_sorted@meta.data$assignment,"_1")
Set2_sorted@meta.data$assignment <- paste0(Set2_sorted@meta.data$assignment,"_1")
Set3_sorted@meta.data$assignment <- paste0(Set3_sorted@meta.data$assignment,"_2")
Set4_sorted@meta.data$assignment <- paste0(Set4_sorted@meta.data$assignment,"_2")


# Addition of the TCR data
Set1_sorted_TCR <- read.csv("./multi_Set_1/outs/vdj_t/filtered_contig_annotations.csv")
Set2_sorted_TCR <- read.csv("./multi_Set_2_Sorted/outs/vdj_t/filtered_contig_annotations.csv")
Set3_sorted_TCR <- read.csv("./multi_Set_3/outs/vdj_t/filtered_contig_annotations.csv")
Set4_sorted_TCR <- read.csv("./multi_Set_4_Sorted/outs/vdj_t/filtered_contig_annotations.csv")

TCRcombined1 <- combineTCR(Set1_sorted_TCR,samples = "A1",ID="A1",cells ="T-AB", filterMulti = T)
TCRcombined2 <- combineTCR(Set2_sorted_TCR,samples = "A1",ID="A1",cells ="T-AB", filterMulti = T)
TCRcombined3 <- combineTCR(Set3_sorted_TCR,samples = "A1",ID="A1",cells ="T-AB", filterMulti = T)
TCRcombined4 <- combineTCR(Set4_sorted_TCR,samples = "A1",ID="A1",cells ="T-AB", filterMulti = T)

TCRcombined1$A1_A1$barcode <- substring(TCRcombined1$A1_A1$barcode, 7)
TCRcombined2$A1_A1$barcode <- substring(TCRcombined2$A1_A1$barcode, 7)
TCRcombined3$A1_A1$barcode <- substring(TCRcombined3$A1_A1$barcode, 7)
TCRcombined4$A1_A1$barcode <- substring(TCRcombined4$A1_A1$barcode, 7)

Set1_sorted <- combineExpression(TCRcombined1, Set1_sorted, cloneCall="aa", proportion = F, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
Set2_sorted <- combineExpression(TCRcombined2, Set2_sorted, cloneCall="aa", proportion = F, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
Set3_sorted <- combineExpression(TCRcombined3, Set3_sorted, cloneCall="aa", proportion = F, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
Set4_sorted <- combineExpression(TCRcombined4, Set4_sorted, cloneCall="aa", proportion = F, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

remove(Set1_sorted_TCR,Set2_sorted_TCR,Set3_sorted_TCR,Set4_sorted_TCR)
remove(TCRcombined1,TCRcombined2,TCRcombined3,TCRcombined4)

# Are there same clones based on "CTaa" in both timepoints? - These will be defined as common clones.
common_Set1.2 <- intersect(Set1_sorted@meta.data$CTaa, Set2_sorted@meta.data$CTaa)
common_Set3.4 <- intersect(Set3_sorted@meta.data$CTaa, Set4_sorted@meta.data$CTaa)
common_Set1.2 <- common_Set1.2[!is.na(common_Set1.2)]
common_Set3.4 <- common_Set3.4[!is.na(common_Set3.4)]

# Are there overlaps of clones between datasets?
common_Set1.3 <- intersect(Set1_sorted@meta.data$CTaa, Set3_sorted@meta.data$CTaa)
common_Set1.4 <- intersect(Set1_sorted@meta.data$CTaa, Set4_sorted@meta.data$CTaa)
common_Set2.3 <- intersect(Set2_sorted@meta.data$CTaa, Set3_sorted@meta.data$CTaa)
common_Set2.4 <- intersect(Set2_sorted@meta.data$CTaa, Set4_sorted@meta.data$CTaa)

common_Set1.3 <- common_Set1.3[!is.na(common_Set1.3)]
common_Set1.4 <- common_Set1.4[!is.na(common_Set1.4)]
common_Set2.3 <- common_Set2.3[!is.na(common_Set2.3)]
common_Set2.4 <- common_Set2.4[!is.na(common_Set2.4)]

Overlaps <- c(common_Set1.3,common_Set1.4,common_Set2.3)
Overlaps <- unique(Overlaps)
remove(Overlaps)
remove(common_Set1.3,common_Set1.4,common_Set2.3,common_Set2.4)

# The data.frames that are created here, include from each set the cells which have assigned clonotypes that are found 
# in the respective first and second timepoint.
Set1_commons <- filter(Set1_sorted@meta.data, grepl(paste(common_Set1.2, collapse="|"), CTaa))
Set2_commons <- filter(Set2_sorted@meta.data, grepl(paste(common_Set1.2, collapse="|"), CTaa))
Set3_commons <- filter(Set3_sorted@meta.data, grepl(paste(common_Set3.4, collapse="|"), CTaa))
Set4_commons <- filter(Set4_sorted@meta.data, grepl(paste(common_Set3.4, collapse="|"), CTaa))

remove(common_Set1.2,common_Set3.4)

# Adding a colum to the metadata that specifies whether a cell is among the once that have assigned clonotypes that are 
# found in the respective first and second timepoint (- basically adding a "yes" to all cells that appear in
# Set1_commons 1-4).

Set1_sorted@meta.data$common.clone <- "no"
Set2_sorted@meta.data$common.clone <- "no"
Set3_sorted@meta.data$common.clone <- "no"
Set4_sorted@meta.data$common.clone <- "no"


Set1_sorted@meta.data[rownames(Set1_commons),19] <- "yes"
Set2_sorted@meta.data[rownames(Set2_commons),19] <- "yes"
Set3_sorted@meta.data[rownames(Set3_commons),19] <- "yes"
Set4_sorted@meta.data[rownames(Set4_commons),19] <- "yes"

remove(Set1_commons,Set2_commons,Set3_commons,Set4_commons)

#Perform FindVariableFeatures
Set1_sorted <- FindVariableFeatures(Set1_sorted, selection.method = "vst", nfeatures = 2000)
Set2_sorted <- FindVariableFeatures(Set2_sorted, selection.method = "vst", nfeatures = 2000)
Set3_sorted <- FindVariableFeatures(Set3_sorted, selection.method = "vst", nfeatures = 2000)
Set4_sorted <- FindVariableFeatures(Set4_sorted, selection.method = "vst", nfeatures = 2000)

# Rename cells to keep dataset info
Set1_sorted <- RenameCells(object = Set1_sorted, add.cell.id = "Set1")
Set2_sorted <- RenameCells(object = Set2_sorted, add.cell.id = "Set2")
Set3_sorted <- RenameCells(object = Set3_sorted, add.cell.id = "Set3")
Set4_sorted <- RenameCells(object = Set4_sorted, add.cell.id = "Set4")

# Add columns to keep dataset info
Set1_sorted$Dataset <- "Set1"
Set2_sorted$Dataset <- "Set2"
Set3_sorted$Dataset <- "Set3"
Set4_sorted$Dataset <- "Set4"

Set1_sorted$Timepoint <- "Acute"
Set2_sorted$Timepoint <- "Recovered"
Set3_sorted$Timepoint <- "Acute"
Set4_sorted$Timepoint <- "Recovered"




###########################################################################################################################################
# Part 2 - Loading and preparation of the pilot data (take healthy patient cells only)
###########################################################################################################################################
# Creating links to the data.
data_dir <- "./Pilot data/filtered_feature_bc_matrix"
SNPs_clusters <- read.csv("./Pilot data/souporcell_5samples.TCR_Boyman_mar2021/snp_demux_Tcells_Boyman_allCustom/clusters.tsv", sep = "\t")

# Loading the data into R
Tcells <- Read10X(data.dir = data_dir)
remove(data_dir)

# Here the Seurat object is created and the counts are normalized
Tcell <- CreateSeuratObject(counts = Tcells$`Gene Expression`)
Tcell <- NormalizeData(Tcell)

# Here, the HTO data and the TotalSeq data and the dCODE data are added as independent assays.
Tcell[["HTO"]] <- CreateAssayObject(counts = Tcells$Custom[c(1:10),])
Tcell[["dCODEs"]] <- CreateAssayObject(counts = Tcells$Custom[c(12:14),])
Tcell[["SurfaceProtein"]] <- CreateAssayObject(counts=as.matrix(Tcells$Custom)[11,, drop = FALSE])

# Also, the these assay objects are normalized.
# Here, it is not clear to me yet, for which assays this is necessary.
Tcell <- NormalizeData(Tcell, assay = "HTO", normalization.method = "CLR")
Tcell <- NormalizeData(Tcell, assay = "SurfaceProtein", normalization.method = "CLR")

# Data filtering
Tcell[["percent.mt"]] <- PercentageFeatureSet(Tcell, pattern = "^MT-")

# This is the step where the filtering by mitochondrial genes and nFeature counts happens.
# For this dataset, I have decided to go with a filtering step that exludes cells with >10% mitochondrial genes.
Tcell <- subset(Tcell, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# In this step, the cells are demultiplexed based on HTO enrichment.
Tcell <- HTODemux(Tcell, assay = "HTO", positive.quantile = 0.99)
table(Tcell@meta.data$hash.ID)

# Now the vdj data is added.
# First, the clonotypedata output file from cellranger is loaded into R
vdj <- read.csv("./Pilot data/filtered_contig_annotations.csv")

# As the output of CellRanger are quantifications of both chains, 
# the next step is to create a single list object with the by cell barcode.
vdjcombined1 <- combineTCR(vdj,samples = "A1",ID="A1")

# Then, the clonotypedata is paired with the seurat object that we already have.
# Importantly, the major requirement for the attachment is matching contig cell 
# barcodes and barcodes in the row names of the meta data of the seurat object. 
# If these do not match, the attachment will fail.
vdjcombined1$A1_A1$barcode <- substring(vdjcombined1$A1_A1$barcode, 7)
Tcell <- combineExpression(vdjcombined1, Tcell, cloneCall="CTaa")

# Now I combine the HTO demultiplexed data with the SNP demultiplexed data.
shortSNPs_clusters <- SNPs_clusters [,c(1:3)]
Tcell@meta.data$barcode <- rownames(Tcell@meta.data)
Tcell@meta.data <- merge(Tcell@meta.data, shortSNPs_clusters, by = 'barcode')
rownames(Tcell@meta.data) <- Tcell@meta.data$barcode
Tcell@meta.data

df <- as.data.frame(count(Tcell@meta.data, HTO_classification.global))
plot1 <- ggplot(df, aes(HTO_classification.global,n))+geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90))+
  geom_text(aes(label = n), size = 3, hjust = 0.5, vjust = -1, position ="stack")
plot1
df2 <- as.data.frame(count(Tcell@meta.data, status))
plot2 <- ggplot(df2, aes(status,n))+geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90))+
  geom_text(aes(label = n), size = 3, hjust = 0.5, vjust = -1, position ="stack")
plot2
plot1 + plot2

Tcell2 <- subset(Tcell, subset = HTO_classification.global=='Singlet')
Tcell2 <- subset(Tcell2, subset = status=='singlet')
nrow(Tcell2@meta.data)
head(Tcell2@meta.data)
tbl <- with(Tcell2@meta.data, table(assignment, HTO_maxID))
tbl
barplot(tbl, beside = TRUE, legend = TRUE)

# The plot above shows that the SNPs and the HTO data resemble each other very well. I am going to keep all healthy donor cells,
# which are SNPs clusters 1,2,3,4
keep.clusters <- c(1,2,3,4)
Tcell <- subset(Tcell,  subset = assignment %in% keep.clusters)
Tcell@meta.data$assignment <- "Healthy"

# FindVariableFeatures is performed.
Tcell <- FindVariableFeatures(Tcell, selection.method = "mean.var.plot")

#Preparing for integration
Tcell <- RenameCells(object = Tcell, add.cell.id = "Healthy")
Tcell$Dataset <- "Healthy"
Tcell$Timepoint <- "Healthy"

# Removal of unnessecary files
remove(df, df2, plot1, plot2, shortSNPs_clusters,SNPs_clusters,Tcell2,Tcells,vdj,vdjcombined1,keep.clusters,tbl)




###########################################################################################################################################
# Part 3 - Integration of the datasets to generate one dataset that contains all cells
###########################################################################################################################################
# I am following the instructions for integration from:
# https://satijalab.org/seurat/articles/integration_introduction.html
# Create lists containing the datasets
Dataset_list <- list(Set1_sorted,Set2_sorted,Set3_sorted,Set4_sorted,Tcell)


# Select features that are repeatedly variable across datasets for integration
Integreation_features <- SelectIntegrationFeatures(object.list = Dataset_list)


# Perform integration --> Takes long!
Integration_anchors <- FindIntegrationAnchors(object.list = Dataset_list, anchor.features = Integreation_features)


# This command creates an 'integrated' data assay
Integrated <- IntegrateData(anchorset = Integration_anchors)

remove(Dataset_list)
remove(Integration_anchors)
remove(Integreation_features)
remove(Set1_sorted,Set2_sorted,Set3_sorted,Set4_sorted)
remove(Tcell)

Integrated@meta.data$assignment[is.na(Integrated@meta.data$Patient)] <- "Healthy"
Integrated@meta.data$status[is.na(Integrated@meta.data$Patient)] <- "Healthy"

###########################################################################################################################################
# Part 4 - Dimensional reduction and clustering of the Integrated dataset
###########################################################################################################################################
Integrated <- ScaleData(Integrated)
Integrated <- RunPCA(Integrated, features = VariableFeatures(object = Integrated))
Integrated <- FindNeighbors(Integrated, dims = 1:15,reduction = "pca")
Integrated <- FindClusters(Integrated, resolution = 0.5)
Integrated <- RunUMAP(Integrated, dims = 1:15)

# We enable here the identification of heterogeneously expressed genes across clusters. Commented out for running speed reasons.
#markers <- FindAllMarkers(Integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,assay = "RNA")
#markers.table <- markers %>%group_by(cluster) %>%top_n(n = 25, wt = avg_log2FC)
#markers.table <- markers.table[,c(2,5,6,7)]
#write.xlsx(markers.table, 'All_cluster_markers.xlsx')


###########################################################################################################################################
# Part 5 - General analysis comparing cells from Healthy, Acute and Recovered patients and Severity conditions
###########################################################################################################################################
Integrated$Timepoint <- factor(x = Integrated$Timepoint, levels = c("Healthy", "Acute","Recovered"))
coul2 <- c("dodgerblue2", "#E31A1C","green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2","maroon", "deeppink1", "blue1", "steelblue4","green1","yellow3")
DimPlot(Integrated, reduction = "umap",label = F,pt.size = 0.1, cols = coul2)+labs(x = "UMAP 1", y="UMAP 2")+xlim(-10, 15)+ylim(-10, 15)

# We noticed that clusters 12 and 13 are very small and not interesting for the analysis, therefore we exclude them at this point
Integrated <- subset(Integrated, idents = c(0:11))
Integrated@meta.data$seurat_clusters <- droplevels(Integrated@meta.data$seurat_clusters)
DimPlot(Integrated, reduction = "umap",label = F,pt.size = 0.1, cols = coul2)+labs(x = "UMAP 1", y="UMAP 2")+xlim(-10, 15)+ylim(-10, 15)

# Creating a dot plot
features <- c("HLA-DRB5","HLA-DPB1","HLA-DQB1","EOMES","TIGIT","RRM2","TK1","CENPF","CENPM","MKI67","MCM4","CENPU",
              "IFIT1","IFIT2","IFIT3","IFIH1","TNF","LTA","IFNG","GZMB","GZMH","GNLY","PRF1","NKG7",
              "TCF7", "SELL", "IL7R", "CCR7")
DotPlot(Integrated, features = features,assay = "RNA",cols = c("lightgrey", "orangered1"))+ xlab("Marker genes")+ylab("Cluster")+coord_flip()
remove(features)

# Creating a UMAP where cells are colored by SNP cluster
# Prepare this also as percentage barplot
table(Integrated@meta.data$assignment)
Integrated_w.o_unkown <- subset(Integrated,subset = assignment !="unknown_1" & assignment !="unknown_2" & assignment !="Healthy")
table(Integrated_w.o_unkown@meta.data$assignment)
DimPlot(Integrated_w.o_unkown, reduction = "umap",label = F,pt.size = 0.1,group.by = "assignment")+labs(x = "UMAP 1", y="UMAP 2")+xlim(-10, 15)+ylim(-10, 15)+ggtitle("")

df <- Integrated_w.o_unkown@meta.data
df <- df %>% group_by(seurat_clusters,assignment) %>% add_count(seurat_clusters,name = "cells.percluster.and.assignment") %>% ungroup()
df <- df %>% group_by(assignment) %>% add_count(assignment,name = "cells.per.assignment") %>% ungroup()
df$percent.cells.per.cluster <- round(100*(df$cells.percluster.and.assignment/df$cells.per.assignment),1)
colnames(df)
df <- df[,c(12,35,38)]
df <- unique(df)

# Every row in this bar plot is one of 20 SNP clusters (the two timepoints are not separated).
ggplot(df, aes(fill=factor(seurat_clusters, levels=c(0,1,2,3,4,5,6,7,8,9,10,11)), y=assignment, x=percent.cells.per.cluster)) + 
  geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ theme(text = element_text(size = 15))+
  labs(x = "", y="SNP assignment")+scale_x_continuous(labels=scales::percent)+labs(fill = "Cluster")+
  scale_fill_manual(values=c(coul2))+ theme(text = element_text(size = 17))

remove(df,Integrated_w.o_unkown)


# Now I am creating individual UMAPs for the three timepoints, only one color per plot
Healthy <- subset(Integrated, subset = Timepoint=="Healthy")
Acute <- subset(Integrated, subset = Timepoint=="Acute")
Recovered <- subset(Integrated, subset = Timepoint=="Recovered")

DimPlot(Healthy, reduction = "umap",group.by = 'orig.ident',label = F,pt.size = 0.25, cols = 'gray40')+labs(x = "UMAP 1", y="UMAP 2")+
  ggtitle("Healthy")+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.position="none")+xlim(-10, 10)+ylim(-10, 15)

DimPlot(Acute, reduction = "umap",group.by = 'orig.ident',label = F,pt.size = 0.25, cols = 'palevioletred4')+labs(x = "UMAP 1", y="UMAP 2")+
  ggtitle("Acute")+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.position="none")+xlim(-10, 10)+ylim(-10, 15)

DimPlot(Recovered, reduction = "umap",group.by = 'orig.ident',label = F,pt.size = 0.25, cols = 'steelblue4')+labs(x = "UMAP 1", y="UMAP 2")+
  ggtitle("Recovered")+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.position="none")+xlim(-10, 10)+ylim(-10, 15)


# Show the cluster contributions per timepoint as a stacked percentage barplot
cluster.contributions.df <- as.data.frame(with(Integrated@meta.data, table(seurat_clusters, Timepoint)))
cluster.contributions.df <- cluster.contributions.df %>% group_by(Timepoint) %>% add_tally(Freq,name = "total.cells.pertimepoint") %>% ungroup()
cluster.contributions.df <- cluster.contributions.df %>% mutate(percentage.contribution= 100*(Freq/total.cells.pertimepoint))
cluster.contributions.df$percentage.contribution <- round(cluster.contributions.df$percentage.contribution,digits = 1)
cluster.contributions.df$percentage.contribution <- as.numeric(cluster.contributions.df$percentage.contribution)

cluster.contributions.df$Timepoint <- factor(cluster.contributions.df$Timepoint, levels = c("Recovered","Acute","Healthy"))

ggplot(cluster.contributions.df, aes(fill=factor(seurat_clusters, levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)), y=Timepoint, x=percentage.contribution)) + 
  geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ theme(text = element_text(size = 15))+
  labs(x = "Cluster contribution", y="Timepoints")+scale_x_continuous(labels=scales::percent)+labs(fill = "Cluster")+
  scale_fill_manual(values=c(coul2))+ theme(text = element_text(size = 17))

remove(cluster.contributions.df)
remove(Acute, Healthy, Recovered)


###########################################################################################################################################
# Part 6 - Dextramer part --> Define whether a cell is specific for a Dextramer
###########################################################################################################################################
# The first task is to remove the cells for which we don't have a TCR now.
Integrated_NA_filtered <- Integrated
Integrated_NA_filtered@meta.data$CTaa.known <- "yes"
Integrated_NA_filtered@meta.data$CTaa.known[is.na(Integrated_NA_filtered@meta.data$CTaa)] <- "no"
Integrated_NA_filtered <- subset(Integrated_NA_filtered, subset= CTaa.known == "yes")

Dextramer_counts <- as.data.frame(Integrated_NA_filtered@assays$Dextramer@counts)
Dextramer_counts_transposed <- data.table::transpose(Dextramer_counts)
rownames(Dextramer_counts_transposed) <- colnames(Dextramer_counts)
colnames(Dextramer_counts_transposed) <- rownames(Dextramer_counts)
Dextramer_counts <- Dextramer_counts_transposed
remove(Dextramer_counts_transposed)

# UMI count needs to be higher than 10 AND more than 5x higher than its negative control.
# Calculating ratios for each Dextramer vs. its negative control
Dextramer_counts$FTSDYYQLY.corrected <- Dextramer_counts$FTSDYYQLY/Dextramer_counts$`negative-Ctrl-A01-01`
Dextramer_counts$TTDPSFLGRY.corrected <- Dextramer_counts$TTDPSFLGRY/Dextramer_counts$`negative-Ctrl-A01-01`
Dextramer_counts$ATEGALNTPK.corrected <- Dextramer_counts$ATEGALNTPK/Dextramer_counts$`negative-Ctrl-Universal`
Dextramer_counts$KTFPPTEPK.corrected <- Dextramer_counts$KTFPPTEPK/Dextramer_counts$`negative-Ctrl-Universal`

# Classifying if a cell is Positive for a specific Dextramer
Dextramer_counts$FTSDYYQLY.classification <- "Negative"
Dextramer_counts$TTDPSFLGRY.classification <- "Negative"
Dextramer_counts$ATEGALNTPK.classification <- "Negative"
Dextramer_counts$KTFPPTEPK.classification <- "Negative"

Dextramer_counts$FTSDYYQLY.classification[Dextramer_counts$FTSDYYQLY>10 & (Dextramer_counts$FTSDYYQLY.corrected>5 | Dextramer_counts$FTSDYYQLY.corrected==Inf)] <- "Positive"
Dextramer_counts$TTDPSFLGRY.classification[Dextramer_counts$TTDPSFLGRY>10 & (Dextramer_counts$TTDPSFLGRY.corrected>5 | Dextramer_counts$TTDPSFLGRY.corrected==Inf)] <- "Positive"
Dextramer_counts$ATEGALNTPK.classification[Dextramer_counts$ATEGALNTPK>10 & (Dextramer_counts$ATEGALNTPK.corrected>5 | Dextramer_counts$ATEGALNTPK.corrected==Inf)] <- "Positive"
Dextramer_counts$KTFPPTEPK.classification[Dextramer_counts$KTFPPTEPK>10 & (Dextramer_counts$KTFPPTEPK.corrected>5 | Dextramer_counts$KTFPPTEPK.corrected==Inf)] <- "Positive"

# For how many Dextramers is a cell Positive?
Dextramer_counts$positive.for.n.dextramers <- rowSums(Dextramer_counts == "Positive") 

# Add the information whether a cell is Dextramer specific to the datasets
colnames(Integrated_NA_filtered@meta.data)
colnames(Dextramer_counts)
Integrated_NA_filtered@meta.data[37:41] <- Dextramer_counts[12:16]

# What about double Positive cells? --> Kick them out --> Then the frequency column has to be updated!
Integrated_NA_filtered <- subset(Integrated_NA_filtered, subset = positive.for.n.dextramers<2)
Integrated_NA_filtered@meta.data$Row.names <- rownames(Integrated_NA_filtered@meta.data)
Integrated_NA_filtered@meta.data <- Integrated_NA_filtered@meta.data %>% group_by(Dataset, CTaa) %>% add_count(CTaa,name = "updated.frequency") %>% ungroup()
Integrated_NA_filtered@meta.data$Frequency <- Integrated_NA_filtered@meta.data$updated.frequency
Integrated_NA_filtered@meta.data <- as.data.frame(Integrated_NA_filtered@meta.data)
rownames(Integrated_NA_filtered@meta.data) <- Integrated_NA_filtered@meta.data$Row.names

# Scatter plot examples
FeatureScatter(Integrated_NA_filtered, feature1 = "ATEGALNTPK", feature2 = "negative-Ctrl-A01-01", pt.size = 1,group.by = 'ATEGALNTPK.classification')+
  ggtitle("ATEGALNTPK classification of cells")+xlab("Dex-ATEGALNTPK UMI counts")+ylab("Dex-control UMI counts")+
  labs(color = "CoV2-Dex classification")
FeatureScatter(Integrated_NA_filtered, feature1 = "ATEGALNTPK", feature2 = "negative-Ctrl-A01-01", pt.size = 1,group.by = 'ATEGALNTPK.classification')+xlim(0,30)+
  ggtitle("ATEGALNTPK classification of cells")+xlab("Dex-ATEGALNTPK UMI counts")+ylab("Dex-control UMI counts")+
  labs(color = "CoV2-Dex classification")
FeatureScatter(Integrated_NA_filtered, feature1 = "ATEGALNTPK", feature2 = "KTFPPTEPK", pt.size = 1,group.by = 'positive.for.n.dextramers')+
  labs(color = "CoV2-Dex classification")+ggtitle("ATEGALNTPK vs. KTFPPTEPK counts of cells")+
  scale_color_manual(labels = c("CoV2-Dex-", "CoV2-Dex+"), values = c("#F8766D", "#00BFC4"))

remove(Dextramer_counts)

# What about TCRs that appear to be specific for more than one peptide?
# First, how many of them do we have?
colnames(Integrated_NA_filtered@meta.data)
Investigation <- Integrated_NA_filtered@meta.data[,c(15,37:41)]
Investigation <- Investigation[Investigation$FTSDYYQLY.classification=="Positive"|Investigation$TTDPSFLGRY.classification=="Positive"|Investigation$ATEGALNTPK.classification=="Positive"|Investigation$KTFPPTEPK.classification=="Positive",]
Investigation <- unique(Investigation)
Investigation$CTaa[duplicated(Investigation$CTaa)] 
remove(Investigation)

# This tells us that TCRs "CAVPKNTGNQFYF_CAISVGNEQFF" and "CALSEGNNDMRF_CATSRGLASTDTQYF" are both Positive for more than one peptide.
# Which peptides are that?
colnames(Integrated_NA_filtered@meta.data)
unique(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa=="CAVPKNTGNQFYF_CAISVGNEQFF",c(37:42)])
unique(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa=="CALSEGNNDMRF_CATSRGLASTDTQYF",c(37:42)])

# We see that both are Positive for TTDPSFLGRY and for ATEGALNTPK
# I am using Featurescatter plots to investigate further.
df <- subset(Integrated_NA_filtered, subset = CTaa =="CAVPKNTGNQFYF_CAISVGNEQFF")
FeatureScatter(df, feature1 = "TTDPSFLGRY", feature2 = "ATEGALNTPK", pt.size = 2,group.by = 'positive.for.n.dextramers')+ggtitle("TCR CAVPKNTGNQFYF_CAISVGNEQFF specificity")

df <- subset(Integrated_NA_filtered, subset = CTaa =="CALSEGNNDMRF_CATSRGLASTDTQYF")
FeatureScatter(df, feature1 = "TTDPSFLGRY", feature2 = "ATEGALNTPK", pt.size = 2,group.by = 'positive.for.n.dextramers')+ggtitle("TCR CALSEGNNDMRF_CATSRGLASTDTQYF specificity")
remove(df)

# From this analyis we decided to manually curate all cells with TCR CAVPKNTGNQFYF_CAISVGNEQFF as being positive for no peptide
# Furthermore, we decided to manually curate cells with TCR CALSEGNNDMRF_CATSRGLASTDTQYF as being positive for peptide ATEGALNTPK
Integrated_NA_filtered@meta.data$positive.for.n.dextramers[Integrated_NA_filtered@meta.data$CTaa=="CAVPKNTGNQFYF_CAISVGNEQFF"] <- 0
Integrated_NA_filtered@meta.data$TTDPSFLGRY.classification[Integrated_NA_filtered@meta.data$CTaa=="CAVPKNTGNQFYF_CAISVGNEQFF"] <- "Negative"
Integrated_NA_filtered@meta.data$ATEGALNTPK.classification[Integrated_NA_filtered@meta.data$CTaa=="CAVPKNTGNQFYF_CAISVGNEQFF"] <- "Negative"
Integrated_NA_filtered@meta.data$positive.for.n.dextramers[Integrated_NA_filtered@meta.data$CTaa=="CALSEGNNDMRF_CATSRGLASTDTQYF" & Integrated_NA_filtered@meta.data$TTDPSFLGRY.classification=="Positive"] <- 0
Integrated_NA_filtered@meta.data$TTDPSFLGRY.classification[Integrated_NA_filtered@meta.data$CTaa=="CALSEGNNDMRF_CATSRGLASTDTQYF"] <- "Negative"

# Now the investigation from above should not find anything anymore.
colnames(Integrated_NA_filtered@meta.data)
Investigation <- Integrated_NA_filtered@meta.data[,c(15,37:41)]
Investigation <- Investigation[Investigation$FTSDYYQLY.classification=="Positive"|Investigation$TTDPSFLGRY.classification=="Positive"|Investigation$ATEGALNTPK.classification=="Positive"|Investigation$KTFPPTEPK.classification=="Positive",]
Investigation <- unique(Investigation)
Investigation$CTaa[duplicated(Investigation$CTaa)] 
remove(Investigation)

# Adding the information for which percentage a clone is Dextramer positive.
# Group by clonotype, count cells in that clonotype and cells in that clonotype that have "positive.for.n.dextramers">0 and report percentage of those.
# This has to be done timepoint by timepoint
df_Acute <- as.data.frame(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$Timepoint=="Acute",])
df_Recovered <- as.data.frame(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$Timepoint=="Recovered",])

dt_Acute <- data.table(df_Acute)
dt_Recovered <- data.table(df_Recovered)

dt_Acute <- dt_Acute[,.(dextramer.sum=sum(positive.for.n.dextramers)),by=CTaa] 
dt_Recovered <- dt_Recovered[,.(dextramer.sum=sum(positive.for.n.dextramers)),by=CTaa] 

# Perform first the join of the created data frames
df_Acute <- inner_join(df_Acute,dt_Acute,by = "CTaa")
df_Recovered <- inner_join(df_Recovered,dt_Recovered,by = "CTaa")

rownames(df_Acute) <- df_Acute$Row.names
rownames(df_Recovered) <- df_Recovered$Row.names

df_Acute <- df_Acute[,c(1,43)]
df_Recovered <- df_Recovered[,c(1,43)]

# Perform the second join to the integrated dataset
Integrated_NA_filtered@meta.data <- merge(Integrated_NA_filtered@meta.data,df_Acute,by=0,all.x = TRUE)
rownames(Integrated_NA_filtered@meta.data) <- Integrated_NA_filtered@meta.data$Row.names
colnames(Integrated_NA_filtered@meta.data)
Integrated_NA_filtered@meta.data <- Integrated_NA_filtered@meta.data[,c(1,3:43,45)]

head(Integrated_NA_filtered@meta.data)
colnames(df_Recovered) <- c("Row.names.recovered","dextramer.sum.recovered")
Integrated_NA_filtered@meta.data <- merge(Integrated_NA_filtered@meta.data,df_Recovered,by=0,all.x = TRUE)
rownames(Integrated_NA_filtered@meta.data) <- Integrated_NA_filtered@meta.data$Row.names

Integrated_NA_filtered@meta.data <- Integrated_NA_filtered@meta.data[,c(1,3:44,46)]
Integrated_NA_filtered@meta.data$dextamer_sum <- coalesce(Integrated_NA_filtered@meta.data$dextramer.sum,Integrated_NA_filtered@meta.data$dextramer.sum.recovered)
Integrated_NA_filtered@meta.data <- Integrated_NA_filtered@meta.data[,c(1:42,45)]

remove(df_Acute,df_Recovered,dt_Acute,dt_Recovered)

Integrated_NA_filtered@meta.data$percent.of.clonotype.positive <- 100*(Integrated_NA_filtered@meta.data$dextamer_sum/Integrated_NA_filtered@meta.data$Frequency)
Integrated_NA_filtered@meta.data$percent.with.values <- paste(100*(Integrated_NA_filtered@meta.data$dextamer_sum/Integrated_NA_filtered@meta.data$Frequency),"(",Integrated_NA_filtered@meta.data$dextamer_sum,"of",Integrated_NA_filtered@meta.data$Frequency,")")

df <- Integrated_NA_filtered@meta.data
only_positives <- df[df$positive.for.n.dextramers>0,]
colnames(only_positives)
only_positives <- only_positives[,c(1,13,14,15,17,19,20,21,41,44,45)]
write.xlsx(only_positives, 'percentages of clone that is dextramer positive.xlsx')
remove(df,only_positives)


# I enable the option to manually classify clones as Dextramer negative for each single Dextramer and also put the "positive.for.n.dextramers" to 0
# This can be done to filter clones with a chosen "percent.of.clonotype.positive" cutoff.
# Default is no filtering.
sum(Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0)
table(Integrated_NA_filtered@meta.data$percent.of.clonotype.positive[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0])

Integrated_NA_filtered@meta.data$FTSDYYQLY.classification[Integrated_NA_filtered@meta.data$percent.of.clonotype.positive<0.1 &
                                                             Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0] <- "Negative"
Integrated_NA_filtered@meta.data$TTDPSFLGRY.classification[Integrated_NA_filtered@meta.data$percent.of.clonotype.positive<0.1 &
                                                            Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0] <- "Negative"
Integrated_NA_filtered@meta.data$ATEGALNTPK.classification[Integrated_NA_filtered@meta.data$percent.of.clonotype.positive<0.1 &
                                                              Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0] <- "Negative"
Integrated_NA_filtered@meta.data$KTFPPTEPK.classification[Integrated_NA_filtered@meta.data$percent.of.clonotype.positive<0.1 &
                                                              Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0] <- "Negative"

Integrated_NA_filtered@meta.data$positive.for.n.dextramers[Integrated_NA_filtered@meta.data$percent.of.clonotype.positive<0.1 &
                                                             Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0] <- 0

sum(Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0)
table(Integrated_NA_filtered@meta.data$percent.of.clonotype.positive[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0],
      Integrated_NA_filtered@meta.data$Frequency[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0])



Integrated@meta.data$Row.names <- rownames(Integrated@meta.data)

# At this point, the preliminary analysis is finished. We use the Seurat objects from this point on to prepare the remaining plots.
# We enable here the option to save and load the  Seurat object for easy loading
#SaveH5Seurat(Integrated_NA_filtered, overwrite = TRUE)
#Integrated_NA_filtered <- LoadH5Seurat("SeuratProject.h5Seurat")



###########################################################################################################################################
# Part 7 - Creating plots for the analysis of COVID specific cells
###########################################################################################################################################
# Keep in mind from here on, that there are two Seurat objects: "Integrated" and "Integrated_NA_filtered".
# Most analysis can be done in "Integrated_NA_filtered" alone. 
# But if we want to highlight cells on the UMAP of the full, integrated dataset, we will use "Integrated".

# We now progress to highlighting Dextramer positive cells in acute and recovered disease
Acute_subset <- subset(Integrated_NA_filtered, subset = Timepoint=="Acute")
Recovered_subset <- subset(Integrated_NA_filtered, subset = Timepoint=="Recovered")
acute_positives <- rownames(Acute_subset@meta.data[Acute_subset@meta.data$positive.for.n.dextramers>0,])
recovered_positives <- rownames(Recovered_subset@meta.data[Recovered_subset@meta.data$positive.for.n.dextramers>0,])


DimPlot(Integrated, label=F,repel = T, group.by="seurat_clusters", cells.highlight= acute_positives, cols.highlight = 'palevioletred4', cols= "grey92")+
  ggtitle("SARS-CoV-2 peptide specific cells in acute disease")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),legend.position="none")+xlim(-10, 10)+ylim(-10, 15)+labs(x = "UMAP 1", y="UMAP 2")

DimPlot(Integrated, label=F,repel = T, group.by="seurat_clusters", cells.highlight= recovered_positives, cols.highlight = 'steelblue4', cols= "grey92")+
  ggtitle("SARS-CoV-2 peptide specific cells in recovered disease")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),legend.position="none")+xlim(-10, 10)+ylim(-10, 15)+labs(x = "UMAP 1", y="UMAP 2")



remove(Acute_subset,acute_positives,Recovered_subset,recovered_positives)



# A bar plot showing the absolute number of positive cells per cluster and timepoint
df <- Integrated_NA_filtered@meta.data
df <- df[df$Timepoint!="Healthy",]
df <- df %>% group_by(seurat_clusters,Timepoint) %>% add_tally(positive.for.n.dextramers,name = "positive.cells.percluster.and.timepoint") %>% ungroup()
ggplot(df, aes(x=seurat_clusters,y=positive.cells.percluster.and.timepoint, fill=Timepoint))+geom_bar(position="dodge", stat="identity",colour="black")+theme_classic()+
  scale_fill_manual(values=c("palevioletred4", "steelblue4"))+ scale_y_continuous(expand = c(0, 0), limits = c(0,110))+labs(x = "Clusters", y="SARS-CoV-2 peptide specific cells")+
  theme(text = element_text(size = 15))+ggtitle("Numbers of Dextramer positive cells per cluster")

remove(df)

# Do the same as percentage plot
df2 <- Integrated_NA_filtered@meta.data
df2 <- df2[df2$Timepoint!="Healthy",]
df2 <- df2 %>% group_by(seurat_clusters) %>% add_count(seurat_clusters,name = "cells.percluster.and.timepoint") %>% ungroup()
df2 <- df2 %>% group_by(seurat_clusters,Timepoint) %>% add_tally(positive.for.n.dextramers,name = "positive.cells.percluster.and.timepoint") %>% ungroup()
df2$percent.pos.cells.per.cluster <- round(100*(df2$positive.cells.percluster.and.timepoint/df2$cells.percluster.and.timepoint),1)

ggplot(df2, aes(x=seurat_clusters,y=percent.pos.cells.per.cluster, fill=Timepoint))+geom_bar(position="dodge", stat="identity",colour="black")+theme_classic()+
  scale_fill_manual(values=c("palevioletred4", "steelblue4"))+ scale_y_continuous(expand = c(0, 0), limits = c(0,50))+labs(x = "Clusters", y="% SARS-CoV-2 peptide specific cells")+
  theme(text = element_text(size = 15))+ggtitle("Perentage of Dextramer positive cells per cluster")


remove(df2)

# Create a plot showing in which clusters Dextramer positive cells are in the Acute and Recovered phase.
Contributions <- subset(Integrated_NA_filtered, subset = positive.for.n.dextramers > 0)
Contributions@meta.data$Timepoint <- droplevels(Contributions@meta.data$Timepoint)
cluster.contributions.df <- as.data.frame(with(Contributions@meta.data, table(seurat_clusters, Timepoint)))
cluster.contributions.df <- cluster.contributions.df %>% group_by(Timepoint) %>% add_tally(Freq,name = "total.cells.pertimepoint") %>% ungroup()
cluster.contributions.df <- cluster.contributions.df %>% mutate(percentage.contribution= 100*(Freq/total.cells.pertimepoint))
cluster.contributions.df$percentage.contribution <- round(cluster.contributions.df$percentage.contribution,digits = 1)
cluster.contributions.df$percentage.contribution <- as.numeric(cluster.contributions.df$percentage.contribution)

cluster.contributions.df$Timepoint <- factor(cluster.contributions.df$Timepoint, levels = c("Recovered","Acute","Healthy"))

ggplot(cluster.contributions.df, aes(fill=factor(seurat_clusters, levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)), y=Timepoint, x=percentage.contribution)) + 
  geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ theme(text = element_text(size = 15))+
  labs(x = "Cluster contribution", y="Timepoints")+scale_x_continuous(labels=scales::percent)+labs(fill = "Cluster")+
  scale_fill_manual(values=c(coul2))+ theme(text = element_text(size = 17))+ggtitle("Distribution of CoV2-Dex+ cells among clusters")+ylab("")

remove(cluster.contributions.df,Contributions)



###########################################################################################################################################
# Part 8 - Creating Donut charts for clonal analysis
###########################################################################################################################################
# The next section creates donut charts showing clone sizes and number of clones for each peptide
# The section can be used to create plots for all peptides, just the peptide need to selected below. Don't forget to change the plot title!
# We are including clones that are found as being Dextramer positive in one disease state, even if they appear as Dex specific only in the other disease state.
Acute_subset <- subset(Integrated_NA_filtered, subset = Timepoint=="Acute")
Recovered_subset <- subset(Integrated_NA_filtered,subset = Timepoint=="Recovered")
Dex.positive.CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$KTFPPTEPK.classification=="Positive"])
Acute.subset.pos.CTaas <- subset(Acute_subset,subset = CTaa %in% Dex.positive.CTaas) 
Recovered.subset.pos.CTaas <- subset(Recovered_subset, subset = CTaa %in% Dex.positive.CTaas)

Pie_subset_Acute_df <- as.data.frame(table(Acute.subset.pos.CTaas@meta.data$CTaa))
Pie_subset_Recovered_df <- as.data.frame(table(Recovered.subset.pos.CTaas@meta.data$CTaa))

colnames(Pie_subset_Acute_df) <- c("Clonotype","Freq")
colnames(Pie_subset_Recovered_df) <- c("Clonotype","Freq")

Acute_CTaa_df <- Pie_subset_Acute_df
Recovered_CTaa_df <- Pie_subset_Recovered_df

# Now create the donut
# Compute percentages
Acute_CTaa_df$fraction <- Acute_CTaa_df$Freq / sum(Acute_CTaa_df$Freq)
Recovered_CTaa_df$fraction <- Recovered_CTaa_df$Freq / sum(Recovered_CTaa_df$Freq)

# Compute the cumulative percentages (top of each rectangle)
Acute_CTaa_df$ymax = cumsum(Acute_CTaa_df$fraction)
Recovered_CTaa_df$ymax = cumsum(Recovered_CTaa_df$fraction)

# Compute the bottom of each rectangle
Acute_CTaa_df$ymin = c(0, head(Acute_CTaa_df$ymax, n=-1))
Recovered_CTaa_df$ymin = c(0, head(Recovered_CTaa_df$ymax, n=-1))

Acute_CTaa_df$percentage <- Acute_CTaa_df$fraction*100
Recovered_CTaa_df$percentage <- Recovered_CTaa_df$fraction*100

Acute_CTaa_df$rounded <- round(Acute_CTaa_df$percentage,digits = 1)
Recovered_CTaa_df$rounded <- round(Recovered_CTaa_df$percentage,digits = 1)

# Compute label position
Acute_CTaa_df$labelPosition <- (Acute_CTaa_df$ymax + Acute_CTaa_df$ymin) / 2
Recovered_CTaa_df$labelPosition <- (Recovered_CTaa_df$ymax + Recovered_CTaa_df$ymin) / 2

Acute_CTaa_df$Legend <- paste0(Acute_CTaa_df$Clonotype," (",Acute_CTaa_df$rounded,"%)")
Recovered_CTaa_df$Legend <- paste0(Recovered_CTaa_df$Clonotype," (",Recovered_CTaa_df$rounded,"%)")

# Make the plot
ggplot(Acute_CTaa_df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Legend)) +
  geom_rect(color='black') +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) + theme_void()+
  geom_text(x = 2, y = 0.5, label = nrow(Acute_CTaa_df), size = 10)+guides(
    fill = guide_legend(
      title = "Clonotypes",
      override.aes = aes(label = "")))+ggtitle("Acute \n KTFPPTEPK")+
  theme(plot.title = element_text(hjust = 0.5,size=22),legend.position="none")

ggplot(Recovered_CTaa_df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Legend)) +
  geom_rect(color='black')  +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) + theme_void()+
  geom_text(x = 2, y = 0.5, label = nrow(Recovered_CTaa_df), size = 10)+guides(
    fill = guide_legend(
      title = "Clonotypes",
      override.aes = aes(label = "")))+ggtitle("Recovered \n KTFPPTEPK")+
  theme(plot.title = element_text(hjust = 0.5,size=22),legend.position="none")



remove(Acute_CTaa_df, Acute_subset,Acute.subset.pos.CTaas,Pie_subset_Acute,Pie_subset_Acute_df,Pie_subset_Recovered,Pie_subset_Recovered_df,
       Recovered_CTaa_df,Recovered_subset,Recovered.subset.pos.CTaas,Dex.positive.CTaas)






###########################################################################################################################################
# Part 9 - Comparing persistent vs non-persistent clones
###########################################################################################################################################
# The first plot creates a piechart that shows all acute clones and whether they are also found in the 6 month timepoint or not
Acute_subset <- subset(Integrated_NA_filtered, subset = Timepoint=="Acute")
unique_acute_CTaas <- unique(Acute_subset@meta.data$CTaa)
Acute_commons <- Acute_subset@meta.data[Acute_subset@meta.data$common.clone=="yes",]
unique_acute_commons <- unique(Acute_commons$CTaa)

persistent_df <- data.frame(
  Occurance=c("persistent","non-persistent"),
  count=c(length(unique_acute_commons),(length(unique_acute_CTaas)-length(unique_acute_commons))))

persistent_df <- persistent_df %>% 
  arrange(desc(Occurance)) %>%
  mutate(prop = count / sum(persistent_df$count) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

ggplot(persistent_df, aes(x="", y=prop, fill=Occurance)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  geom_text(aes(y = ypos, label = count), color = "white", size=6) +
  scale_fill_manual(values = c("dodgerblue3", "deeppink4")) +ggtitle("Persistent vs. non-persistent clones \n (All unique clones)")+
  theme(plot.title = element_text(hjust = 0.5,size=22),text = element_text(size = 26))

remove(Acute_commons,Acute_subset,persistent_df,unique_acute_commons,unique_acute_CTaas)


# Create the same plot only for COVID specific clones (clones that are specific in acute phase but not necessarily specific in recovered phase)!
# New, include also those that are found to be positive in the recovered phase but not necessarily positive in the acute phase. 
Acute_but_non_common_CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0 & 
                                                                      Integrated_NA_filtered@meta.data$Timepoint=="Acute" &
                                                                        Integrated_NA_filtered@meta.data$common.clone=="no"]) 
Common_positive_CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0 &
                                                              Integrated_NA_filtered@meta.data$common.clone=="yes"])

persistent_df2 <- data.frame(
  Occurance=c("persistent","non-persistent"),
  count=c(length(Common_positive_CTaas),length(Acute_but_non_common_CTaas)))

persistent_df2 <- persistent_df2 %>% 
  arrange(desc(Occurance)) %>%
  mutate(prop = count / sum(persistent_df2$count) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

ggplot(persistent_df2, aes(x="", y=prop, fill=Occurance,colo)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  geom_text(aes(y = ypos, label = count), color = "white", size=12) +
  scale_fill_manual(values = c("dodgerblue3", "deeppink4")) +ggtitle("Persistent vs. non-persistent clones \n (Specific unique clones in Acute phase)")+
  theme(plot.title = element_text(hjust = 0.5,size=22),text = element_text(size = 26))

remove(persistent_df2,Acute_but_non_common_CTaas,Common_positive_CTaas)


# Highlight specific persistent and non persistent clones in the UMAP
Acute_but_non_common_CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0 & 
                                                                             Integrated_NA_filtered@meta.data$Timepoint=="Acute" &
                                                                             Integrated_NA_filtered@meta.data$common.clone=="no"]) 
Common_positive_CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0 &
                                                                        Integrated_NA_filtered@meta.data$common.clone=="yes"])

Non_persistents <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa %in% Acute_but_non_common_CTaas& 
                                                               Integrated_NA_filtered@meta.data$Timepoint=="Acute",])
Persistents <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa %in% Common_positive_CTaas & 
                                                           Integrated_NA_filtered@meta.data$Timepoint=="Acute",])


list <- list(Persistents,Non_persistents)
names(list) <- c("Persistent", "Non-persistent")
DimPlot(Integrated, label=F,repel = T, group.by="seurat_clusters", cells.highlight= list, cols.highlight = c("deeppink4", "dodgerblue3"), cols= "grey92")+
  labs(x = "UMAP 1", y="UMAP 2")+xlim(-10, 10)+ylim(-10, 15)+ggtitle("")

remove(list,Acute_but_non_common_CTaas,Common_positive_CTaas)


# To which clusters do the persistent and non persistent clones belong?
ps_and_nps <- subset(Integrated_NA_filtered,subset = Row.names %in% Persistents |Row.names %in% Non_persistents)
ps_and_nps@meta.data <- ps_and_nps@meta.data %>% group_by(common.clone,seurat_clusters) %>% add_count(seurat_clusters,name = "cells.per.cluster") %>% ungroup()

sum(ps_and_nps@meta.data$common.clone=="yes")
sum(ps_and_nps@meta.data$common.clone=="no")

ps_and_nps@meta.data$percent <- 0
ps_and_nps@meta.data$percent[ps_and_nps@meta.data$common.clone =="yes"] <- round(100*(ps_and_nps@meta.data$cells.per.cluster[ps_and_nps@meta.data$common.clone =="yes"]/sum(ps_and_nps@meta.data$common.clone=="yes")),1)
ps_and_nps@meta.data$percent[ps_and_nps@meta.data$common.clone =="no"] <- round(100*(ps_and_nps@meta.data$cells.per.cluster[ps_and_nps@meta.data$common.clone =="no"]/sum(ps_and_nps@meta.data$common.clone=="no")),1)

df <- ps_and_nps@meta.data
colnames(df)
df <- df[,c(19,35,47)]
df <- unique(df)
ggplot(df, aes(x=seurat_clusters,y=percent, fill=common.clone))+geom_bar(position="dodge", stat="identity",colour="black")+theme_classic()+
  scale_fill_manual(values=c("palevioletred4", "steelblue4"))+ scale_y_continuous(expand = c(0, 0), limits = c(0,50))+labs(x = "Clusters", y="% SARS-CoV-2 peptide specific cells")+
  theme(text = element_text(size = 15))+labs(fill='Persistent clone')+ggtitle("Distribution of cells belonging to Persistent \n or Non-persistent clones across clusters")

ggplot(df, aes(fill=factor(seurat_clusters, levels=c(0,1,2,3,4,5,6,7,8,9,10,11)), y=common.clone, x=percent)) + 
  geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ theme(text = element_text(size = 15))+
  labs(x = "", y="Persistent")+scale_x_continuous(labels=scales::percent)+labs(fill = "Cluster")+
  scale_fill_manual(values=c(coul2))+ theme(text = element_text(size = 17))

remove(df,ps_and_nps)


# Violin plot comparing clone sizes specific vs non-specific
Violin.plot.persistents <- Integrated_NA_filtered@meta.data[Persistents,]
Violin.plot.nonpersistents <- Integrated_NA_filtered@meta.data[Non_persistents,]

colnames(Violin.plot.persistents)
colnames(Violin.plot.nonpersistents)
Violin.plot.persistents$occurence <- "Persistents"
Violin.plot.nonpersistents$occurence <- "Non-persistents"

Violin.plot.df <- rbind(Violin.plot.persistents,Violin.plot.nonpersistents)
colnames(Violin.plot.df)
Violin.plot.df <- Violin.plot.df[,c(15,17,46)]
Violin.plot.df <- unique(Violin.plot.df)

ggplot(Violin.plot.df, aes(x=occurence, y=Frequency, fill=occurence)) + 
  geom_violin()+theme_classic()+geom_jitter(shape=16, position=position_jitter(0.0))+
  scale_fill_manual(values=c("dodgerblue3", "deeppink4"))+ labs(x="",y="Clonal Frequency")+ theme(legend.title = element_blank(),text = element_text(size = 20))+
  scale_y_continuous(trans='log2')

remove(Violin.plot.df,Violin.plot.nonpersistents,Violin.plot.persistents)



# Highlighting the 10 most expanded persistent and nonpersistent clones in the UMAP -> This has to be adapted based on the cutoff chosen in Part 6.
# First persistent
Specifics_persistent <- Integrated_NA_filtered
Specifics_persistent@meta.data$Row.names <- rownames(Specifics_persistent@meta.data)
Specifics_persistent@meta.data$Specific <- "no"
Specifics_persistent@meta.data$Specific[rownames(Specifics_persistent@meta.data)%in%Persistents] <- "yes" # this has to be changed depending on the group
Specifics_persistent@meta.data$Clone <- "a"
Specifics_persistent@meta.data$Clone[Specifics_persistent@meta.data$Specific=="yes"] <- Specifics_persistent@meta.data$CTaa[Specifics_persistent@meta.data$Specific=="yes"]
Specifics_persistent <- subset(Specifics_persistent, subset = Specific =="yes")
Specifics_persistent@meta.data$Clone <- as.numeric(as.factor(Specifics_persistent@meta.data$Clone))
Specifics_persistent@meta.data <- Specifics_persistent@meta.data %>% group_by(CTaa) %>% add_count(Clone, name = "fake.Clonesize") %>% ungroup()
Specifics_persistent@meta.data <- as.data.frame(Specifics_persistent@meta.data)
rownames(Specifics_persistent@meta.data) <- Specifics_persistent@meta.data$Row.names
table(Specifics_persistent@meta.data$fake.Clonesize)

First <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==83, ])
Second <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==77, ])
Third <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==40, ])
Fourth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==27, ])
Fifth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==26, ])
Sixth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==22, ])
Seventh <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==16, ])
Eigth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==12, ])
Nineth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==10, ])
Tenth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==9, ])

Integrated@meta.data$Highlight.column <- "AAA"
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%First] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%First])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Second] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Second])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Third] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Third])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Fourth] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Fourth])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Fifth] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Fifth])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Sixth] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Sixth])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Seventh] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Seventh])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Eigth] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Eigth])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Nineth] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Nineth])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Tenth][c(1,3,6,7,8,9,10,12,15)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Tenth])[1]

highlight.colors <- c("grey92",coul2[1:10])
DimPlot(Integrated, reduction = "umap",label = F,pt.size = 0.8, group.by = "Highlight.column",order = T,cols = highlight.colors)+
  labs(x = "UMAP 1", y="UMAP 2")+xlim(-10, 15)+ylim(-10, 15)+ggtitle("Top 10 expanded persistent clones")+theme(text = element_text(size = 10))


# Then the non-persistent
Specifics_persistent <- Integrated_NA_filtered
Specifics_persistent@meta.data$Row.names <- rownames(Specifics_persistent@meta.data)
Specifics_persistent@meta.data$Specific <- "no"
Specifics_persistent@meta.data$Specific[rownames(Specifics_persistent@meta.data)%in%Non_persistents] <- "yes" # this has to be changed depending on the group
Specifics_persistent@meta.data$Clone <- "a"
Specifics_persistent@meta.data$Clone[Specifics_persistent@meta.data$Specific=="yes"] <- Specifics_persistent@meta.data$CTaa[Specifics_persistent@meta.data$Specific=="yes"]
Specifics_persistent <- subset(Specifics_persistent, subset = Specific =="yes")

Specifics_persistent@meta.data$Clone <- as.numeric(as.factor(Specifics_persistent@meta.data$Clone))
Specifics_persistent@meta.data <- Specifics_persistent@meta.data %>% group_by(CTaa) %>% add_count(Clone, name = "fake.Clonesize") %>% ungroup()
table(Specifics_persistent@meta.data$fake.Clonesize)
Specifics_persistent@meta.data <- as.data.frame(Specifics_persistent@meta.data)
rownames(Specifics_persistent@meta.data) <- Specifics_persistent@meta.data$Row.names
table(Specifics_persistent@meta.data$fake.Clonesize)

First <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==11, ])
Second <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==9, ])
Third <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==5, ])
Fourth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==5, ])
Fifth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==4, ])
Sixth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==4, ])
Seventh <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==4, ])
Eigth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==4, ])
Nineth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==3, ])
Tenth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==3, ])

Integrated@meta.data$Highlight.column <- "AAA"
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%First] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%First])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Second] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Second])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Third][c(1,2,3,4,5)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Third])[1]
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Fourth][c(6,7,8,9,10)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Fourth])[2]
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Fifth][c(1,2,4,7)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Fifth])[1]
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Sixth][c(3,6,10,13)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Sixth])[2]
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Seventh][c(5,8,9,15)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Seventh])[3]
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Eigth][c(11,12,14,16)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Eigth])[4]
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Nineth][c(1,2,3)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Nineth])[1]
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Tenth][c(4,5,6)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Tenth])[2]

highlight.colors <- c("grey92",coul2[1:10])
DimPlot(Integrated, reduction = "umap",label = F,pt.size = 0.8, group.by = "Highlight.column",order = T,cols = highlight.colors)+
  labs(x = "UMAP 1", y="UMAP 2")+xlim(-10, 15)+ylim(-10, 15)+ggtitle("Top 10 expanded non-persistent clones")+theme(text = element_text(size = 10))

remove(Specifics_persistent,First,Second,Third,Fourth,Fifth,Sixth,Seventh,Eigth,Nineth,Tenth,highlight.colors,)


###########################################################################################################################################
# Part 10 - GSEA
###########################################################################################################################################
# We want to compare two groups of cells with each other: Persistent vs Nonpersistent cells
##### Input files 
list <- list(Persistents,Non_persistents)
names(list) <- c("persistent", "nonpersistent")

GSEA_persistent.vs.nonpersistent <- subset(Integrated_NA_filtered, 
                                           subset = Row.names %in% list$persistent | Row.names %in% list$nonpersistent)
GSEA_persistent.vs.nonpersistent$persistence <- "nonpersistent"
GSEA_persistent.vs.nonpersistent$persistence[GSEA_persistent.vs.nonpersistent$Row.names %in% list$persistent] <- "persistent"

##### Choose the collection(s) of gene sets to be tested for enrichment 
# see https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# we choose the Hallmark collection:
m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#### Obtain differential expression results with Seurat's FindMarkers 
# We lower the log2 fold change threshold to 0 because
# we want to account also for small changes,
# as GSEA considers all of the genes in an experiment, 
# not only those above an arbitrary cutoff in terms of fold-change or significance.
DefaultAssay(GSEA_persistent.vs.nonpersistent) <- "RNA"
Idents(GSEA_persistent.vs.nonpersistent) <- "persistence"
p_vs_np.results <- FindMarkers(GSEA_persistent.vs.nonpersistent, 
                               logfc.threshold = -Inf,  # equivalent to using 0
                               min.pct = 0.1, 
                               ident.1 = "persistent", ident.2 = "nonpersistent")

#### Define the ranking metric and obtain the pre-ranked feature list 
# For each feature (i.e., each gene) tested, we define the ranking metric as
# the product between the -log of the nominal p-value and the sign of the log2 fold change.
p_vs_np.results <- p_vs_np.results %>% 
  mutate(ranking_metric = -log(p_val)*sign(avg_log2FC) ) %>%
  arrange(desc(ranking_metric)) %>% rownames_to_column(var = "feature")

preranked_list <- p_vs_np.results %>% arrange(desc(ranking_metric)) %>% 
  dplyr::select(feature, ranking_metric) %>% deframe

#### Run GSEA 
# NB: we set the seed before running (F)GSEA as the n permutations (nperm) are random, and we want to make the results reproducible
# The original pre-ranked GSEA (from Subramanian et al., 2005) used 1000 permutations.
# We decide to use more permutations, as FGSEA is fast and allows for that.
# We also decide to run FGSEA-Simple (with permutations),
# rather than FGSEA-multilevel (default option since v3 of the FGSEA pre-print),
# as it's the method consistent with the reference implementation (from Subramanian et al., 2005).
set.seed(42)
fgseaResults <- fgsea(fgsea_sets, stats = preranked_list, nperm=100000) %>%
  as_tibble() %>%
  arrange(desc(NES))

# To print the results to a tab-delimited file, first restructure the list of genes in the leading edge: 
# fgseaResults %>% mutate(leadingEdge = leadingEdge %>% vapply(., paste, collapse = ", ", character(1L))) %>% write.table("gsea_hallmarks_results.txt", sep="\t", quote=F, row.names = F)

#### Plot the NES of the gene sets 
# Only for those found significantly enriched with
# adjusted p-value below a certain threshold (here, padj < 0.1)
p.nes <- ggplot(fgseaResults %>% filter(padj < 0.1), aes(reorder(pathway, NES), NES)) +
  geom_col(fill="blue") +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways") + theme_cowplot()

p.nes

#### Save the individual gene set plots 
# Select the gene sets to plot (by p-value adjusted)
genesets_to_plot <- fgseaResults %>% filter(padj < 0.1) %>% select(pathway) %>% unlist %>% as.character

pdf("Fig4f.pdf", onefile=TRUE)
for (i in 1:length(genesets_to_plot)){
  print(plotEnrichment(fgsea_sets[[ genesets_to_plot[i] ]], preranked_list) + 
          # make the gene set names more readable and remove the 'HALLMARK' string at their beginning
          labs(title=stringr::str_replace_all(stringr::str_remove(genesets_to_plot[i], pattern = "HALLMARK_"), 
                                              pattern = "_", replacement = " "), hjust=0.5) +
          theme_cowplot() + labs(x = "Rank", y = "Enrichment Score") +
          geom_line(color = "blue") + geom_point(color = "blue", size=0.1) +
          # center the title
          theme(plot.title = element_text(hjust = 0.5))
  )
}
dev.off()

remove(fgsea_sets,fgseaResults,GSEA_persistent.vs.nonpersistent,list,m_df,p_vs_np.results,p.nes,genesets_to_plot,i,preranked_list)

###########################################################################################################################################
# Part 11 - DGEA
###########################################################################################################################################
# DGEA Persistent vs non persistent cells
persistent.vs.non.persistent.markers <- as.data.frame(FindMarkers(Integrated_NA_filtered ,ident.1 = Persistents, ident.2 = Non_persistents, assay = "RNA"))
persistent.vs.non.persistent.markers <- persistent.vs.non.persistent.markers[,c(2,5)]
persistent.vs.non.persistent.markers <- persistent.vs.non.persistent.markers[1:170,]
write.xlsx(persistent.vs.non.persistent.markers, 'persistent.vs.non.persistent.markers.xlsx',row.names = TRUE)

# DGEA for selected genes
features <- c('IL7R','GZMB','CX3CR1','IL2RG','NR4A1','ZAP70','STAT3')
persistent.vs.non.persistent.markers.selected.genes <- as.data.frame(FindMarkers(Integrated_NA_filtered ,ident.1 = Persistents, ident.2 = Non_persistents,
                                                                                 logfc.threshold = 0, assay = "RNA",min.pct = 0,features = features))
persistent.vs.non.persistent.markers.selected.genes <- persistent.vs.non.persistent.markers.selected.genes[,c(2,5)]
write.xlsx(persistent.vs.non.persistent.markers.selected.genes, 'persistent.vs.non.persistent.markers.selected.genes.xlsx',row.names = TRUE)

# FGEA for TotalSeq proteins
features_found.in.Protein <- c('CD45RA','CCR7.1')
persistent.vs.non.persistent.proteins <- as.data.frame(FindMarkers(Integrated_NA_filtered ,ident.1 = Persistents, ident.2 = Non_persistents,
                                                                   logfc.threshold = 0, assay = "Protein",min.pct = 0,features = features_found.in.Protein))
persistent.vs.non.persistent.proteins <- persistent.vs.non.persistent.proteins[,c(2,5)]
write.xlsx(persistent.vs.non.persistent.proteins, 'persistent.vs.non.persistent.proteins.xlsx',row.names = TRUE)

remove(persistent.vs.non.persistent.markers,persistent.vs.non.persistent.markers.selected.genes,persistent.vs.non.persistent.proteins,features,features_found.in.Protein)


###########################################################################################################################################
# Part 12 - Alluvial plots, testing if there are clones which are found in more than one SNP cluster, venn diagram
###########################################################################################################################################
# To create alluvial plots, scRepertoire expects as input:
# - a list
# - the list elements should have the names of the samples you want to compare
# I have changed the function "compareClonotypes' of package "Repertoire" so that instead of relative frequencies, absolute frequencies are shown in the alluvial plot.
# trace(compareClonotypes, edit=TRUE)
# Alluvial plot including all persistent clones
persistent_CTaas <- Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$Row.names %in% Persistents]
persistent_CTaas <-  unique(persistent_CTaas)

Cells_with_persistent_CTaas <- Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa %in% persistent_CTaas,]
Acute_Cells_with_persistent_CTaas <- Cells_with_persistent_CTaas[Cells_with_persistent_CTaas$Timepoint=="Acute",]
Recovered_Cells_with_persistent_CTaas <- Cells_with_persistent_CTaas[Cells_with_persistent_CTaas$Timepoint=="Recovered",]

List <- list(Acute_Cells_with_persistent_CTaas,Recovered_Cells_with_persistent_CTaas)
names(List) <- c("Acute","Recovered")

compareClonotypes(List, samples = c("Acute", "Recovered"), 
                  cloneCall="aa", graph = "alluvial")+ggtitle("SARS-CoV-2 common specific clones")+
  theme(plot.title = element_text(hjust = 0.5))

remove(Acute_Cells_with_persistent_CTaas,Cells_with_persistent_CTaas,List,Recovered_Cells_with_persistent_CTaas,persistent_CTaas)



### Testing if there are clones which are found in more than one SNP clusters
# Split by Set1+2 and Set3+4
head(Integrated_NA_filtered@meta.data)
table(Integrated_NA_filtered@meta.data$Dataset)
table(Integrated_NA_filtered@meta.data$assignment)
Cross_patient_clones <- subset(Integrated_NA_filtered, subset = assignment != 'Healthy' & assignment != 'unknown_1' & assignment != 'unknown_2')
Cross_patient_clones_Set1_2 <- subset(Cross_patient_clones, subset = Dataset=="Set1" |Dataset=="Set2")
Cross_patient_clones_Set3_4 <- subset(Cross_patient_clones, subset = Dataset=="Set3" |Dataset=="Set4")

table(Cross_patient_clones_Set1_2@meta.data$assignment)
table(Cross_patient_clones_Set3_4@meta.data$assignment)

Cross_patient_clones_Set1_2 <- subset(Cross_patient_clones_Set1_2,subset = positive.for.n.dextramers>0)
colnames(Cross_patient_clones_Set1_2@meta.data)
Cross_patient_clones@meta.data <- Cross_patient_clones@meta.data[,c(12,15)]
test <- table(Cross_patient_clones@meta.data$CTaa,Cross_patient_clones@meta.data$assignment)
test <- as.data.frame.matrix(test) 
test$zeros <- unname(rowSums(test==0))
table(test$zeros)



# Make a pie chart out of this
pie_df <- data.frame(
  Assigned.to.n.SNP.clusters=c("1","2",">2"),
  count=c(4075,159,17))

pie_df <- pie_df %>% 
  arrange(desc(Assigned.to.n.SNP.clusters)) %>%
  mutate(prop = count / sum(pie_df$count) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

pie_df$percent <- paste0(pie_df$prop,"%")
pie_df$percent <- c("","95.9%","")
pie_df <- pie_df[order(pie_df$count,decreasing=T),]
pie_df$Assigned.to.n.SNP.clusters <- factor(pie_df$Assigned.to.n.SNP.clusters,levels = c("1", "2", ">2"))
levels(pie_df$Assigned.to.n.SNP.clusters)

ggplot(pie_df, aes(x="", y=prop, fill=Assigned.to.n.SNP.clusters)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  geom_text(aes(x=1, y = ypos, label = percent), color = "black", size=8) + ggtitle("Number of SNP clusters to which a clone is assigned")+
  theme(plot.title = element_text(hjust = 0.5,size=22),text = element_text(size = 20))+ labs(fill = "")

remove(Cross_patient_clones,Cross_patient_clones_Set1_2,Cross_patient_clones_Set3_4,pie_df,test)



### Venn Diagram
Positive_CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0])

Positive_CTaas_in_Acute <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$Timepoint=="Acute"&
                                                                          Integrated_NA_filtered@meta.data$CTaa %in% Positive_CTaas])
Positive_CTaas_in_Recovered <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$Timepoint=="Recovered"&
                                                                              Integrated_NA_filtered@meta.data$CTaa %in% Positive_CTaas])
Venn_list <- list('Acute' = Positive_CTaas_in_Acute,'Recovered' = Positive_CTaas_in_Recovered)



# We have noticed that the number of overlapping clones is higher than expected after previous analysis steps. 41 would be expected, but there seem
# to be 44 overlapping clones between the acute and the recovered phase.
# What about the three overlapping clones that should not be there?
# I suspect they are clones that are shared between acute and common but come from different sets! So they are not common clones.
overlap <- intersect(Venn_list[[1]],Venn_list[[2]])
Investigation <- Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa %in% overlap & Integrated_NA_filtered@meta.data$common.clone=="no",]
unique(Investigation$CTaa) # --> This confirms, that there are 3 non.common clones which add to the 41 "real" overlapping ones.
kicks <- unique(Investigation$CTaa)

# To correct this, I am annotating them in the Positive_CTaas_in_Acute and Positive_CTaas_in_Recovered
Positive_CTaas_in_Acute[grepl(kicks[1],Positive_CTaas_in_Acute)] <- paste0(kicks[1],"_Acute")
Positive_CTaas_in_Acute[grepl(kicks[2],Positive_CTaas_in_Acute)] <- paste0(kicks[2],"_Acute")
Positive_CTaas_in_Acute[grepl(kicks[3],Positive_CTaas_in_Acute)] <- paste0(kicks[3],"_Acute")

Positive_CTaas_in_Recovered[grepl(kicks[1],Positive_CTaas_in_Recovered)] <- paste0(kicks[1],"_Recovered")
Positive_CTaas_in_Recovered[grepl(kicks[2],Positive_CTaas_in_Recovered)] <- paste0(kicks[2],"_Recovered")
Positive_CTaas_in_Recovered[grepl(kicks[3],Positive_CTaas_in_Recovered)] <- paste0(kicks[3],"_Recovered")

Venn_list <- list('Acute' = Positive_CTaas_in_Acute,'Recovered' = Positive_CTaas_in_Recovered)

# Finally, we are creating two different versions of venn digrams, the plot in the figure is a mix of the two.
ggvenn(Venn_list,c('Acute','Recovered'))

Venn_object <- euler(Venn_list)
plot <-plot(Venn_object, alpha=0.6, fill=c("blue", "red"),labels=list(font=1, cex=2))
plot

remove(Investigation, plot,Venn_list,Venn_object,kicks,overlap,Positive_CTaas,Positive_CTaas_in_Acute,Positive_CTaas_in_Recovered)


