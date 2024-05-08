# InGene
Find relevant genes from tSNE and UMAP projections of scRNA-seq  data

---
"InGene Tutorial (R)"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## To install the package 

```{r}
Install devtools or remotes package (if not already installed):

install.packages("devtools")

OR

install.packages("remotes")

Then use either to install the package from github:

devtools::install_github("cgiiitd/InGene/R/InGene")

OR

remotes::install_github("cgiiitd/InGene/R/InGene")
```

## The first step is to load the dataset. Here we are loading the Darmanis dataset, which contains neurons, astrocytes and oligodendrocytes cells after the cell types were filtered out. There are 227 cells left. The data set can be accessed inside the data folder

```{r}

raw_Data = read.table("data/Darmanis/Darmanis_Raw_Count.csv",header=TRUE,sep = ",")
samples = raw_Data[1,][-1]
genes = raw_Data$X
raw_Data$X = NULL
```

## Next, we read the annotations, and match it to the cell types, and drop the duplicate and NA genes as well.

```{r}

annotation = read.table("data/Darmanis/Darmanis_Annotation.csv",header = TRUE, sep = ",")
annot = annotation[,2]
colnames(raw_Data) = annot

index_ToDrop = which(duplicated(genes) == TRUE | is.na(genes) == TRUE)
genes = genes[-index_ToDrop]
raw_Data = raw_Data[-index_ToDrop,]
```

## Creating a single cell experiment object out of the raw data. The steps involved are filtering out poorly expressed cells and genes, normalizing and scaling the data.

```{r}
sce_Darmanis <- preprocess_Data(raw_data = raw_Data,gene_list = genes,
                                annot = annot,nfeatures="all", min_Cell = 0.1, min_Reads=0.01, min_Gene = 200)

```

## Next we create the 2D nonlinear embeddings. Here we are using UMAP. 
## Plot the 2D embeddings to visualize the cell types and clusters

```{r}

DR_umap_Darmanis = dim_Red(sce  = sce_Darmanis ,annot = sce_Darmanis$orig.ident, num_threads = 1,resolution = 0.3,reduction = 'umap',verbose = TRUE, perplex=perplex, n.neighbors = 30, uncert = 0.5, spread = 1.5)


Darmanis_umap_allGenes_labels = vis_2D_truth(DR_umap_Darmanis,x='UMAP 1',y='UMAP 2') 
Darmanis_umap_allGenes_clusters = vis_2D_pred(DR_umap_Darmanis,x='UMAP 1',y='UMAP 2') 

plot(Darmanis_umap_allGenes_labels) #umap with the ground truth labels
plot(Darmanis_umap_allGenes_clusters) #umap with the predicetd clusters


```

## Selecting representative cells from each cluster

```{r}

cl_res_conf_umap_Darmanis = rele_Cells(dim_Red_List = DR_umap_Darmanis,
                                             sce = sce_Darmanis,
                                             min_Cell=20, n_Trees = 1000,num_threads = 12,
                                             choose_random_cells = TRUE)

```

## Ranking the genes using the representative cells. 

```{r}

ori_conf_1000_umap_Darmanis = getHighestCorrelation(sce =  sce_Darmanis,reduced = DR_umap_Darmanis, cl_res = cl_res_conf_umap_Darmanis, n_gene_seq = seq(10,1000,10),
                                                    sample_cells=TRUE)



top_Genes = getRelevantGenes(cl_res_conf_umap_Darmanis, nTop =500, ori = ori_conf_1000_umap_Darmanis, genes = rownames(sce_Darmanis)) #Selecting the top 500 ranked genes for further analysis
```

## Plots.  UMAP with all genes (A) ---> UMAP with top 1000 InGene genes (B)
(A) ![UMAP with all genes](/Plots/MB_AllGenesUMAP.png)   (B)  ![UMAP with all genes](/Plots/MB_InGene_UMAP.png)   

