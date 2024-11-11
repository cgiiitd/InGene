#' Filter out poor quality cells and genes
#' @description Filter out poor quality cells and genes using given thresholds
filter_cells_genes <- function (raw_Data, min_Reads = 1, min_Cell = 0.1, min_Gene = 1000,
                                log = T)
{
  
  #cell Filtering
  MIN = min(raw_Data)
  good_Cells = apply(raw_Data, 2, function(x) sum(x > MIN)) >= min_Gene
  temp_Data = raw_Data[, good_Cells]
  
  #Gene filtering
  C = floor(dim(temp_Data)[2] * min_Cell)
  exprs_Genes = apply(temp_Data, 1, function(x) sum(x > min_Reads)) >= C
  temp_Data = temp_Data[exprs_Genes, ]
  
  cat(paste("Remaining genes:", dim(temp_Data)[1], "\n", sep = ""))
  cat(paste("Remaining cells:", dim(temp_Data)[2],  "\n", sep = ""))
  
  return(temp_Data)
}


#' Preprocess the raw data
#' @description Filters out poor quality cells and genes by calling the \code{filter_cells_genes} function.Then normalizes and scales the data.
#' @import SingleCellExperiment
#' @import Seurat
#' @param raw_data Count or FPKM data for single cell transcriptomes (genes should be on rows)
#' @param min_Reads A gene is considered detected if it has at least \code{min_Reads} number of reads in certain percentage of all cells
#' @param min_Cell Minimum \code{min_Cells} fraction of all cells should have at least \code{min_Reads} number of reads to be considered as reliably expressed
#' @param min_Gene A cell should have at least \code{min_Gene} number of expressed genes
#' @param gene_list the gene names for the dataset
#' @param annot the cell labels/annotations/ground truth. If it is unknown, then consecutive numbers starting from 1 wil be used as cell names.
#' @param nfeatures The number of variable features that Seurat will use to scale the normalized the data. It can be a number or "all", in which case, all genes after filtering will be retained. We recommend keeping the value as "all".
#' @return SingleCellExperiment object with the bad cells removed
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' counts <- matrix(rpois(5000000, lambda = 10), ncol=5000, nrow=5000)
#' genes <- 1:nrow(counts)
#' annot <- paste0("CellType", rep(1:2, each = 2500)) 
#' sce <- preprocess_Data(counts, gene_List = genes, 
#' annot = annot, nfeatures = "all",  
#' min_Cell = 0.05, 
#' min_Reads=0.001, 
#' min_Gene = 1000)

preprocess_Data <- function(raw_data,gene_list,annot, nfeatures="all",
                            min_Reads = 5, min_Cell = 0.05, min_Gene = 1000)
{
  
  # Percentage of cells in which a gene should be expressed to be considered for further analysis
  cells_before_filtering = 1:ncol(raw_data)
  colnames(raw_data) = cells_before_filtering
  rownames(raw_data) <- gene_list
  
  
  if(min_Reads != 0){
    raw_data = filter_cells_genes(raw_Data = raw_data,min_Reads = min_Reads,
                                  min_Cell = min_Cell, min_Gene = min_Gene)
  }
  else{stop("Must enter a minimum threshold to remove bad cells. If you dont want to remove any cells please enter 0 for min_Reads")}
  
  common_ind = match(colnames(raw_data), cells_before_filtering)
  
  gene_list <- rownames(raw_data)
  annot = annot[common_ind]
  
  sce <- SingleCellExperiment(assays = list(counts = raw_data)) #transpose it to cells*genes
  normcounts(sce) = assay(sce)
  
  sce = CreateSeuratObject(counts = (counts(sce)))
  
  
  sce <- SingleCellExperiment(assays = list(counts = raw_data)) #transpose it to cells*genes
  normcounts(sce) = assay(sce)
  rownames(sce) <- gene_list
  
  number_genes = length(gene_list)
  
  sce = CreateSeuratObject(counts = (counts(sce)))
  sce$orig.ident <- annot
  
  
  sce <- NormalizeData(object = sce)
  if(min_Reads == 0){
    sce[["RNA"]]@data = raw_data
  }
  
  
  if(nfeatures != "all"){
    sce <- FindVariableFeatures(object = sce,nfeatures=min(number_genes, nfeatures))
  }
  
  else{
    sce <- FindVariableFeatures(object = sce,nfeatures=nrow(sce))
  }
  
  sce <- ScaleData(object = sce)
  sce$orig.ident <- annot
  sce <- RenameCells(sce, new.names = c(1:ncol(sce)))
  
  return(sce)
}

#'Fits GMM distribution to each cluster. Used by \code{dim_Red} function.

fitGMM <- function(embeddings,seed = 222, annot = FALSE, uncert_thresh = 0.1)
{
  
  set.seed(seed)
  xyMclust <- MclustDA(data.frame(embeddings[, 1], embeddings[, 2]),
                       class = embeddings$leiden_clusters)
  
  probs = c()
  ind = c()
  pred= c()
  
  cluster_names = sort(unique(embeddings$leiden_clusters))
  probabilities = list()
  index = c()
  probs = c()
  grps = c()
  for(i in 1:length(cluster_names)){
    grp = cluster_names[i]
    index = c(index,xyMclust$models[[grp]]$observations)
    probs = c(probs, xyMclust$models[[grp]]$uncertainty)
    grps = c(grps,rep(grp,length(xyMclust$models[[grp]]$observations)))
  }
  
  probabilities$rows = index
  probabilities$uncert = probs
  probabilities$clust = grps
  probabilities = as.data.frame(probabilities)
  
  probability = probabilities[match(rownames(embeddings),probabilities$rows),]
  
  if (is.logical(annot)) {
    
    embedding_class <- data.frame(X = embeddings[, 1], Y = embeddings[, 2],
                                  pred = embeddings$leiden_clusters,
                                  uncert = probability$uncert)
  }
  else {
    
    embedding_class <- data.frame(X = embeddings[, 1], Y = embeddings[, 2],
                                  pred = embeddings$leiden_clusters,
                                  uncert = probability$uncert,
                                  truth = annot)
    cat("Confusion matrix")
    print(table(embedding_class$pred, embedding_class$truth))
  }
  
  
  uncertainty = ifelse(embedding_class$uncert <= uncert_thresh,
                       "uncertain", "certain")
  embeddings = data.frame(embedding_class, uncertainty)
  return(embeddings)
}


#'Perform PCA. Used by \code{dim_Red} function.
do_pca_seurat <- function (sce, annot = FALSE, num_threads=1)
{
  cells = ncol(sce)
  genes = nrow(sce)
  approx = TRUE
  if(cells <=100 || genes <=100){
    approx = FALSE
  }
  
  sce <- RunPCA(object = sce, approx = approx)
  return(sce)
  
}

#'Perform UMAP. Used by \code{dim_Red} function.
do_umap_seurat <- function (sce, annot = FALSE,   n.neighbors = 20,
                            num_threads=1,scale=FALSE, seed = 42, min.dist = 0.5,
                            spread = 1, a = NULL, b = NULL,
                            resolution = 0.5, npcs = 50)
{
  print("Here doing UMAP")
  set.seed(seed)
  cells = ncol(sce)
  genes = nrow(sce)
  approx = TRUE
  if(cells <=200 || genes <=200){
    approx = FALSE
  }
  
  sce <- RunPCA(object = sce,features=sce@assays$RNA@var.features,seed.use=seed,
                approx = approx, npcs=npcs)
  
  if(is.null(a) && is.null(b)){
    
    sce <- RunUMAP(object = sce,umap.method = "uwot",
                   seed.use=seed, n.neighbors = n.neighbors,  dims = c(1:min(npcs,nrow(sce))), min.dist = min.dist,
                   spread = spread, metric = "euclidean")
  }
  else{
    sce <- RunUMAP(object = sce,umap.method = "uwot",
                   seed.use=seed, n.neighbors = n.neighbors,  dims = c(1:min(npcs,nrow(sce))),
                   a=a, b=b, metric = "euclidean")
    
  }
  
  sce$orig.ident <- annot
  
  embeddings = list()
  embeddings$X = c(sce@reductions$umap@cell.embeddings[,1])
  embeddings$Y = c(sce@reductions$umap@cell.embeddings[,2])
  
  embeddings = as.data.frame(embeddings)
  rownames(embeddings) = colnames(sce)
  
  nn = FindNeighbors(sce, k.param =n.neighbors)
  clusters = FindClusters(nn, resolution = resolution,
                          method="igraph", algorithm = 4)
  
  embeddings$leiden_clusters = clusters$seurat_clusters
  
  return(embeddings)
}


#'Perform tSNE. Used by \code{dim_Red} function.
do_tSNE_seurat <- function (sce, annot = FALSE, num_threads=1,scale=FALSE, seed = 42 ,
                            npcs = 50, perplex = -1, n.neighbors = 20, resolution = 0.5)
{
  set.seed(seed)
  if(perplex == -1){
    Perp = exp(-0.179 + 0.51*log(ncol(sce)))
  }
  else{
    Perp = perplex
  }
  
  cat(paste("Perplexity set for tSNE is", Perp, "\n"))
  cells = ncol(sce)
  genes = nrow(sce)
  approx = TRUE
  if(cells <=200 || genes <=200){
    approx = FALSE
  }
  
  set.seed(seed)
  
  sce <- RunPCA(object = sce,features = sce@assays$RNA@var.features, npcs = npcs, seed.use=seed, approx=approx)
  set.seed(seed)
  
  sce <- RunTSNE(object = sce,features = sce@assays$RNA@var.features,dims = c(1:min(npcs,nrow(sce))),
                 seed.use=seed,  check_duplicates = FALSE, perplexity = Perp)
  
  sce$orig.ident <- annot
  
  embeddings = list()
  embeddings$X = c(sce@reductions$tsne@cell.embeddings[,1])
  embeddings$Y = c(sce@reductions$tsne@cell.embeddings[,2])
  embeddings = as.data.frame(embeddings)
  rownames(embeddings) = colnames(sce)
  
  embeddings = as.data.frame(embeddings)
  rownames(embeddings) = colnames(sce)
  set.seed(seed)
  
  nn = FindNeighbors(sce, k.param =n.neighbors)
  clusters = FindClusters(nn, resolution = resolution,
                          method="igraph", algorithm = 4)
  
  embeddings$leiden_clusters = clusters$seurat_clusters
  
  return(embeddings)
  
}




#' Dimension reduction
#' @description Performs dimension reduction using UMAP/tSNE/PCA. Then fits GMM on the clusters in latent dimensions.
#' @param sce sce SingleCellExperiment object returned by \code{preprocess_Data} function
#' @param annot Annotation (eg. cell type) of each single cell if available, default is FALSE
#' @param uncert Threshold of uncertainty in modelling / class assignment. Beyong the set threhold, cells will be marked as uncertain. Recommended value is 0.5
#' @param perplex Perpexity parameter for tSNE Set it to -1 if you want the value to be calculated according to the formula described in paper. Recommended value is -1
#' @param reduction Value can be 'umap' or 'tsne' or 'pca'. Set value according to the dimension reduction techniqe you want to use.
#' @param n.neighbors Parameter to create neighborhood graph for UMAP/tSNE
#' @param min.dist min.dist Parameter for UMAP
#' @param spread Parameter for UMAP. Recommended value is 1.5
#' @param a Parameter for UMAP
#' @param b Parameter for UMAP
#' @param resolution Resolution value to create clusters for Leiden clustering. 
#' @param npcs Number of Principal Components to be used for UMAP/tSNE
#' 
#' @return Returns the Embeddings list containing UMAP/tSNE/PCA derived latent dimensions plus class labels obtained using GMM
#' plus marking of cells as certain or uncertain plus supplied annotations along with the probability values
#' @export
#'
#' @import Rtsne
#' @import mclust
#' @import ggplot2
#' @import gplots
#' @examples
#' # Load necessary libraries
#' library(InGene)
#' library(SingleCellExperiment)
#' 
#' # Generate example data
#' set.seed(42)
#' counts <- matrix(rpois(500000, lambda = 10), ncol = 5000, nrow = 5000)
#' genes <- paste0("Gene", 1:nrow(counts))
#' annot <- paste0("CellType", rep(1:2, each = 2500)) 
#' 
#' # Preprocess data
#' sce <- preprocess_Data(
#'   raw_data = counts,
#'   gene_list = genes,
#'   annot = annot
#' )
#' 
#' # Provide cell-wise annotations if you require plots with ground truth
#' # (Already provided above)
#' 
#' # Perform dimension reduction
#' DR_list <- dim_Red(
#'   sce = sce,
#'   annot = annot,
#'   num_threads = 1,
#'   resolution = 0.5,
#'   reduction = 'umap',
#'   verbose = FALSE,
#'   seed = 42,
#'   n.neighbors = 20,
#'   uncert = 0.2,
#'   spread = 1
#' )


dim_Red <- function (sce, annot = FALSE, uncert = 0.2, perplex = 30,
                     num_threads = 1, reduction='umap', n.neighbors = 20, 
                     min.dist = 0.5, spread = 1, verbose = FALSE, seed = 42,
                     a = NULL, b = NULL, resolution = 0.5, npcs=50)
{
  
  
  if(reduction == 'tsne')
  {
    
    
    embeddings = do_tSNE_seurat(sce = sce, annot = annot,
                                num_threads = num_threads,
                                seed = seed, perplex = perplex, n.neighbors = n.neighbors,
                                resolution = resolution, npcs=npcs)
  }
  
  
  if(reduction == 'umap')
  {
    embeddings = do_umap_seurat(sce = sce, annot = annot,
                                min.dist = min.dist,spread = spread,npcs=npcs,
                                num_threads = num_threads, seed = seed, a= a, b = b,
                                n.neighbors = n.neighbors, resolution = resolution)
  }
  
  
  if(reduction == 'pca')
  {
    embeddings = do_pca_seurat(sce = sce, annot = annot, num_threads = num_threads)
    embeddings = list()
    embeddings$X = c(sce@reductions$pca@cell.embeddings[,1])
    embeddings$Y = c(sce@reductions$pca@cell.embeddings[,2])
    embeddings = as.data.frame(embeddings)
    rownames(embeddings) = colnames(sce)
    
    embeddings = as.data.frame(embeddings)
    rownames(embeddings) = colnames(sce)
  }
  
  embeddings = fitGMM(embeddings = embeddings,seed = seed, 
                             annot = annot, uncert_thresh = uncert)
  
  result = list(Embeddings = embeddings)
  
}

#' UMAP/tSNE plot to be visualized with cell annotations
#'
#' @param dim_Red_List The list returned by \code{dim_Red} function
#' @param x The label for X axis. Can be UMAP 1, tSNE 1 or PCA 1 (or any other label according to o user's choice)
#' @param y The label for Y axis. Can be UMAP 2, tSNE 2 or PCA 2 (or any other label according to o user's choice)
#' @return plot
#' @export
#'
#' @import gplots
#' @import ggplot2
#' @import scales

vis_2D_truth <- function(dim_Red_List,x='UMAP 1',y='UMAP 2')
{
  identities <- levels(as.factor(dim_Red_List$Embeddings$truth))
  cc <- hue_pal()(length(identities))
  
  Embeddings = list()
  Embeddings$X = dim_Red_List$Embeddings$X
  Embeddings$Y = dim_Red_List$Embeddings$Y
  
  Embeddings$Labels <- as.factor(dim_Red_List$Embeddings$truth)
  Embeddings = as.data.frame(Embeddings)
  
  
  ggplot(Embeddings, aes(X,Y)) + geom_point(aes(color = Labels),size=0.25)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          axis.title.y = element_text(color = "black", size = 15),
          axis.title.x = element_text(color = "black", size = 15),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +  scale_colour_manual(values = cc)+ scale_shape_manual(values = sample(c(1:20)))+xlab(x)+ylab(y)
}



#' UMAP/tSNE plot to be visualized with cell clusters as decided by Leiden clustering
#' @param dim_Red_List The list returned by \code{dim_Red} function
#' @param x The label for X axis. Can be UMAP 1, tSNE 1 or PCA 1 (or any other label according to o user's choice)
#' @param y The label for Y axis. Can be UMAP 2, tSNE 2 or PCA 2 (or any other label according to o user's choice)
#' @return plot
#' @export
#'
#' @import gplots
#' @import ggplot2
#' @import scales


vis_2D_pred <- function(dim_Red_List,x='UMAP1',y='UMAP2')
{
  identities <- levels(as.factor(dim_Red_List$Embeddings$pred))
  cc <- hue_pal()(length(identities))
  
  
  Embeddings = list()
  Embeddings$X = dim_Red_List$Embeddings$X
  Embeddings$Y = dim_Red_List$Embeddings$Y
  
  Embeddings$Labels <- as.factor(dim_Red_List$Embeddings$pred)
  Embeddings = as.data.frame(Embeddings)
  
  ggplot(Embeddings, aes(X,Y)) + geom_point(aes(color = Labels),size=0.25)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          axis.title.y = element_text(color = "black", size = 15),
          axis.title.x = element_text(color = "black", size = 15),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +  scale_colour_manual(values = cc)+ scale_shape_manual(values = sample(c(1:20)))+xlab(x)+ylab(y)
  
  
}



#'  Subsamples relevant cells from each cluster 
#' @description Feature importance is learned from subsampled cells of each cluster. Finds relevant genes(features) from 2D projections by posing it as a classification problem.
#' @param dim_Red_List The list returned by \code{dim_Red} function
#' @param sce SingleCellExperiment object returned by \code{preprocess_Data} function
#' @param min_Cell No of cells to be sampled from each unsupervised group found by GMM modelling
#' @param choose_random_cells Can be TRUE or FALSE. Set it to true if you need cells to sampled randomly in addition to the cells selected by GMM Recommended value is TRUE.
#' @param min_randomCell_perCluster Number of cells to be selected randomly from each cluster
#' @param n_Trees No of trees to be populated by the Random Forest algorithm
#' @param seed seed value to be set for random sampling and for Random Forest classifier. Default value is 42
#' @param num_threads Number of threads to run Random Forest
#' @return A list containing both gene indices (var_ind), ordered by importance score, and associated important scores (importance). Aslo returns list of cells selected from each cluster (cells_ind).
#'
#' @import ranger
#' @import e1071
#' @export
#' @examples
#' rele_Cells = rele_Cells(dim_Red_List  = DR_list ,
#' sce = sce, min_Cell=20,
#'n_Trees = 1000,num_threads = 4,  
#'choose_random_cells = TRUE, 
#'min_randomCell_perCluster = 20, seed = 42)
rele_Cells <- function(dim_Red_List, sce, min_Cell = 20, n_Trees = 1000,num_threads=1, 
                       choose_random_cells = TRUE, min_randomCell_perCluster = 20, seed = 42)
{
  processed_original_Data = sce[["RNA"]]@data
  cat("You may have to wait for several seconds depending on size of your input data! \n")
  res = list()
  temp_projections <- dim_Red_List$Embeddings
  if (min_Cell == "all") {
    min_Cell = nrow(temp_projections)
  }
  
  colnames(processed_original_Data) <- rownames(temp_projections)
  
  temp_projections <- temp_projections[which(temp_projections$uncertainty == "certain"),]
  clusters_total <- sort(unique(dim_Red_List$Embeddings$pred))
  clusters_projections <- sort(unique(temp_projections$pred))
  uncertain_clusters <- setdiff(clusters_total, clusters_projections)
  uncertain_clus = dim_Red_List$Embeddings[which(dim_Red_List$Embeddings$pred %in%
                                                   uncertain_clusters), ]
  temp_projections_certain = data.frame()
  
  for (i in clusters_total) {
    ind <- which(dim_Red_List$Embeddings$pred == i)
    cells_i = dim_Red_List$Embeddings[ind, ]
    count_cells = min(min_Cell, nrow(cells_i))
    cells_i <- cells_i[order(-cells_i$uncert), ][1:count_cells,
    ]
    temp_projections_certain <- rbind(temp_projections_certain, cells_i)
  }
  
  temp_projections <- rbind(temp_projections, uncertain_clus)
  
  for (i in clusters_projections) {
    ind <- which(dim_Red_List$Embeddings$pred == i)
    count_cells = min(as.integer(min_Cell), length(ind))
    cells_i = dim_Red_List$Embeddings[ind, ]
    cells_i <- cells_i[order(cells_i$uncert), ][1:count_cells,
    ]
    a = setdiff(rownames(cells_i), rownames(temp_projections))
    if (length(a) > 0) {
      ind_rows = match(a, rownames(cells_i))
      cells_i = cells_i[ind_rows, ]
      temp_projections <- rbind(temp_projections, cells_i)
    }
  }
  
  ind_cells = match(rownames(temp_projections), colnames(processed_original_Data))
  
  good_clus <- as.numeric(clusters_total)
  good_clus <- as.data.frame(good_clus)
  U <- good_clus[, 1][order(good_clus[, 1])]
  
  projections_comb_intra = c()
  ind_used = c()
  x <- cbind(X = temp_projections$X, Y = temp_projections$Y)
  y = temp_projections$pred
  x = as.matrix(x)
  rownames(x) = rownames(temp_projections)
  
  dat <- data.frame(x, y)
  limit = as.integer(min_Cell/2)
  annot_selected = c()
  
  #Subsample cells
  for (i in U) {
    
    ind <- which(temp_projections$pred == i)
    cells_i = temp_projections[ind, ]
    ind_others <- which(temp_projections_certain$pred != i)
    data_others = temp_projections_certain[ind_others, ]
    class_others = data_others$pred
    ind_i <- which(temp_projections$pred == i)
    data_i = temp_projections[ind_i, ]
    class_i = data_i$pred
    
    data_combined = rbind(data_others, data_i)
    class_combined = rbind(c(class_others, class_i))
    
    cells_i = temp_projections[ind, ]
    y = as.numeric(class_combined)
    y[y != i] = 0
    y[y == i] = 1
    x = as.matrix(cbind(X = data_combined$X, Y = data_combined$Y))
    dat <- data.frame(x, y)
    svm_indices_i = c()
    count_cells <- min(min_Cell, length(rownames(cells_i)))
    
    #one-vs-rest SVM to get support vectors
    if (i %in% clusters_total) {
      dat$y = as.factor(y)
      set.seed(seed)
      svmfit <- svm(y ~ ., data = dat, kernel = "radial",
                    scale = T, class.weights = "inverse")
      svm_indices <- svmfit$index
      svm_indices_i <- intersect(svm_indices, rownames(cells_i))
      ind_svm = match(svm_indices_i, rownames(temp_projections))
      cells_svm_i <- temp_projections[ind_svm, ]
      count_svm_cells = min(limit, length(rownames(cells_svm_i)))
      cells_svm_i <- cells_svm_i[order(-cells_svm_i$uncert),
      ]
      s_svm = rownames(cells_svm_i)[1:count_svm_cells]
      s_svm_ascend = tail(rownames(cells_svm_i), count_svm_cells)
    }
    
    cells_i <- cells_i[order(-cells_i$uncert), ]
    
    s = rownames(cells_i)[1:count_cells]
    s_ascend = tail(rownames(cells_i), min(limit, count_cells))
    
    if (length(svm_indices_i) == 0) {
      s_svm = s
      s_svm_ascend = s_ascend
    }
    all_s <- c(unique(c(s, s_svm, s_ascend, s_svm_ascend)))
    matching_ind <- match(all_s, rownames(temp_projections))
    ind_used <- c(ind_used, all_s)
    annot_selected = c(annot_selected, all_s)
    if (length(matching_ind) >= 2) {
      projections_comb_intra <- rbind(projections_comb_intra, t(combn(matching_ind,
                                                                      2)))
    }
    else {
      projections_comb_intra <- rbind(projections_comb_intra, t(combn(c(matching_ind,
                                                                        matching_ind), 2)))
    }
  }
  
  
  projections_comb_inter <- c()
  matching_ind <- match(ind_used, rownames(temp_projections))
  all_comb <- data.frame(t(combn((matching_ind), 2)))
  
  
  for (i in 1:nrow(all_comb)) {
    if (is.null(projections_comb_inter) == TRUE || nrow(projections_comb_inter) <
        nrow(projections_comb_intra)) {
      if (temp_projections[all_comb[i, 1], 3] != temp_projections[all_comb[i,2], 3]) {
        projections_comb_inter <- rbind(projections_comb_inter, all_comb[i,])
      }
    }
    else {
      break
    }
  }
  
  #Train Random Forest classifier to find feature importance
  classification_data = c()
  class = c()
  check_data <- processed_original_Data[, ind_cells]
  
  mat_ind_1 = check_data[, projections_comb_inter[, 2]]
  mat_ind_2 = check_data[, projections_comb_inter[, 1]]
  
  numerator = (abs(t(mat_ind_1) - t(mat_ind_2)))
  denominator = ((t(mat_ind_1) + t(mat_ind_2)))
  numerator <- numerator + 1
  denominator <- denominator + 1
  
  Man_mat1 = ((numerator/denominator))
  mat_ind_1 = as.matrix(check_data[, projections_comb_intra[, 2]])
  mat_ind_2 = as.matrix(check_data[, projections_comb_intra[, 1]])
  
  
  numerator = (abs(t(mat_ind_1) - t(mat_ind_2)))
  denominator = ((t(mat_ind_1) + t(mat_ind_2)))
  numerator <- numerator + 1
  denominator <- denominator + 1
  Man_mat2 = ((numerator/denominator))
  
  
  classification_data = rbind(Man_mat1, Man_mat2)
  class = c(rep(0, nrow(projections_comb_inter)), rep(1, nrow(projections_comb_intra)))
  MAT = as.matrix(cbind(classification_data, class))
  rownames(MAT) = NULL
  MAT <- data.frame(MAT)
  MAT$class <- as.factor(MAT$class)
  
  rf_model <- ranger(dependent.variable.name = "class", data = MAT,
                     write.forest = TRUE, importance = "impurity", num.trees = n_Trees,
                     save.memory = T, num.threads = num_threads, seed = seed)
  
  ind = order(rf_model$variable.importance, decreasing = T)
  imp = rf_model$variable.importance[ind]
  
  sample_cells = rownames(dim_Red_List$Embeddings)
  
  if(choose_random_cells == TRUE){
    
    set.seed(seed)
    random_cells_data = c()
    clus = sort(unique( dim_Red_List$Embeddings$pred))
    for(i in clus){
      ind_clust <- which(dim_Red_List$Embeddings$pred == i) #indices of all data points that have the label of the U[i]th cluster.
      size = min(min_randomCell_perCluster, length(ind_clust))
      set.seed(seed)
      ind_i = sample(ind_clust, size)
      random_cells_data = c(random_cells_data,ind_i)
    }
    
    sample_cells = unique(c(ind_used, random_cells_data))
  }
  
  res = list(var_ind = ind, importance = imp, cells_ind = ind_used, sample_cells = sample_cells)
  print(rf_model)
  return(res)
  
}




#' Calculate Construction Error of 2D Plot (tSNE/UMAP/PCA)
#' Used by \code{getHighestCorrelation} function
constructionError <- function(original_data, est_dist)
{
  
  
  # original distance matrix
  DR_new = t(original_data)
  ori_dist = as.matrix(dist((DR_new)))
  return(cor(as.vector(ori_dist),as.vector(est_dist),method="spearman"))
}


#' Finds correlation between original projection and latent dimensions
#' @description  Calculates Spearman correlation coefficient of genes in latent dimension with the original projection based on euclidean distance.
#' @import Rtsne
#' @import mclust
#' @param sce sce SingleCellExperiment object returned by \code{preprocess_Data} function
#' @param cl_res List object returned by \code{rele_Cells} function
#' @param reduced List object of latent dimensions returned by \code{dim_Red} function
#' @param n_gene_seq Sequence of number of genes from 10 to desired number (default 1000), by a periodic gap of 10
#' @param sample_cells Boolean value parameter. TRUE indicates that sampled cells (by GMM and random selection) will be used to calculate correlation. FALSE indicates all cells will be used. Recommended value is TRUE.
#' 
#' @return Returns the list with correlation value for the genes. Considering 'n' genes with the highest correlation value will give us the set of candidate genes that explain the latent projection in the best way
#' @export
getHighestCorrelation <- function(sce, reduced, cl_res, n_gene_seq = seq(10,1000,by = 10), sample_cells=TRUE)
{
  data = sce[["RNA"]]@data
  ori = c()
  imp_var_ind = cl_res$var_ind
  reduced = reduced$Embeddings[,c(1:2)]
  
  if(sample_cells == TRUE){
    
    reduced = reduced[cl_res$sample_cells,]
    data = data[,cl_res$sample_cells]
  }
  
  est_dist = as.matrix(dist((reduced)))
  
  for(i in 1:length(n_gene_seq))
  {
    
    rows = rownames(data)[imp_var_ind[1:n_gene_seq[i]]]
    
    
    error = constructionError(original_data = data[rows,],est_dist=est_dist )
    
    print(error)
    ori = append(ori,error)
  }
  
  print(paste("Highest correlation is ",ori[which(ori==max(ori))] ,"obtained with top",n_gene_seq[which(ori==max(ori))],"genes"))
  
  return(ori)
}

#' Finds top N relevant genes of the latent dimension
#' @description Uses \code{getHighestCorrelation} to find top N genes influencing the UMAP/tSNE
#' @param cl_res List object returned by \code{rele_Cells} function
#' @param nTop Number of top relevant genes to be considered. Set it to 0 if you want to find it using using maximum correlation. Recommended value is 500.
#' @param n_gene_seq Sequence of number of genes from 10 to desired number (default 1000), by a periodic gap of 10
#' @param ori List with correlation values returned by \code{getHighestCorrelation} function
#' @return Returns list with top relevant genes
#' @export

getRelevantGenes <- function(cl_res, nTop = 500, 
                             ori=c(), n_gene_seq = seq(10,1000,by = 10), genes)
{
  if((is.na(nTop) | nTop == 0) & (length(ori) == 0)){
    stop("Need either top number of genes to select or list with correlation values")
  }
  
  if(nTop == 0 | is.na(nTop)){
    nTop = n_gene_seq[which(ori==max(ori))]
    
  }
  
  if(length(ori) > length(genes)){
    stop("List of correlation values cannot be greater than number of genes")
  }
  
  topk_genes = cl_res$var_ind[1:nTop]
  topGenes = genes[topk_genes]
  
  return(topGenes)
  
}

