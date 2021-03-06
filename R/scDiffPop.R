
setClass("scDiffPop", slots = list(results = "data.frame",
                                   tree = "matrix",
                                   counts = "matrix",
                                   meta.data = "list",
                                   markers = "list",
                                   pathways = "list"))

scDiffPop <- function(sco, annotation = sco@meta.data$cellType,
                      patient_id = sco@meta.data$patient,
                      response = sco@meta.data$response,
                      nmarkers = 25, use_seurat_clusters = FALSE,
                      find_pathways = FALSE,
                      nperm = 250,
                      nmarker = 25,
                      ncores = 1) {

  ## If using seurat clusters, then the cell types are the cluster IDs.
  if(use_seurat_clusters) {
    #CHECK THIS EXISTS
    cell_types <- as.integer(sco$seurat_clusters)
    sco@meta.data$cellType <- sapply(cell_types, function(x) {paste("s", x, sep ="")})
  }
  else {
    sco@meta.data$cellType <- annotation
  }

  sco@meta.data$patient <- patient_id
  sco@meta.data$response <- response


  ### Make the cell Tree
  Tree <- makeCellTree(sco)

  cell_types <- unique(sco@meta.data$cellType)
  phenotypes <- sort(unique(sco@meta.data$response))
  results <- data.frame(group = Tree[,1], enrichment = rep(0, nrow(Tree)), effect = rep(0, nrow(Tree)), pval = rep(0, nrow(Tree)),
                        padj = rep(0, nrow(Tree)), hetero = rep(0, nrow(Tree)),
                        lm_coef=rep(0,nrow(Tree)), lm_pval=rep(0,nrow(Tree)),
                        rpca_coef=rep(0,nrow(Tree)),rpca_pval=rep(0,nrow=Tree))
  marker_list <- list()
  counts <- matrix(0, nrow = nrow(Tree)+1, ncol = 2); colnames(counts) <- c(phenotypes[1], phenotypes[2])
  counts[1,] <- c(length(which(sco@meta.data$response == phenotypes[1])), length(which(sco@meta.data$response == phenotypes[2])))
  for(i in 1:nrow(Tree)) {
    cat("ITERATION: ", i, "\n")
    subtree <- as.vector(DFS(Tree, i, cell_types))
    print(subtree)
    counts[i+1,1] <- sum(sco@meta.data$cellType[sco@meta.data$response == phenotypes[1]] %in% subtree)
    counts[i+1,2] <- sum(sco@meta.data$cellType[sco@meta.data$response == phenotypes[2]] %in% subtree)

    #Find Markers of this node
    Idents(sco) = as.factor(ifelse(sco$cellType %in% subtree, 1 ,2))
    markers <- Seurat::FindMarkers(sco, min.pct = 0.1, only.pos = TRUE, logfc.threshold = 0.25, ident.1 = 1)
    markers <- markers[1:min(nmarkers, nrow(markers)),]

    markers$p_val_adj[markers$p_val_adj == 0] <- .Machine$double.xmin
    markers$l2FC
    x <- rep(1, nmarkers); names(x) <- rownames(markers)
    dds <- getPseudoBulkCounts(sco, subtree, design)
    res <- results(dds)
    print("DESEQ Complete")
    y <- res$stat; names(y) <- rownames(res); y <- y[names(y) %in% names(x)]; y <- y[names(x)]
    dds <- dds[rownames(dds) %in% names(x),]

    genes_use <- list()
    genes_use$main <- Tree[i,1]
    genes_use$x <- x; genes_use$y <- y
    marker_list[[i]] <- genes_use

    print(x); print(y)
    stat <- sum(x*y)
    print("STAT:"); print(stat/nmarkers)
    pval <- permutation_test(x, nperm, dds, stat, ncores)
    results$effect[i] <- stat/nmarkers
    results$pval[i] <- pval

    print("Running mixed effects model")
    y <- ifelse(sco@meta.data$cellType %in% subtree, 1, 0)
    lm.data <- data.frame(y=y,response=sco@meta.data$response,patient=sco@meta.data$patient)
    fit <- glmer(y~response + (1|patient),family="binomial",data=lm.data)
    print(summary(fit))
    results$lm_coef[i] <- summary(fit)$coefficients[2,1]
    results$lm_pval[i] <- summary(fit)$coefficients[2,4]

    ## PCA MODEL
    embeddings <- sco@reductions$pca@cell.embeddings[sco@meta.data$cellType %in% subtree,]
    npcs <- 50
    coefs <- rep(0, npcs)
    tvals <- rep(1,npcs)
    for(j in 1:npcs) {
      print(j)
      lm.data <- data.frame(y=embeddings[,j],response=sco@meta.data$response[sco@meta.data$cellType %in% subtree],
                            patient=sco@meta.data$patient[sco@meta.data$cellType %in% subtree])
      try({fit <- lmer(y ~ response + (1|patient),data=lm.data)
      sumfit <- summary(fit)
      coefs[j] <- sumfit$coefficients[2,1]
      tvals[j] <- sumfit$coefficients[2,5]})
    }
    results$rpca_coef[i] <- sqrt(sum(coefs^2))
    results$rpca_pval[i] <- min(tvals)

    ifelse(stat > 0, results$enrichment[i] <- phenotypes[2], results$enrichment[i] <- phenotypes[1])
  }
  results$padj <- p.adjust(results$pval, method = "fdr")
  meta.data <- list()
  meta.data$phenotypes <- c(phenotypes[2], phenotypes[1])
  out <- new("scDiffPop", results = results, tree = Tree, meta.data = meta.data, markers = marker_list, counts = counts)
  return(out)
}

makeCellTree <- function(sco) {
  cell_types <- unique(sco@meta.data$cellType)
  group <- rep(1, length(cell_types))
  TreeMat <- matrix(1, nrow = length(cell_types), ncol = 1)
  counter <- 1

  while(length(unique(group)) != length(cell_types)) {
    #Get group with most clusters
    ixs <- which(group == Mode(group))
    print(group)
    cat("Current group: ", cell_types[ixs], "\n")
    split <- splitGroup(sco[,sco$cellType %in% cell_types[ixs]], cell_types[ixs])

    if(length(unique(split)) > 1) {
      counter <- counter + 1
      subgroup <- group[ixs]
      subgroup[split == 2] <- counter
      counter <- counter + 1
      subgroup[split == 1] <- counter
      group[ixs] <- subgroup

      TreeMat <- cbind(TreeMat, group)
    }
    else {
      for(j in ixs) {
        counter <- counter + 1
        group[j] <- counter
        TreeMat <- cbind(TreeMat, group)
      }
    }
  }

  Tree <- matrix(nrow = 0, ncol = 2)

  counter <- 1

  for(i in 2:ncol(TreeMat)) {
    newvals <- unique(TreeMat[TreeMat[,i] > counter,i])
    for(j in newvals) {
      counter <- counter + 1
      ixs <- which(TreeMat[,i] == counter)
      newvec <- c(counter, TreeMat[ixs[1], i-1])
      Tree <- rbind(Tree, newvec)
    }
  }

  ### Process Tree

  #Change names of leaves
  cntr <- 1
  for(ct in TreeMat[,ncol(TreeMat)]) {
    ix <- which(Tree[,1] == ct)
    Tree[ix,1] <- cell_types[cntr]
    cntr <- cntr+1
  }

  # Also give the columns a name
  colnames(Tree) <- c("Child", "Parent")
  rownames(Tree) <- c(1:nrow(Tree))

  return(Tree)
}

splitGroup <- function(sco_sub, ixs) {
  downsample <- min(100, nrow(sco_sub@meta.data))
  down_ixs <- c()
  for(i in ixs) {
    rxs <- which(sco_sub@meta.data$cellType == i)
    downsample <- min(100, length(rxs))
    rxs <- rxs[sample(length(rxs), downsample, replace = FALSE)]
    down_ixs <- c(down_ixs, rxs)
  }
  sco_sub <- sco_sub[,down_ixs]
  if(length(ixs) == 2) {
    return(c(1,2))
  }

  sco_sub <- Seurat::NormalizeData(sco_sub)

  #Find variable features
  sco_sub <- Seurat::FindVariableFeatures(sco_sub, verbose = FALSE)

  #Scale data
  sco_sub <- Seurat::ScaleData(sco_sub, verbose = FALSE)

  #Run PCA on the variable features. Get 50 dimensional embeddings
  sco_sub <- Seurat::RunPCA(sco_sub, verbose = FALSE)
  embeddings <- Seurat::Embeddings(object = sco_sub, reduction = 'pca')[,1:50]

  pseudobulk <- matrix(0, nrow = 0, ncol = 50)
  for(i in ixs) {
    rxs <- which(sco_sub@meta.data$cellType == i)
    pseudobulk <- rbind(pseudobulk, colMeans(embeddings[rxs,]))
  }
  #Cluster via k means
  km <- kmeans(pseudobulk, centers = 2, iter.max = 10)$cluster

  return(km)
}

Mode <- function(x) {
  xu <- unique(x)
  return(xu[which.max(tabulate(match(x, xu)))])
}

DFS <- function(Tree, node, cell_types) {
  if(node == 0) {
    return(cell_types)
  }

  leaves <- c()

  row <- Tree[node,]

  newparent <- row[1]

  allchild <- which(Tree[,2] == row[1])

  for(child in allchild) {
    leaves <- c(leaves, DFS(Tree, child))
  }

  if(length(allchild) == 0) {
    return(row[1])
  }

  return(leaves)
}

getPseudoBulkCounts <- function(sco, subtree, response) {
  counts <- sco@assays$RNA@counts
  pseudobulk <- matrix(0, nrow = nrow(counts), ncol = length(unique(sco@meta.data$patient)))
  rownames(pseudobulk) <- rownames(counts)
  response <- rep(0, ncol(pseudobulk))
  pseudoMetaData <- matrix(0, nrow = ncol(pseudobulk), ncol = ncol(sco@meta.data))
  pseudoMetaData <- as.data.frame(pseudoMetaData); colnames(pseudoMetaData) <- colnames(sco@meta.data)
  for(i in 1:ncol(pseudobulk)) {
    ixs <- which(sco@meta.data$patient == unique(sco@meta.data$patient)[i])
    response[i] <- sco@meta.data$response[ixs[1]]
    pseudoMetaData[i,] <- sco@meta.data[ixs[1],]
  }
  sco_sub <- sco[,sco$cellType %in% subtree]
  counts <- sco_sub@assays$RNA@counts
  for(i in 1:ncol(pseudobulk)) {
    ixs <- which(sco_sub@meta.data$patient == unique(sco@meta.data$patient)[i])
    if(length(ixs) == 0) {
      pseudobulk[,i] <- 1
    }
    else if(length(ixs) == 1) {
      pseudobulk[,i] <- counts[,ixs] + 1
    }
    else {
      pseudobulk[,i] <- Matrix::rowSums(counts[,ixs]) + 1
    }
  }

  colData <- data.frame(gene = c(1:length(response)), response = as.factor(response))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = pseudobulk, colData = colData, design = ~response)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  return(dds)
}

runDESeq <- function(counts, colData, response) {
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~response)
  dds <- DESeq2::DESeq(dds)
  results <- DESeq2::results(dds)
  return(results)
}

permutation_test <- function(x, nperm, dds, stat, ncores) {
  vec <- c(1:ncol(dds))
  print(dds)
  null_dist <- unlist(parallel::mclapply(c(1:nperm), function(i) {
    vec <- sample(vec, size = length(vec), replace = FALSE)
    dds$response <- dds$response[vec]
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    y <- res$stat; names(y) <- rownames(res); y <- y[names(y) %in% names(x)]; y <- y[names(x)]
    return(sum(x*y))
  }, mc.cores = ncores))
  print(null_dist)
  pval <- (sum(abs(null_dist) > abs(stat)) + 1)/(nperm + 1)
  print("PVAL:"); print(pval)
  return(pval)
}
