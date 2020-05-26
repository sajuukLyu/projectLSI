
#' Calculate Latent Semantic Indexing (LSI)
#'
#' Perform svd and then LSI to the input count matrix.
#'
#' @param mat the normalized sparse count matrix
#' @param ndims how many dimensions will be calculated, default 50
#' @param seed.use the random seed used, default 42
#'
#' @return a list containing the result of svd and parameters used
#'
#' @import Matrix
#' @importFrom irlba irlba
#'
#' @export
calcLSI <- function(mat, ndims = 50, seed.use = 42){

  set.seed(seed.use)

  # filter genes
  mat <- mat[rowSums(mat) > 0, ]

  # tf-idf transform
  colSm <- colSums(mat)
  rowSm <- rowSums(mat)

  tf <- t(t(mat) / colSm)
  idf <- log(1 + ncol(mat) / rowSm)

  tfidf <- as(Diagonal(x = idf), "sparseMatrix") %*% tf

  # calc LSI
  svd <- irlba(tfidf, ndims)

  svdDiag <- diag(svd$d)

  matSVD <- svd$v %*% svdDiag
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC_", seq_len(ncol(matSVD)))

  fLoad <- svd$u
  rownames(fLoad) <- rownames(mat)
  colnames(fLoad) <- colnames(matSVD)

  sdev <- svd$d / sqrt(max(1, ncol(mat) - 1))

  out <- list(
    matSVD = matSVD,
    fLoad = fLoad,
    rowSm = rowSm,
    colSm = colSm,
    sdev = sdev,
    svd = svd,
    ndims = ndims,
    seed.use = seed.use)

  out
}


#' Project to Given Latent Semantic Indexing (LSI)
#'
#' Project the input count matrix to a given LSI space using pre-calculated svd output.
#'
#' @param mat the normalized sparce count matrix
#' @param lsi a list of the output of calcLSI()
#'
#' @return a dense matrix of svd cell loading
#'
#' @import Matrix
#'
#' @export
projectLSI <- function(mat, lsi){

  # filter genes
  ref.gene <- rownames(lsi$fLoad)
  used.gene <- intersect(rownames(mat), ref.gene)
  zero.gene <- setdiff(ref.gene, used.gene)

  used.mat <- rbind(mat[used.gene, ],
                    matrix(0, nrow = length(zero.gene), ncol = ncol(mat),
                           dimnames = list(zero.gene, colnames(mat))))
  used.mat <- used.mat[ref.gene, ]

  # tf-idf transform
  tf <- t(t(used.mat) / colSums(used.mat))
  idf <- log(1 + length(lsi$colSm) / lsi$rowSm)

  tfidf <- as(Diagonal(x = idf), "sparseMatrix") %*% tf

  # projection to given svd
  matSVD <- t(tfidf) %*% lsi$svd$u
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC_", seq_len(ncol(matSVD)))

  as.matrix(matSVD)
}

#' Calculate Psudo-single-cell Count Data
#'
#' Calculate the psudo-single-cell RNA-seq data by down-sampling bulk RNA-seq counts.
#'
#' @param bulk.data the bulk RNA-seq counts matrix, genes in row and samples in column.
#' @param n how many times to down-sample per bulk sample, default 100
#' @param depth how many counts down-sample to, default 3000
#'
#' @return a sparse matrix of down-sampled count matrix with genes in row and n*sample psudo-cells in column
#'
#' @import Matrix
#'
#' @export
psudoSC <- function(bulk.data, n = 100, depth = 3000){

  bulk.data <- ceiling(bulk.data)
  psudo.bulk <- list()

  message("downsampling counts...")
  for(i in colnames(bulk.data)){

    all.count <- rep(rownames(bulk.data), bulk.data[, i])

    psudo.matrix <- matrix(0, nrow = nrow(bulk.data), ncol = n,
                           dimnames = list(rownames(bulk.data), paste0(i, "_", 1:n)))

    for(j in 1:n){
      psudo.count <- table(sample(all.count, depth))
      psudo.matrix[names(psudo.count), j] <- psudo.count
    }

    psudo.bulk[[i]] <- as(psudo.matrix, "sparseMatrix")
  }

  message("merging all samples...")
  psudo.all <- psudo.bulk[[1]]
  for(i in names(psudo.bulk)[2:ncol(bulk.data)]){
    psudo.all <- cbind(psudo.all, psudo.bulk[[i]])
  }

  psudo.all
}


