#' Estimate the Number of Principal Components (Internal)
#'
#' Determines the optimal number of principal components (PCs) to retain
#' from a data matrix using either the "elbow" or "permutation" method.
#'
#' @param data A numeric matrix, list, or SVD/rSVD object. If a matrix is provided,
#'   the function performs singular value decomposition (SVD). If a list or `rsvd`
#'   object with a component `d` is provided, its singular values are used directly.
#' @param method Character string specifying the method used to estimate
#'   the number of PCs. One of:
#'   \itemize{
#'     \item `"elbow"` — uses the curvature ("elbow") in the singular value spectrum.
#'     \item `"permutation"` — compares the observed singular values to those
#'       from permuted data to estimate significance.
#'   }
#'   Defaults to `"elbow"`.
#' @param B Integer specifying the number of permutations to use for the
#'   `"permutation"` method. Ignored if `method = "elbow"`. Default is 20.
#' @param seed Optional integer for reproducibility. Default is `NULL`.
#'
#' @return An integer giving the estimated number of principal components (`nsv`).
#'
#' @details
#' The function provides two strategies for estimating the number of significant
#' components in a dataset:
#' \enumerate{
#'   \item **Elbow method** — identifies the point of diminishing returns in the
#'     singular value spectrum by analyzing the second derivative of the singular values.
#'   \item **Permutation method** — compares observed singular values against
#'     a null distribution obtained by permuting the data `B` times. Components
#'     whose variance exceeds the null expectation (p ≤ 0.1) are considered significant.
#' }
#'
#' For large datasets, the function uses randomized SVD (`rsvd`) for computational efficiency.
#' If the input is already an `rsvd` or SVD object, the singular values are used directly.
#'
#' @examples
#' \dontrun{
#' mat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' num.pc(mat, method = "elbow")
#' num.pc(mat, method = "permutation", B = 50, seed = 123)
#' }
#'
#' @seealso \code{\link{svd}}, \code{\link{rsvd}}
#'
#' @keywords internal


num.pc = function (data, method="elbow", B = 20, seed = NULL)
{

  method=match.arg(method, c("elbow", "permutation"))
  if (!is.null(seed)) {
    set.seed(seed)
  }
  warn <- NULL
  if((class(data)!="list")&(class(data)!="rsvd")){
    message("Computing svd")
    n <- ncol(data)
    m <- nrow(data)
    data=rowNorm(data)
    if(n<500){
      k=n
    }
    else{
      k=max(200,n/4)
    }
    if(k==n){
      uu <- svd(data)
    }
    else{
      set.seed(123456);uu <- rsvd(data,k, q=3)
    }
  }
  else if (!is.null(data[["d"]])){
    if(method=="permutation"){
      message("Original data is needed for permutation method.\nSetting method to elbow")
      method="elbow"
    }

    uu=data

  }



  if(method=="permutation"){
    #nn = min(c(n, m))
    dstat <- uu$d[1:k]^2/sum(uu$d[1:k]^2)
    dstat0 <- matrix(0, nrow = B, ncol = k)
    for (i in 1:B) {
      dat0 <- t(apply(data, 1, sample, replace = FALSE))
      if(k==n){
        uu0 <- svd(dat0)
      }
      else{
        set.seed(123456);
        uu0 <- rsvd(dat0,k, q=3)
      }
      dstat0[i, ] <- uu0$d[1:k]^2/sum(uu0$d[1:k]^2)
    }
    psv <- rep(1, k)
    for (i in 1:k) {
      psv[i] <- mean(dstat0[, i] >= dstat[i])
    }
    for (i in 2:k) {
      psv[i] <- max(psv[(i - 1)], psv[i])
    }

    nsv <- sum(psv <= 0.1)
  }
  else if (method=="elbow"){
    x=smooth(xraw<-abs(diff(diff(uu$d))), twiceit = T)
    #plot(x)


    nsv=which(x<=quantile(x, 0.5))[1]+1

  }
  return(nsv)
}



#
#' Compute Ridge-Regularized Pseudoinverse
#'
#' Computes the Moore–Penrose pseudoinverse of a matrix with optional
#' ridge (Tikhonov) regularization for numerical stability.
#'
#' @param m A numeric matrix to invert.
#' @param alpha A non-negative numeric scalar specifying the ridge
#'   regularization strength. Use `alpha = 0` for the standard pseudoinverse.
#'
#' @return A matrix representing the (regularized) pseudoinverse of `m`.
#'

#' @keywords internal

pinv.ridge=function (m, alpha = 0)
{
  msvd = svd(m)
  if (length(msvd$d) == 0) {
    return(array(0, dim(m)[2:1]))
  }
  else {
    if (alpha > 0) {
      ss = (msvd$d^2) + alpha^2
      msvd$d = ss/msvd$d
    }
    out = msvd$v %*% (1/msvd$d * t(msvd$u))
    rownames(out) = rownames(m)
    colnames(out) = colnames(m)
    out
  }
}


######
#' generate a binary gene-pathway matrix from a gmt file
#'
#'
#'
#'
#' @param gmt file; the gmt file name
#'
#'
#' @return A matrix representing the gene-pathway, row is gene, the column is pathway.
#'

#' @export

gmt_to_binary_matrix <- function(gmt_file) {
  # Step 1: Read lines
  gmt_lines <- readLines(gmt_file)

  # Step 2: Parse each line into list of vectors
  gmt_list <- lapply(gmt_lines, function(line) strsplit(line, "\t")[[1]])

  # Step 3: Create named list of gene sets
  gmt_named <- setNames(
    lapply(gmt_list, function(x) x[-c(1, 2)]),  # remove name and description
    sapply(gmt_list, function(x) x[1])          # pathway names
  )

  # Step 4: Get unique gene list
  all_genes <- unique(unlist(gmt_named))
  pathway_names <- names(gmt_named)

  # Step 5: Create and fill binary matrix
  gene_by_pathway <- matrix(0, nrow = length(all_genes), ncol = length(pathway_names),
                            dimnames = list(all_genes, pathway_names))

  for (path in pathway_names) {
    gene_set <- gmt_named[[path]]
    gene_by_pathway[gene_set, path] <- 1
  }

  return(gene_by_pathway)
}




#####
#' Compute pathway-LF matrix based on gene-LF and gene-pathway matrices
#'
#'
#'
#' @param Z matrix; gene-LF matrix with gene in rows and LF in columns.
#' @param Chat matrix; regularized inverse of C (gene-pathway matrix), which is used to derive a prior U and select candidate pathways
#' @param pathMat; a binary gene-pathway matrix with gene in rows and pathways in columns
#' @param glm_alpha; alpha parameter in glmnet
#' @param maxPath; the maximum number of pathway to be considered
#' @param target.frac; the proportion of latent factors with pathway enriched
#' @param L3; the lambada parameter to control the sparsity of pathway
#' @return U matrix; pathway-LF matrix with pathway in rows and LF in columns
#'

#' @export
solveU=function(Z,  Chat, pathMat,  glm_alpha=0.9, maxPath=10, target.frac=0.7, L3=NULL){


  Ur=Chat%*%Z #get U by OLS
  Ur=apply(-Ur,2,rank) #rank
  Urm=apply(Ur,1,min)

  U=matrix(0,nrow=ncol(pathMat), ncol=ncol(Z))
  if(is.null(L3)){

    lambdas=exp(seq(-4,-12,-0.125))
    results=list()
    lMat=matrix(nrow=length(lambdas), ncol=ncol(Z))
    for(i in 1:ncol(Z)){

      iip=which(Ur[,i]<=maxPath)
      gres=glmnet::glmnet(y=Z[,i], x=pathMat[,iip], alpha=glm_alpha, lower.limits=0, lambda = lambdas,intercept=T,  standardize=F )

      gres$iip=iip
      lMat[,i]=colSums(as.matrix(gres$beta)>0)
      results[[i]]=gres
    }
    fracs=rowMeans(lMat>0)
    iibest=which.min(abs(target.frac-fracs))
    iibest


    for(i in 1:ncol(Z)){
      U[results[[i]]$iip,i]=results[[i]]$beta[,iibest]
    }#for i
    rownames(U)=colnames(pathMat)
    colnames(U)=1:ncol(Z)

    Utmp=solveU(Z,  Chat, pathMat, glm_alpha=0.9, maxPath=10,  L3=lambdas[iibest])

    #stop()
    return(list(U=U, L3=lambdas[iibest]))
  }
  else{ #do one fit with a given lambda
    for(i in 1:ncol(Z)){

      iip=which(Ur[,i]<=maxPath)

      gres=glmnet::glmnet(y=Z[,i], x=pathMat[,iip], alpha=glm_alpha, lower.limits=0, lambda = L3,intercept=T,  standardize=F )
      U[iip,i]=as.numeric(gres$beta)
    }

    return(U)
  }
}

######


#' Copy a matrix
#'
#' @keywords internal
copyMat <- function(mat, zero = FALSE) {
  # Create new matrix with same dimensions
  matnew <- matrix(nrow = nrow(mat), ncol = ncol(mat))

  # Copy row and column names
  rownames(matnew) <- rownames(mat)
  colnames(matnew) <- colnames(mat)

  # Initialize with zeros if requested
  if (zero) {
    matnew[] <- 0
  }

  return(matnew)
}




#####
#' Cross validation of results from mcCoordPath


#' @keywords internal

crossVal=function(coPathRes, data, pathMat, pathMatcv){

  out=matrix(ncol=4, nrow=0)
  ii=which(colSums(coPathRes$U)>0)
  Uauc=copyMat(coPathRes$U,T)
  Up=copyMat(coPathRes$U,T)
  Up[]=1
  for ( i in ii){

    iipath=which(coPathRes$U[,i]>0)

    if (length(iipath) > 1){
      for(j in iipath){
        iiheldout=which((rowSums(pathMat[,iipath, drop=F])==0) |(pathMat[,j]>0&pathMatcv[,j]==0))
        aucres=AUC(pathMat[iiheldout,j], coPathRes$Z[iiheldout,i])
        out=rbind(out,c(colnames(pathMat)[j], i, aucres$auc, aucres$pval))
        Uauc[j,i]=aucres$auc
        Up[j,i]=aucres$pval
      }}else{
        j <- iipath
        iiheldout=which((rowSums(matrix(pathMat[,iipath],ncol=1))==0) |(pathMat[,j]>0&pathMatcv[,j]==0))
        aucres=AUC(pathMat[iiheldout,j], coPathRes$Z[iiheldout,i])
        out=rbind(out,c(colnames(pathMat)[j], i, aucres$auc, aucres$pval))
        Uauc[j,i]=aucres$auc
        Up[j,i]=aucres$pval
      }#else
  }
  out=data.frame(out,stringsAsFactors = F)
  out[,3]=as.numeric(out[,3])
  out[,4]=as.numeric(out[,4])
  out[,5]=p.adjust(out[,4],method="BH")
  colnames(out)=c("pathway", "LF index", "AUC", "p-value", "FDR")
  return(list(Uauc=Uauc, Upval=Up, summary=out))
}


#####
#' Compute AUC


#' @keywords internal
#'
AUC<-function(labels, values){
  posii=which(labels>0)
  negii=which(labels<=0)
  posn=length(posii)
  negn=length(negii)
  posval=values[posii]
  negval=values[negii]
  myres=list()
  if(posn>0&negn>0){
    testres=wilcox.test(posval, negval, alternative="greater", conf.int=TRUE);

    myres$low=testres$conf.int[1]
    myres$high=testres$conf.int[2]
    myres$auc=(testres$statistic)/(posn*negn)
    myres$pval=testres$p.value
  }
  else{
    myres$auc=0.5
    myres$pval=NA
  }
  return(myres)
}


#####
#' Calculate AUC of results from mcCoordPath


#' @keywords internal
#'
getAUC=function(coPathRes, data, pathMat){
  Y=as.matrix(data)
  B=coPathRes$B
  Z=coPathRes$Z
  Zcv=copyMat(Z)
  k=ncol(Z)
  L1=coPathRes$L1
  L2=coPathRes$L2
  for (i in 1:5){
    ii=(0:(floor(nrow(data)/5)-1))*5+i
    ii=ii[ii<=nrow(Z)]


    Bcv=solve(crossprod(Z[-ii,])+L2*diag(k))%*%t(Z[-ii,])%*%Y[-ii,]

    Zcv[ii,]=Y[ii, ]%*%t(Bcv)%*%solve(tcrossprod(Bcv)+L1*diag(k))
  }

  out=matrix(ncol=4, nrow=0)
  ii=which(colSums(coPathRes$U)>0)
  Uauc=copyMat(coPathRes$U,T)
  Up=copyMat(coPathRes$U,T)
  Up[]=1;
  for ( i in ii){

    iipath=which(coPathRes$U[,i]>0)

    for(j in iipath){
      aucres=AUC(pathMat[,j], Zcv[,i])
      out=rbind(out,c(colnames(pathMat)[j], i, aucres$auc, aucres$pval))
      Uauc[j,i]=aucres$auc
      Up[j,i]=aucres$pval
    }
  }
  out=data.frame(out,stringsAsFactors = F)
  out[,3]=as.numeric(out[,3])
  out[,4]=as.numeric(out[,4])
  out[,5] <- p.adjust(out[,4], method = "BH")
  colnames(out)=c("pathway", "LF index", "AUC", "p-value", "FDR")

  return(list(Uauc=Uauc, Upval=Up, summary=out))
}

