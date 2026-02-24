######
#' Identify multicellular coordinated pathways
#'
#' This function identifies coordinated and cell-type–specific pathway activities
#' across multiple cell types using a multi-view latent factor model constrained
#' by gene–pathway priors.

#' @param data_list A list of gene expression matrices, one per cell type, with genes in rows and samples in columns.
#' @param pathMat A binary gene–pathway matrix (genes × pathways) defining gene–pathway membership.
#' @param scale Logical; whether to scale gene expression within each cell type. Default = FALSE.
#' @param allGenes Logical; whether to use all genes (TRUE) or only those appearing in `pathMat` (FALSE). Default = TRUE.
#' @param k Integer or NULL; number of latent factors. If NULL, chosen automatically.
#' @param L1 Numeric or NULL; penalty controlling the reconstruction of gene loadings from pathway priors (‖Z − C·U‖²). Default = NULL (auto-selected).
#' @param L2 Numeric or NULL; ridge penalty on the shared latent activity matrix `B` (‖B‖²). Default = NULL (auto-selected).
#' @param L3 Numeric or NULL; lasso penalty on pathway weights `U` to encourage pathway sparsity. Default = NULL (auto-selected).
#' @param L4 Numeric or NULL; ridge penalty on cell-type–specific latent factors `V`. Default = NULL (auto-selected).
#' @param L3_method Character; method to combine automatically chosen L3 values across cell types (e.g., `"mean"`, `"median"`). Default = "mean".
#' @param frac Numeric (0–1); proportion of latent factors required to have pathway support. Default = 0.7.
#' @param max.iter Integer; maximum number of iterations for model fitting. Default = 350.
#' @param maxPath Integer; maximum number of pathways per latent factor. Default = 10.
#' @param minGenes Integer; minimum number of genes per pathway. Default = 10.
#' @param doCrossval Logical; whether to perform cross-validation to tune penalties. Default = TRUE.
#' @param glm_alpha Numeric; alpha parameter for elastic-net penalty in glmnet. Default = 0.9.
#' @param trace Logical; whether to display progress information. Default = FALSE.
#' @param tol Numeric; convergence threshold. Default = 1e-6.
#' @param seed Integer; random seed for reproducibility. Default = 123456.
#'

#' @return A list with the following components:
#' \describe{
#'   \item{B}{Matrix of coordinated (shared) latent factors (LFs × samples).}
#'   \item{Z_list}{List of gene–latent factor loading matrices (\(Z_v\)) for each cell type.}
#'   \item{U_list}{List of pathway–latent factor matrices (\(U_v\)) for each cell type.}
#'   \item{V_list}{List of cell-type–specific latent factor matrices (\(V_v\)).}
#'   \item{AUC_list}{List of AUC or pathway activity quality metrics per cell type.}
#'   \item{L1}{The effective L1 penalty value used in the model.}
#'   \item{L2}{The effective L2 penalty value used in the model.}
#'   \item{L3}{The effective L3 penalty value used in the model.}
#'   \item{L4}{The effective L4 penalty value used in the model.}
#'   \item{Vratio}{Variance ratio of cell-type–specific patterns (`V`) relative to coordinated patterns (`B`).}
#'   \item{heldOutGenes}{Vector of genes held out during cross-validation (if `doCrossval = TRUE`).}
#' }
#'
#' @export



mcCoordPath <- function(data_list, pathMat,scale = FALSE,allGenes = TRUE,
                        k = NULL, L1 = NULL, L2 = NULL, L3 = NULL, L4 = NULL,L3_method=c("mean","min","max","median"),
                        frac = 0.7,  maxPath = 10,minGenes = 10,
                        doCrossval = TRUE,  glm_alpha = 0.9,
                          max.iter = 350, trace = FALSE, tol = 1e-6, seed = 123456
) {


  num_views <- length(data_list)


  data_list <- lapply(data_list, function(x) {
    t(scale(t(x), center = TRUE, scale = scale))
  })

  Y_list <- list()
  C_list <- list()
  pathMat_C_list <- list()
  heldOutGenes <- list()

  # Prepare Y and C matrices per view
  for (i in seq_len(num_views)) {
    if (!allGenes) {
      # Keep only genes common to both
      common_genes <- intersect(rownames(data_list[[i]]), rownames(pathMat))
      message("View ", i, ": Selecting common genes: ", length(common_genes))

      Y_list[[i]] <- data_list[[i]][common_genes, , drop = FALSE]
      C_list[[i]] <- pathMat[common_genes, , drop = FALSE]
      pathMat_C_list[[i]] <- C_list[[i]]

    } else {
      # Expand pathMat to include all genes from data_list[[i]]
      extra_genes <- setdiff(rownames(data_list[[i]]), rownames(pathMat))
      if (length(extra_genes) > 0) {
        eMat <- matrix(0,
                       nrow = length(extra_genes),
                       ncol = ncol(pathMat),
                       dimnames = list(extra_genes, colnames(pathMat)))
        pathMat <- rbind(pathMat, eMat)
      }

      # Remove genes in pathMat not present in data_list[[i]]
      common_genes <- intersect(rownames(data_list[[i]]), rownames(pathMat))
      Y_list[[i]] <- data_list[[i]][common_genes, , drop = FALSE]
      C_list[[i]] <- pathMat[common_genes, , drop = FALSE]
      pathMat_C_list[[i]] <- C_list[[i]]
    }
  }

  # Cross-validation mask for pathMat per view
  if (doCrossval) {
    if (!is.null(seed)) set.seed(seed)
    for (i in seq_len(num_views)) {
      C_cv <- C_list[[i]]
      for (j in seq_len(ncol(C_cv))) {
        pos <- which(C_cv[, j] > 0)
        if (length(pos) >= minGenes) {
          held_out <- sample(pos, floor(length(pos) / 5))
          C_cv[held_out, j] <- 0
          heldOutGenes[[paste0("View", i, "_", colnames(C_cv)[j])]] <- rownames(C_cv)[held_out]
        }
      }
      C_list[[i]] <- C_cv
    }
  }

  # Compute  Chat (pseudoinverse ridge or solve)

  #Chat_list <- lapply(C_list, function(C) pinv.ridge(crossprod(C), 5) %*% t(C))
  Chat_list <- lapply(C_list, function(C) solve(crossprod(C) + 5 * diag(ncol(C)), t(C)))

  # Rename rows to avoid duplication when stacking
  for (i in seq_len(num_views)) {
    rownames(Y_list[[i]]) <- paste0("view", i, "_", rownames(Y_list[[i]]))
  }

  # Combine all views vertically
  Y_all <- do.call(rbind, Y_list)
  ns <- ncol(Y_all)

  # Compute SVD

  message("Computing SVD on combined data")
  if (ns > 500) {
    set.seed(123456)
    svdres <- rsvd::rsvd(Y_all, k = min(ns, max(200, ns / 4)), q = 3)
  } else {
    svdres <- svd(Y_all)
  }


  if (is.null(k)) {
    k <- min(num.pc(svdres) * 2, floor(ns * 0.9))
    message("k is set to ", k)
  }
  if (is.null(L2)) {
    L2 <- svdres$d[k]
    message("L2 is set to ", L2)
  }
  if (is.null(L1)) {
    L1 <- svdres$d[k] / 2
    message("L1 is set to ", L1)
  }
  if (is.null(L4)) {
    L4 <- svdres$d[k]
    message("L4 is set to ", L4)
  }

  # Initialize B and Z
  B <- t(svdres$v[, 1:k] %*% diag(svdres$d[1:k]))
  Z_all <- (Y_all %*% t(B)) %*% solve(tcrossprod(B) + L1 * diag(k))

  # Avoid all-negative columns in Z and flip sign accordingly
  for (j in seq_len(ncol(Z_all))) {
    if (all(Z_all[, j] <= 0)) {
      Z_all[, j] <- -Z_all[, j]
      B[j, ] <- -B[j, ]
    }
  }
  Z_all[Z_all < 0] <- 0

  # Split Z_all back into list per view
  Z_list <- vector("list", num_views)
  for (i in seq_len(num_views)) {
    rows <- grep(paste0("^view", i, "_"), rownames(Z_all))
    Z_i <- Z_all[rows, , drop = FALSE]
    # Remove view prefix from rownames to match original gene names
    new_rownames <- sub(paste0("^view", i, "_"), "", rownames(Z_i))
    rownames(Z_i) <- new_rownames
    rownames(Y_list[[i]]) <- new_rownames
    Z_list[[i]] <- Z_i
  }

  # Initialize U and V lists
  U_list <- lapply(C_list, function(C) matrix(0, ncol = k, nrow = ncol(C),
                                              dimnames = list(colnames(C), paste0("LF", 1:k))))
  V_list <- lapply(Y_list, function(Y) matrix(0, nrow = k, ncol = ncol(Y),
                                              dimnames = list(paste0("LF", 1:k), colnames(Y))))

  L3_given <- !is.null(L3)
  iter_full_start <- 20
  iter_full <- iter_full_start
  ratio <- NA

  for (it in 1:max.iter) {

    if (it >= iter_full_start) {
      if (it == iter_full && !L3_given) {
        # Update U and estimate L3
        U_tmp <- lapply(seq_len(num_views), function(i) {
          solveU(Z_list[[i]], Chat_list[[i]], C_list[[i]],  glm_alpha, maxPath,
                 target.frac = frac)
        })
        U_list <- lapply(U_tmp, function(u) u$U)
        #L3 <- mean(sapply(U_tmp, function(u) u$L3))
        L3_values <- sapply(U_tmp, function(u) u$L3)
        L3_method <- match.arg(L3_method)
        L3 <- switch(L3_method,
                     min    = min(L3_values),
                     mean   = mean(L3_values),
                     max    = max(L3_values),
                     median = median(L3_values))



        message("New L3 estimated: ", L3)

        # Update V and L4
        Vsq_sum <- 0
        for (i in seq_len(num_views)) {
          Z <- Z_list[[i]]
          Y <- Y_list[[i]]
          lambda <- 1e-6
          denom_V <- (1 + L4) * crossprod(Z) + lambda * diag(k)
          num_V <- crossprod(Z, Y - Z %*% B)
          V_list[[i]] <- solve(denom_V, num_V)
          Vsq_sum <- Vsq_sum + sum(V_list[[i]]^2)
        }

        ratio <- sqrt(Vsq_sum) / sqrt(sum(B^2))
        if (ratio < 0.01) {
          L4 <- max(L4 / 2, 1e-6)
          message("New L4 is estimated to ", L4)
        } else if (ratio > 0.5) {
          L4 <- L4 * 2
          message("New L4 is estimated to ", L4)
        }

        iter_full <- iter_full + iter_full_start

      } else {
        # Update U with fixed L3
        U_list <- lapply(seq_len(num_views), function(i) {
          solveU(Z_list[[i]], Chat_list[[i]], C_list[[i]],  glm_alpha, maxPath,
                 L3 = L3)
        })

        # Update V with fixed L4 every iteration after warm-up
        for (i in seq_len(num_views)) {
          Z <- Z_list[[i]]
          Y <- Y_list[[i]]
          lambda <- 1e-6
          denom_V <- (1 + L4) * crossprod(Z) + lambda * diag(k)
          num_V <- crossprod(Z, Y - Z %*% B)
          V_list[[i]] <- solve(denom_V, num_V)
        }
      }
    }

    # Update Z per view
    for (i in seq_len(num_views)) {
      B_plus_V <- B + V_list[[i]]
      if (it >= iter_full_start) {
        Z_num <- Y_list[[i]] %*% t(B_plus_V) + L1 * C_list[[i]] %*% U_list[[i]]
      } else {
        Z_num <- Y_list[[i]] %*% t(B)
      }
      denom_Z <- solve(tcrossprod(B_plus_V) + L1 * diag(k))
      Z_new <- Z_num %*% denom_Z
      Z_new[Z_new < 0] <- 0
      Z_list[[i]] <- Z_new
    }

    # Update B across all views
    B_num <- Reduce(`+`, lapply(seq_len(num_views), function(i) crossprod(Z_list[[i]], Y_list[[i]])))
    B_denom <- Reduce(`+`, lapply(seq_len(num_views), function(i) crossprod(Z_list[[i]]))) + L2 * diag(k)
    B_new <- solve(B_denom, B_num)

    Bdiff <- sum((B - B_new)^2) / sum(B^2)
    if (trace) message(sprintf("iter %d: Bdiff = %.6f", it, Bdiff))

    B <- B_new

    if (it > iter_full_start && Bdiff < tol) {
      message("Converged at iteration ", it)
      break
    }
  }

  # Compute AUC for each view
  AUC_list <- vector("list", num_views)
  for (i in seq_len(num_views)) {
    if (doCrossval) {
      AUC_list[[i]] <- crossVal(list(U = U_list[[i]], Z = Z_list[[i]]), Y_list[[i]], pathMat_C_list[[i]], C_list[[i]])
    } else {
      AUC_list[[i]] <- getAUC(list(U = U_list[[i]], Z = Z_list[[i]],L1=L1,L2=L2), Y_list[[i]], pathMat_C_list[[i]])
    }
    message(sprintf("View %d: %d LFs with AUC > 0.70", i, sum(apply(AUC_list[[i]]$Uauc, 2, max) > 0.70)))
  }

  return(list(
    B = B,
    Z_list = Z_list,
    U_list = U_list,
    V_list = V_list,
    AUC_list = AUC_list,
    L1=L1,
    L2=L2,
    L3 = L3,
    L4=L4,
    Vratio = ratio,
    heldOutGenes = heldOutGenes,
    scale=scale
  ))
}



#####


#' Compute variance explained by latent factors
#'
#' This function calculates the proportion of variance in each cell type’s gene expression
#' that is explained by each coordinated (shared) and cell-type–specific latent factor
#' identified by `mcCoordPath`.
#'
#' @param data_list A list of gene expression matrices, one per cell type, with genes in rows and samples in columns.
#' @param coPathRes The result object from \code{mcCoordPath}, containing at least the elements \code{Z_list}, \code{B}, and \code{V_list}.
#' @param view_names Character vector specifying the names of cell types corresponding to the elements in \code{data_list}.
#'
#' @return A numeric matrix with latent factors in rows and cell types in columns.
#' Each entry represents the proportion of variance in that cell type’s expression data
#' explained by the corresponding latent factor.
#'
#' @export
compute_variance_explained <- function(data_list, coPathRes,  view_names=NULL){

  Z_list<-coPathRes$Z_list
  V_list<-coPathRes$V_list
  B<-coPathRes$B
  scale<-coPathRes$scale


  k <- nrow(B)  # number of LFs

  if (is.null(view_names)){
    if(!is.null(names(data_list))){
      view_names=names(data_list)}else {view_names <- paste0("View", seq_along(data_list))}
  }




  var_explained_mat <- matrix(0, nrow = k, ncol = length(data_list),
                              dimnames = list(paste0("LF", 1:k), view_names))


  data_list <- lapply(data_list, function(x) {
    t(scale(t(x), center = TRUE, scale = scale))
  })



  for (i in seq_along(data_list)) {
    Y <- data_list[[i]]
    Z <- Z_list[[i]]
    V <- V_list[[i]]

    # Align by common genes
    common_genes <- intersect(rownames(Y), rownames(Z))
    if (length(common_genes) == 0) {
      warning(paste("No overlapping genes in view", i))
      next
    }

    Y_sub <- Y[common_genes, , drop = FALSE]
    Z_sub <- Z[common_genes, , drop = FALSE]
    B_plus_V <- B + V

    total_var <- sum(Y_sub^2)

    for (l in 1:k) {
      recon_l <- Z_sub[, l, drop = FALSE] %*% B_plus_V[l, , drop = FALSE]
      #var_l <- sum(recon_l^2)
      sse_l<-sum((Y_sub-recon_l)^2)
      var_explained_mat[l, i] <- 1-sse_l / total_var
    }
  }

  return(var_explained_mat)
}



#' Infer ligand activity to explain changes in a gene list
#'
#' This function evaluates how well each ligand's predicted target gene profile
#'  can explain the observed changes in a given set of genes.
#' It computes ligand activity scores using either a Wilcoxon rank-sum test or
#' a permutation-based approach.
#'
#' @param genelist Character vector of genes of interest (e.g., differentially expressed genes in receiver cells).
#' @param ligand_target_matrix Matrix of predicted transcriptional influences
#'   of ligands on target genes, typically obtained from \pkg{nichenetr}.
#'   Rows represent ligands, columns represent target genes, and values represent
#'   regulatory potential scores.
#' @param ligands Character vector of potential ligands to test.
#' @param background Character vector of background genes to use as a reference.
#' @param method Character string specifying the method to infer ligand activity;
#'   either `"wilcox"` or `"permutation"`. Default = `"wilcox"`.
#' @param nperm Integer; number of permutations to perform if `method = "permutation"`.
#'   Default = 1000.
#'
#' @return A data frame (or matrix) containing predicted ligand activity statistics,
#'   including:
#'   \itemize{
#'     \item \code{AUC}: Area under the ROC curve for each ligand.
#'     \item \code{AUPR}: Area under the precision–recall curve.
#'     \item \code{mean_diff}: Mean difference in target gene regulatory scores
#'       between target and background genes.
#'     \item \code{p_value}: Significance of the activity score (based on Wilcoxon or permutation test).
#'   }
#'

#'
#' @export

infer_ligand_activity <- function(genelist,
                                  ligand_target_matrix,
                                  ligands,
                                  background,
                                  method = c("wilcox", "permutation"),
                                  n_perm = 1000) {
  method <- match.arg(method)

  # Filter valid targets
  valid_targets <- intersect(rownames(ligand_target_matrix), background)
  geneset <- intersect(genelist, valid_targets)

  if(length(geneset) == 0) stop("No genes in genelist match background/ligand_target_matrix.")

  lig_results <- data.frame(ligand = ligands,
                            auc = NA,
                            aupr = NA,
                            mean_diff = NA,
                            pval = NA,
                            stringsAsFactors = FALSE)

  for (lig in ligands) {
    scores <- ligand_target_matrix[valid_targets, lig]  # ligand in columns
    gene_indicator <- names(scores) %in% geneset

    # --- AUC ---
    auc <- tryCatch({
      pROC::auc(response = gene_indicator,
                predictor = scores,
                quiet = TRUE)
    }, error = function(e) NA)

    # --- auPR ---
    aupr <- tryCatch({
      PRROC::pr.curve(scores.class0 = scores[gene_indicator],
                      scores.class1 = scores[!gene_indicator],
                      curve = FALSE)$auc.integral
    }, error = function(e) NA)

    # --- Mean difference ---
    mean_diff <- mean(scores[gene_indicator], na.rm = TRUE) - mean(scores[!gene_indicator], na.rm = TRUE)

    # --- p-value ---
    if (method == "wilcox") {
      pval <- tryCatch({
        wilcox.test(scores[gene_indicator],
                    scores[!gene_indicator],alternative="greater")$p.value
      }, error = function(e) NA)
    } else if (method == "permutation") {

      set_size <- length(geneset)
      perm_means <- replicate(n_perm, {
        perm_genes <- sample(valid_targets, set_size)
        mean(scores[perm_genes]) - mean(scores[setdiff(valid_targets,perm_genes)])
      })
      #pval <- mean(perm_means >= mean_diff, na.rm = TRUE)
      perm_mean <- mean(perm_means, na.rm = TRUE)
      perm_sd   <- sd(perm_means, na.rm = TRUE)

      if (is.na(perm_sd) || perm_sd == 0) {
        pval <- NA
      } else {
        z <- (mean_diff - perm_mean) / perm_sd
        pval <- 1 - pnorm(z)   # P(mean_diff > perm distribution)
      }


    }

    # Store results
    lig_results[lig_results$ligand == lig, "auc"] <- as.numeric(auc)
    lig_results[lig_results$ligand == lig, "aupr"] <- as.numeric(aupr)
    lig_results[lig_results$ligand == lig, "mean_diff"] <- mean_diff
    lig_results[lig_results$ligand == lig, "pval"] <- pval
  }

  return(lig_results)
}





#####
#' Infer ligand–receptor interactions for a given latent factor
#'
#' This function identifies potential ligand–receptor interactions between cell types
#' that are associated with a specific latent factor from a coordinated multi-view model.
#' Optionally, it can evaluate ligand activity using NicheNet's ligand–target matrix.
#'
#' @param coPathRes A result object from \code{mcCoordPath}
#' @param lf_index Integer specifying the index of the latent factor to analyze.
#' @param topn Integer; number of top genes contributing to the selected latent factor
#'   to be considered for ligand–receptor interaction analysis.
#' @param lr_df Data frame of ligand–receptor pairs, with ligands in the first column
#'   and receptors in the second column.
#' @param ligand_target_matrix Optional matrix of predicted transcriptional influences
#'   of ligands on target genes, typically obtained from \pkg{NicheNetR}. Rows represent
#'   ligands, columns represent target genes, and values represent regulatory potential scores.
#'   If \code{NULL}, ligand activity is not predicted.
#' @param method Character string specifying the method to infer ligand activity when
#'   \code{ligand_target_matrix} is provided. Either \code{"wilcox"} or \code{"permutation"}.
#'   Default = \code{"wilcox"}.
#' @param view_names Character vector specifying the names of cell types corresponding
#'   to the elements in \code{data_list}.
#'
#' @return A data frame (or matrix) containing inferred ligand–receptor interactions, including:
#' \itemize{
#'   \item \code{ligand}: Ligand gene name.
#'   \item \code{receptor}: Receptor gene name.
#'   \item \code{sender}: Cell type expressing the ligand.
#'   \item \code{receiver}: Cell type expressing the receptor.
#'   \item \code{AUC}: Area under the ROC curve for each ligand.
#'   \item \code{AUPR}: Area under the precision–recall curve.
#'   \item \code{mean_diff}: Mean difference in regulatory potential between target and background genes.
#'   \item \code{p_value}: Significance of ligand activity (based on Wilcoxon or permutation test).
#' }
#'

#' @export

infer_lr_interactions <- function(
    coPathRes,
    lf_index = 1,
    topn = 200,
    lr_df = LRpair,
    ligand_target_matrix = NULL,
    method = c("wilcox", "permutation"),
    view_names = NULL
) {
  method <- match.arg(method)

  lr_pairs <- data.frame(
    ligand = as.character(lr_df[[1]]),
    receptor = as.character(lr_df[[2]]),
    stringsAsFactors = FALSE
  )

  Z_list <- coPathRes$Z_list
  gene_sets <- lapply(Z_list, function(Z) {
    Z_lf <- Z[, lf_index]
    names(sort(Z_lf, decreasing = TRUE))[1:min(topn, length(Z_lf))]
  })

  if (!is.null(view_names) && length(Z_list) != length(view_names)) {
    stop("view_names do not match the number of elements in Z_list.")
  }

  if (is.null(view_names)) {
    names(gene_sets) <- paste0("View", seq_along(Z_list))
  } else {
    names(gene_sets) <- view_names
  }

  interaction_results <- list()

  for (sender in names(gene_sets)) {
    for (receiver in names(gene_sets)) {
      ligands <- gene_sets[[sender]]
      receptors <- gene_sets[[receiver]]

      lr_matches <- lr_pairs[lr_pairs$ligand %in% ligands & lr_pairs$receptor %in% receptors, ]

      # If ligand_target_matrix is provided, filter and infer activity
      if (!is.null(ligand_target_matrix)) {
        lr_matches <- lr_matches[lr_matches$ligand %in% colnames(ligand_target_matrix), , drop = FALSE]
      }

      if (nrow(lr_matches) > 0) {
        Z_sender <- Z_list[[which(names(gene_sets) == sender)]]
        Z_receiver <- Z_list[[which(names(gene_sets) == receiver)]]

        lr_matches$ligand_score <- Z_sender[lr_matches$ligand, lf_index]
        lr_matches$receptor_score <- Z_receiver[lr_matches$receptor, lf_index]
        lr_matches$sender <- sender
        lr_matches$receiver <- receiver

        #  Only call infer_ligand_activity if matrix is provided
        if (!is.null(ligand_target_matrix)) {
          lig_act <- infer_ligand_activity(
            receptors,
            ligand_target_matrix,
            unique(lr_matches$ligand),
            rownames(Z_receiver),
            method = method
          )
          lr_matches <- merge(lr_matches, lig_act, by = "ligand", all.x = TRUE)
        }

        interaction_results[[paste(sender, receiver, sep = "_")]] <- lr_matches
      }
    }
  }

  final_lr_result <- do.call(rbind, interaction_results)
  rownames(final_lr_result) <- NULL

  # Only compute adjusted p-values and ranks if pval exists
  if ("pval" %in% colnames(final_lr_result)) {
    final_lr_result$pval_adj <- p.adjust(final_lr_result$pval, method = "BH")
    final_lr_result$logp <- -log10(final_lr_result$pval_adj + 1e-300)
    md_z <- scale(final_lr_result$mean_diff)
    logp_z <- scale(final_lr_result$logp)
    final_lr_result$rank_score <- md_z + logp_z
  }

  return(final_lr_result)
}






####
#' Extend pathway–latent factor associations to new pathways
#'
#' This function updates or extends the pathway–latent factor mappings (\code{U_list})
#' from a previous \code{mcCoordPath} result by incorporating a new pathway–gene
#' annotation matrix. It can be used to test whether additional pathways are associated
#' with existing coordinated latent factors.
#'
#' @param data_list A list of gene expression matrices, one per cell type, with
#'   genes in rows and samples in columns.
#' @param coPathRes The result object from \code{mcCoordPath}
#' @param newpathMat A binary gene–pathway matrix (genes in rows, pathways in columns)
#'   representing the new pathways to be evaluated.
#' @param target.frac Numeric value between 0 and 1; the fraction of latent factors
#'   expected to have pathway associations. Default = \code{0.7}.
#' @param maxPath Integer; maximum number of pathways to retain for each latent
#'   factor after extension. Default = \code{10}.
#' @param minGenes Integer; minimum number of genes required in each pathway
#'   to be considered. Default = \code{10}.
#' @param useSparse Logical; whether to use sparse matrix operations for efficiency.
#'   Default = \code{TRUE}.
#'
#' @return A list with the same structure as the original \code{mcCoordPath} result,
#'   but with updated \code{U_list} matrices reflecting the extended pathway annotations.
#'

#' @export


extendU <- function(data_list, coPathRes, newpathMat, target.frac = 0.7, maxPath = 10, minGenes = 10, useSparse = TRUE) {

  Z_list <- coPathRes$Z_list
  C_list <- list()
  Chat_list <- list()
  U_list <- list()
  AUC_list <- list()
  scale=coPathRes$scale
  data_list <- lapply(data_list, function(x) {
    t(scale(t(x), center = TRUE, scale = scale))
  })


  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = length(Z_list), style = 3)

  for (i in seq_along(Z_list)) {

    # Match gene order
    data_list[[i]] <- data_list[[i]][rownames(Z_list[[i]]), , drop = FALSE]

    # Keep only common genes
    common_genes <- intersect(rownames(Z_list[[i]]), rownames(newpathMat))
    Z_i <- Z_list[[i]][common_genes, , drop = FALSE]
    data_i <- data_list[[i]][common_genes, , drop = FALSE]

    # Subset newpathMat and remove pathways with too few genes
    C_i <- newpathMat[common_genes, , drop = FALSE]
    gene_counts <- colSums(C_i != 0)
    keep_paths <- which(gene_counts >= minGenes)
    C_i <- C_i[, keep_paths, drop = FALSE]

    if (ncol(C_i) == 0) {
      warning(sprintf("No pathways left in view %d after minGenes filter.", i))
      U_list[[i]] <- NULL
      AUC_list[[i]] <- NULL
      setTxtProgressBar(pb, i)
      next
    }

    # Convert to sparse matrix if requested
    if (useSparse) {C_i <- as(C_i, "dgCMatrix")}

    # Compute ridge projection
    CtC <- Matrix::crossprod(C_i)
    Chat_i <- Matrix::solve(CtC + 5 * Matrix::Diagonal(ncol(C_i)), Matrix::t(C_i))

    # Compute U
    U_i <- solveU(Z_i, Chat_i, C_i, target.frac = target.frac, maxPath = maxPath)$U

    # Compute AUC
    AUC_i <- getAUC(list(U = U_i, Z = Z_i, L1 = coPathRes$L1, L2 = coPathRes$L2), data_i, C_i)

    # Save results
    Z_list[[i]] <- Z_i
    C_list[[i]] <- C_i
    Chat_list[[i]] <- Chat_i
    U_list[[i]] <- U_i
    AUC_list[[i]] <- AUC_i

    # Update progress bar
    setTxtProgressBar(pb, i)
  }

  # Close progress bar
  close(pb)

  return(list(
    B = coPathRes$B,
    Z_list = Z_list,
    U_list = U_list,
    AUC_list = AUC_list
  ))
}







