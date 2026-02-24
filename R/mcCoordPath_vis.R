#' Plot multicellular coordinated pathways
#'
#' Generate a heatmap showing pathway scores across multiple cell types/views
#' from a mcCoordPath result.
#' This function visualizes the relationship between latent factors and pathway
#' enrichment, highlighting significant associations.
#'
#' @param coPathRes A result object from \code{mcCoordPath}, which must include
#'   a list \code{AUC_list} with summary tables for each cell type/view.
#' @param auc_cutoff Numeric. Minimum AUC threshold to retain significant
#'   pathway–latent factor associations. Default is \code{0.6}.
#' @param fdr_cutoff Numeric. Maximum FDR threshold to retain significant
#'   pathway–latent factor associations. Default is \code{0.05}.
#' @param view_names Character vector (optional). Names of the cell types to
#'   display in the heatmap. If \code{NULL}, names will default to "View1", "View2", etc.
#' @param LF_v Numeric vector (optional). Specific latent factors (LFs) to
#'   highlight or subset in the heatmap. If \code{NULL}, all latent factors are shown.
#'
#' @return A \code{ggplot} or \code{pheatmap} object showing a heatmap of pathway
#'   enrichment scores. Views are displayed as rows (grouped by latent factors),
#'   and pathways as columns. Significant associations are annotated with stars
#'   based on FDR thresholds:
#'   \itemize{
#'     \item "***" : FDR < 1e-4
#'     \item "**"  : 1e-4 <= FDR < 1e-3
#'     \item "*"   : 1e-3 <= FDR < 1e-2
#'   }
#'   The function returns invisibly a list containing:
#'   \itemize{
#'     \item \code{mat} – the ordered AUC matrix used for plotting
#'     \item \code{fdr} – the corresponding FDR matrix
#'     \item \code{sig_labels} – matrix of significance stars
#'   }
#'
#' @details
#' The function combines results across multiple cell types (views) and filters
#' pathway–latent factor associations based on the AUC and FDR cutoffs. Rows
#' are ordered by latent factor index, and gaps are added between different LFs
#' in the heatmap for clarity. Optionally, specific LFs can be selected via
#' \code{LF_v}.
#'
#' @examples
#' \dontrun{
#' plot_multiCellular_pathways(coPathRes = result)
#'
#' }
#'
#' @export

plot_multicellular_pathways<- function(coPathRes,
                                           auc_cutoff = 0.6,
                                           fdr_cutoff = 0.05,
                                           view_names = NULL,
                                           LF_v = NULL) {
  AUC_list <- coPathRes$AUC_list

  # 1. Prepare view names
  if (is.null(names(AUC_list)) & is.null(view_names)) {
    names(AUC_list) <- paste0("View", seq_along(AUC_list))
  } else if (!is.null(view_names) & length(AUC_list) == length(view_names)) {
    names(AUC_list) <- view_names
  } else {
    warning("The length of view names does not match the number of views in `coPathRes`.")
  }

  # 2. Combine summaries
  AUC_summary_df <- do.call(rbind, lapply(seq_along(AUC_list), function(i) {
    view_name <- names(AUC_list)[i]
    auc_summary <- AUC_list[[i]]$summary
    if (nrow(auc_summary) == 0) return(NULL)
    auc_summary$View <- view_name
    auc_summary
  }))

  # 3. Filter by cutoff
  auc_filt <- AUC_summary_df %>%
    filter(AUC > auc_cutoff, FDR < fdr_cutoff)

  if (nrow(auc_filt) == 0) {
    message("No pathways passed the cutoff (AUC > ", auc_cutoff, ", FDR < ", fdr_cutoff, "). Returning NULL.")
    return(NULL)
  }

  # 4. Optional: filter by LF selection
  if (!is.null(LF_v)) {
    lf_idx <- as.numeric(str_extract(LF_v, "\\d+"))
    auc_filt <- auc_filt %>% filter(`LF index` %in% lf_idx)
    if (nrow(auc_filt) == 0) {
      message("No LFs in ", paste(LF_v, collapse = ", "), " passed the criteria.")
      return(NULL)
    }
  }

  # 5. Add row labels
  auc_filt <- auc_filt %>%
    mutate(RowID = paste0(View, "_LF", `LF index`))

  # 6. Pivot to AUC matrix
  mat <- auc_filt %>%
    select(RowID, pathway, AUC) %>%
    pivot_wider(names_from = pathway, values_from = AUC, values_fill = NA) %>%
    column_to_rownames("RowID") %>%
    as.data.frame()
  mat[is.na(mat)] <- 0

  # 7. FDR matrix
  fdr_mat_df <- auc_filt %>%
    select(RowID, pathway, FDR) %>%
    pivot_wider(names_from = pathway, values_from = FDR, values_fill = NA) %>%
    column_to_rownames("RowID") %>%
    as.data.frame()
  fdr_mat_df <- fdr_mat_df[, colnames(mat), drop = FALSE]

  # 8. Convert FDR to stars
  fdr_to_stars <- function(fdr) {
    stars <- rep("", length(fdr))
    stars[fdr < 1e-4] <- "***"
    stars[fdr >= 1e-4 & fdr < 1e-3] <- "**"
    stars[fdr >= 1e-3 & fdr < 1e-2] <- "*"
    stars[is.na(fdr)] <- ""
    stars
  }

  sig_labels <- matrix(fdr_to_stars(as.vector(as.matrix(fdr_mat_df))),
                       nrow = nrow(fdr_mat_df),
                       ncol = ncol(fdr_mat_df),
                       dimnames = dimnames(fdr_mat_df))

  # 9. Order rows
  lf_order <- as.numeric(str_extract(rownames(mat), "(?<=LF)\\d+"))
  row_order <- order(lf_order)
  mat_ordered <- mat[row_order, ]
  sig_labels_ordered <- sig_labels[row_order, , drop = FALSE]
  lf_order_sorted <- lf_order[row_order]
  gaps <- which(diff(lf_order_sorted) != 0)

  # 10. Plot heatmap
  pheatmap(
    mat_ordered,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    gaps_row = gaps,
    main = paste0("Combined Pathway AUCs (AUC>", auc_cutoff, ", FDR<", fdr_cutoff, ")"),
    color = viridis(100),
    display_numbers = sig_labels_ordered,
    number_color = "black",
    fontsize_number = 10,
    fontsize_row = 8,
    fontsize_col = 6,
    na_col = "grey90"
  )

  # Return matrices invisibly for later use
  invisible(list(
    mat = mat_ordered,
    fdr = fdr_mat_df[row_order, , drop = FALSE],
    sig_labels = sig_labels_ordered
  ))
}







#' Plot all coordinated latent factor scores across samples
#'
#' Visualize the cooridnated latent factor (LF) scores stored in the B matrix
#' from a mcCoordPath result.
#' Each LF can be plotted either as a line plot
#' or boxplot grouped by sample group.
#'
#' @param coPathRes A result object from \code{mcCoordPath}
#' @param sample_group Optional vector indicating sample groups for boxplot
#'   comparison. Must be the same length as the number of samples (columns of B).
#' @param plot_type Character, one of \code{"auto"}, \code{"line"}, or
#'   \code{"boxplot"} (default \code{"auto"}). If \code{"auto"}, line plot is
#'   used when \code{sample_group} is NULL, otherwise boxplot.
#' @param ncol Number of columns when arranging multiple LF plots. Default is 1.
#' @param show_xaxis_text Logical. Whether to show x-axis text (sample names). Default is TRUE.
#' @param p_cutoff Numeric. Significance threshold for displaying pairwise p-values
#'   in boxplots. Default is 0.05.
#' @param add_violin Logical. Whether to add violin plot behind the boxplot. Default is FALSE
#' @return A \code{patchwork} object arranging plots of all latent factors.
#'   Each plot shows the LF scores across samples. Boxplots may display
#'   significance stars for pairwise comparisons.
#'
#' @details
#' This function loops over all latent factors (rows of B) and generates individual
#' plots. If \code{plot_type = "line"}, LF scores are connected across samples.
#' If \code{plot_type = "boxplot"}, samples are grouped by \code{sample_group} and
#' pairwise t-tests are performed to highlight significant differences.
#'
#' @examples
#' \dontrun{
#' plot_all_LFs(coPathRes = result)
#' }
#'
#' @export

plot_all_LFs <- function(coPathRes,
                         sample_group = NULL,
                         plot_type = c("auto", "line", "boxplot"),
                         ncol = 1,
                         show_xaxis_text = TRUE,
                         p_cutoff = 0.05,
                         add_violin = FALSE) {

  plot_type <- match.arg(plot_type)
  B <- coPathRes$B

  if (is.null(colnames(B))) {
    colnames(B) <- paste0("s", 1:ncol(B))
  }

  samples <- colnames(B)

  # Determine plot type
  if (plot_type == "auto") {
    plot_type <- if (is.null(sample_group)) "line" else "boxplot"
  }

  # Build data frame for all LFs
  df_list <- lapply(rownames(B), function(lf) {
    data.frame(
      LF = lf,
      Sample = samples,
      Score = as.numeric(B[lf, ]),
      Group = if (!is.null(sample_group)) sample_group else samples,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, df_list)
  df$Group <- factor(df$Group)

  # Base theme
  base_theme <- theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title.x = element_blank(),
      plot.margin = margin(5, 5, 15, 5)
    )

  if (!show_xaxis_text) {
    base_theme <- base_theme +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }

  # Generate plots per LF
  plots <- lapply(unique(df$LF), function(lf) {
    subdf <- df[df$LF == lf, ]
    subdf$Group <- factor(subdf$Group)
    ptitle <- lf

    ## ---------------- LINE PLOT ---------------- ##
    if (plot_type == "line") {
      return(
        ggplot(subdf, aes(x = Sample, y = Score, group = 1)) +
          geom_line(color = "steelblue") +
          geom_point(color = "steelblue", size = 1.5) +
          base_theme +
          labs(title = ptitle, y = "Score")
      )
    }

    ## ---------------- BOX + (optional) VIOLIN ---------------- ##
    groups <- levels(subdf$Group)
    comparisons <- combn(groups, 2, simplify = FALSE)

    # Compute p-values
    sig_comparisons <- lapply(comparisons, function(x) {
      vals1 <- subdf$Score[subdf$Group == x[1]]
      vals2 <- subdf$Score[subdf$Group == x[2]]
      pval <- t.test(vals1, vals2)$p.value
      if (pval < p_cutoff) return(x) else return(NULL)
    })
    sig_comparisons <- Filter(Negate(is.null), sig_comparisons)

    y_max <- max(subdf$Score, na.rm = TRUE)
    y_range <- diff(range(subdf$Score, na.rm = TRUE))

    # Start plot
    base_plot <- ggplot(subdf, aes(x = Group, y = Score, fill = Group))

    # Add violin if requested
    if (add_violin) {
      base_plot <- base_plot +
        geom_violin(trim = FALSE, alpha = 0.3, width = 0.9, color = NA)
    }

    # Add boxplot
    base_plot <- base_plot +
      geom_boxplot(
        width = if (add_violin) 0.25 else 0.6,
        outlier.size = 0.5
      ) +
      base_theme +
      theme(legend.position = "none") +
      labs(title = ptitle, y = "Score")

    # Add p-values if significant pairs exist
    if (length(sig_comparisons) > 0) {
      base_plot <- base_plot +
        stat_compare_means(
          method = "t.test",
          comparisons = sig_comparisons,
          label = "p.signif",
          step.increase = 0.1,
          size = 4
        ) +
        expand_limits(y = y_max + 0.2 * y_range)
    }

    base_plot
  })

  wrap_plots(plots, ncol = ncol)
}


#######

#' Plot variance explained by latent factors across cell types (views)
#'
#' Generate a heatmap showing the proportion of variance explained by latent
#' factors in each view or cell type.
#'
#' @param data_list A list of data matrices for each cell type (view).
#' @param coPathRes A result object from \code{mcCoordPath} containing the latent
#'   factor matrices.
#' @param view_names Optional character vector specifying the names of the cell types.
#' @param display_numbers Logical. Whether to display numerical values in the heatmap cells. Default is TRUE.
#' @param cols Color palette for the heatmap. Default is a blue gradient from
#'   "#deebf7" to "#3182bd".
#'
#' @return Invisibly returns a list with components:
#'   \itemize{
#'     \item \code{heatmap}: the \code{pheatmap} object
#'     \item \code{values}: matrix of variance explained values
#'   }
#'
#' @details
#' Variance explained is computed for each latent factor across each view. The
#' resulting matrix is visualized as a heatmap, optionally showing numeric
#' values for clarity.
#'
#' @examples
#' \dontrun{
#' plot_variance_explained(data_list = my_data, coPathRes = result)
#' }
#'
#' @export

plot_variance_explained<-function(data_list, coPathRes, view_names=NULL, display_numbers=T,cols=colorRampPalette(c("#deebf7", "#3182bd"))(10)){



  var_explained<-compute_variance_explained(data_list,coPathRes = coPathRes, view_names = view_names)


  var_fig <- pheatmap(var_explained,
                      cluster_rows = FALSE,
                      cluster_cols = FALSE,
                      color = cols,
                      display_numbers = display_numbers)

  invisible(list(
    heatmap = var_fig,
    values = var_explained
  ))
}


####
#' Plot correlations of latent factors across cell types (views)
#'
#' Generate a series of heatmaps showing the correlation of latent factors (LFs)
#' between multiple cell types. Each heatmap corresponds to one latent factor,
#' displaying the pairwise correlations of that factor across all cell types.
#'
#' @param coPathRes A result object from \code{mcCoordPath}
#' @param view_names Optional character vector specifying names of the cell types.
#'   If \code{NULL}, defaults to "View1", "View2", etc.
#' @param ncol Integer. Number of columns when arranging the heatmaps in a grid.
#'   Default is 2.
#'
#' @return A combined \code{grob} object arranging all latent factor correlation
#'   heatmaps.


#' @export
plot_Cor_Celltypes<-function(coPathRes,view_names=NULL,ncol=2){

  nview=length(coPathRes$Z_list)
  nLF=ncol(coPathRes$Z_list[[1]])
  LF_names<-paste0("LF", 1:nLF)
  if (is.null(names(coPathRes$Z_list)) & is.null(view_names)) {
    view_names <- paste0("View", 1:nview)
  }else if(!is.null(view_names) & nview==length(view_names)){
    view_names<-view_names}else {warining("the length of view names does not match the number of views in the result")}


  corresult<-array(1,dim=c(nview,nview,nLF),dimnames=list(view_names,view_names,LF_names))
  for (i in 1:(nview-1)){
    for (j in (i+1):nview){
      comm<-intersect(rownames(coPathRes$Z_list[[i]]),rownames(coPathRes$Z_list[[j]]))
      x1<-coPathRes$Z_list[[i]][comm,]
      y1<-coPathRes$Z_list[[j]][comm,]
      cortmp<-diag(cor(x1,y1))
      corresult[i,j,]<-corresult[j,i,]<-cortmp

    }

  }
  p_list <- lapply(1:nLF, function(m) pheatmap(corresult[,,m],main=paste0("LF",m),treeheight_col = 0,silent = TRUE))
  g_list <- lapply(p_list, function(p) p$gtable)

  # Return combined grob (not plot directly)
  combined <- do.call(grid.arrange, c(g_list, ncol = ncol))
  return(combined)


}






#' Plot top genes in enriched pathways for a latent factor
#'
#' Generate a heatmap of the top genes from pathways enriched in a specific latent
#' factor (LF) across multiple cell types or views. Only pathways passing AUC and
#' FDR thresholds are considered, and the top genes per pathway are selected
#' based on their latent factor loadings.
#'
#' @param coPathRes A result object from \code{mcCoordPath}
#' @param pathMat A binary matrix indicating gene membership in pathways
#' @param lf_index Integer specifying which latent factor to examine. Default is 1.
#' @param topn Integer. Maximum number of top genes per pathway to include in the heatmap. Default is 10.
#' @param auc_threshold Numeric. Minimum AUC threshold to consider a pathway enriched. Default is 0.7.
#' @param fdr_threshold Numeric. Maximum FDR threshold to consider a pathway enriched. Default is 0.05.
#' @param scale_rows Logical. If TRUE, rows of the heatmap are scaled (z-score). Default is FALSE.
#' @param view_names Optional character vector specifying names of the views/cell types.
#'   If NULL, defaults to "View1", "View2", etc.
#'
#' @return The function plots a heatmap using \code{pheatmap} showing the top genes
#'   in enriched pathways for the selected latent factor. Genes are annotated by
#'   the pathways they belong to across all cell types.

#' @export

plot_top_genes_LF_in_pathways<- function(coPathRes, pathMat, lf_index = 1, topn = 10, auc_threshold = 0.7, fdr_threshold = 0.05, scale_rows = FALSE, view_names=NULL) {


  n_views <- length(coPathRes$Z_list)

  AUC_list <- coPathRes$AUC_list


  if (is.null(names(AUC_list)) & is.null(view_names)) {
    names(AUC_list) <- paste0("View", seq_along(AUC_list))
  }else if(!is.null(view_names) & length(AUC_list)==length(view_names)){
    names(AUC_list)<-view_names}else {warining("the length of view names does not match the number of views in the res")}



  all_top_genes <- c()

  # For each view separately
  for (v in seq_len(n_views)) {
    summary_df <- coPathRes$AUC_list[[v]]$summary

    # Select pathways passing thresholds for this view & lf_index
    selected_pathways <- subset(summary_df, `LF index` == lf_index & AUC > auc_threshold & FDR < fdr_threshold)$pathway
    selected_pathways <- intersect(selected_pathways, colnames(pathMat))
    if (length(selected_pathways) == 0) next

    Z <- coPathRes$Z_list[[v]]

    # Collect top genes from each pathway for this view
    top_genes_view <- c()
    for (pw in selected_pathways) {
      genes_in_path <- rownames(pathMat)[pathMat[, pw] == 1]
      genes_in_path <- intersect(genes_in_path, rownames(Z))
      if (length(genes_in_path) == 0) next

      loadings <- Z[genes_in_path, lf_index]
      loadings <- loadings[!is.na(loadings)]
      n_g <- min(topn, length(loadings))
      top_genes_pathway <- names(sort(loadings, decreasing = TRUE))[1:n_g]
      top_genes_view <- unique(c(top_genes_view, top_genes_pathway))
    }

    all_top_genes <- unique(c(all_top_genes, top_genes_view))
  }

  if (length(all_top_genes) == 0) stop("No top genes found in any view's enriched pathways.")

  # Build heatmap matrix across all views and all top genes
  heatmap_mat <- do.call(cbind, lapply(coPathRes$Z_list, function(Zmat) {
    vals <- rep(0, length(all_top_genes))
    names(vals) <- all_top_genes
    common_genes <- intersect(rownames(Zmat), all_top_genes)
    vals[common_genes] <- Zmat[common_genes, lf_index]
    vals
  }))
  rownames(heatmap_mat) <- all_top_genes
  colnames(heatmap_mat) <- names(AUC_list)

  if (scale_rows) heatmap_mat <- t(scale(t(heatmap_mat)))

  # Annotate genes with pathways (across all views)
  # Create a combined list of which gene belongs to which pathways (from any view)
  gene_pathway_map <- list()
  for (v in seq_len(n_views)) {
    summary_df <- coPathRes$AUC_list[[v]]$summary
    selected_pathways <- subset(summary_df, `LF index` == lf_index & AUC > auc_threshold & FDR < fdr_threshold)$pathway
    selected_pathways <- intersect(selected_pathways, colnames(pathMat))
    if (length(selected_pathways) == 0) next

    for (pw in selected_pathways) {
      genes_in_path <- rownames(pathMat)[pathMat[, pw] == 1]
      genes_in_path <- intersect(genes_in_path, all_top_genes)
      for (gene in genes_in_path) {
        gene_pathway_map[[gene]] <- c(gene_pathway_map[[gene]], pw)
      }
    }
  }

  # Format annotation: collapse multiple pathways per gene with "; "
  gene_pathway <- sapply(all_top_genes, function(g) {
    if (!is.null(gene_pathway_map[[g]])) {
      paste(unique(gene_pathway_map[[g]]), collapse = "; ")
    } else {
      "None"
    }
  })
  gene_pathway <- factor(gene_pathway)

  row_annotation <- data.frame(
    pathway = gene_pathway
  )
  rownames(row_annotation) <- all_top_genes

  # Colors
  pw_levels <- levels(gene_pathway)
  n_pw <- length(pw_levels)
  palette <- grDevices::rainbow(n_pw)
  names(palette) <- pw_levels

  ann_colors <- list(
    pathway = palette
  )

  # Plot
  pheatmap(heatmap_mat,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           fontsize_row = 6,
           main = paste0("Top genes in enriched pathways per view (LF", lf_index, ")"),
           annotation_row = row_annotation,
           annotation_colors = ann_colors)
}



#######
#' Plot top genes in a specific pathway for a latent factor
#'
#' Generate a heatmap of the top genes from a specific pathway that is enriched
#' in a given latent factor (LF) across multiple cell types or views. The function filters pathways
#' by AUC and FDR thresholds and selects the top genes per view based on latent factor loadings.
#'
#' @param coPathRes A result object from \code{mcCoordPath}
#' @param pathMat A binary matrix indicating gene membership in pathways
#' @param pathName Character string specifying the pathway to plot. Partial matching
#'   is allowed.
#' @param lf_index Integer specifying which latent factor to examine. Default is 1.
#' @param topn Integer. Maximum number of top genes per view to include. Default is 20.
#' @param auc_threshold Numeric. Minimum AUC threshold to consider a pathway enriched. Default is 0.7.
#' @param fdr_threshold Numeric. Maximum FDR threshold to consider a pathway enriched. Default is 0.05.
#' @param scale_rows Logical. If TRUE, rows of the heatmap are scaled (z-score). Default is FALSE.
#' @param view_names Optional character vector specifying names of the views/cell types.
#'   If NULL, defaults to "View1", "View2", etc.
#' @param ignore_case Logical. If TRUE, pathway name matching ignores case. Default is TRUE.
#'
#' @return The function plots a heatmap of the top genes in the selected pathway
#'   across all views. It also invisibly returns a list containing:
#'   \item{heatmap}{The \code{pheatmap} object.}
#'   \item{values}{Matrix of latent factor loadings used in the heatmap.}
#'   \item{top_genes}{Character vector of selected top genes.}
#'   \item{matched_path}{The exact pathway name matched.}
#'   \item{enriched_views}{Views in which the pathway passed thresholds.}
#'

#' @export


plot_top_genes_LF_in_a_pathway <- function(coPathRes, pathMat,
                                           pathName,
                                           lf_index = 1,
                                           topn = 20,
                                           auc_threshold = 0.7,
                                           fdr_threshold = 0.05,
                                           scale_rows = FALSE,
                                           view_names = NULL,
                                           ignore_case = TRUE) {

  n_views <- length(coPathRes$Z_list)
  AUC_list <- coPathRes$AUC_list

  if (is.null(names(AUC_list)) & is.null(view_names)) {
    names(AUC_list) <- view_names<-paste0("View", seq_along(AUC_list))
  }else if(!is.null(view_names) & length(AUC_list)==length(view_names)){
    names(AUC_list)<-view_names}else {warining("the length of view names does not match the number of views in the res")}



  # --- partial matching for pathway name ---
  matched_paths <- grep(pathName, colnames(pathMat), value = TRUE, ignore.case = ignore_case)
  if (length(matched_paths) == 0) {
    stop(paste0("No pathway matches found for '", pathway_name, "'."))
  } else if (length(matched_paths) > 1) {
    message("Multiple matches found: ", paste(matched_paths, collapse = ", "),
            "\nUsing the first match: ", matched_paths[1])
    matched_path <- matched_paths[1]
  } else {
    matched_path <- matched_paths
  }

  # --- collect top genes from enriched views ---
  top_genes_all <- c()
  enriched_views <- c()

  for (v in seq_len(n_views)) {
    summary_df <- AUC_list[[v]]$summary
    sel <- summary_df$pathway == matched_path &
      summary_df$`LF index` == lf_index &
      summary_df$AUC > auc_threshold &
      summary_df$FDR < fdr_threshold

    if (!any(sel)) next

    enriched_views <- c(enriched_views, view_names[v])

    genes_in_path <- rownames(pathMat)[pathMat[, matched_path] == 1]
    Z <- coPathRes$Z_list[[v]]
    genes_in_path <- intersect(genes_in_path, rownames(Z))
    if (length(genes_in_path) == 0) next

    loadings <- Z[genes_in_path, lf_index]
    loadings <- loadings[!is.na(loadings)]
    n_g <- min(topn, length(loadings))

    top_genes <- names(sort(loadings, decreasing = TRUE))[1:n_g]
    top_genes_all <- unique(c(top_genes_all, top_genes))
  }

  if (length(top_genes_all) == 0) {
    stop(paste0("Pathway '", matched_path, "' does not pass thresholds in any view."))
  }

  # --- build heatmap matrix for top genes across ALL views ---
  heatmap_mat <- do.call(cbind, lapply(coPathRes$Z_list, function(Z) {
    vals <- rep(0, length(top_genes_all))
    names(vals) <- top_genes_all
    common <- intersect(top_genes_all, rownames(Z))
    vals[common] <- Z[common, lf_index]
    vals
  }))

  rownames(heatmap_mat) <- top_genes_all
  # IMPORTANT: set column names to the view labels we computed
  colnames(heatmap_mat) <- view_names

  if (scale_rows) heatmap_mat <- t(scale(t(heatmap_mat)))

  # --- avoid breaks error if identical ---
  range_vals <- range(heatmap_mat, na.rm = TRUE)
  if (diff(range_vals) == 0) {
    message("All heatmap values identical; adjusting range slightly to avoid 'breaks' error.")
    range_vals <- range_vals + c(-1e-6, 1e-6)
  }

  title_text <- paste0("Top genes in ", matched_path,
                       " (LF", lf_index, ") — enriched in: ",
                       if (length(enriched_views) > 0) paste(enriched_views, collapse = ", ") else "None")

  # --- plot: show column names, rotate if long, adjust fontsize if many views ---
  col_fontsize <- ifelse(nchar(paste(view_names, collapse = " ")) > 80 || length(view_names) > 6, 8, 10)
  hm <- pheatmap(
    heatmap_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = TRUE,    # <- now show view names
    angle_col = 45,          # rotate column labels for readability
    fontsize_row = 6,
    fontsize_col = col_fontsize,
    main = title_text
  )

  invisible(list(
    heatmap = hm,
    values = heatmap_mat,
    top_genes = top_genes_all,
    matched_path = matched_path,
    enriched_views = enriched_views
  ))
}


#' Plot ligand-receptor interaction network
#'
#' Visualize ligand-receptor interactions as a directed network, with nodes
#' representing ligands or receptors in specific cell types (views), and edges
#' representing interactions. Node size reflects ligand/receptor scores, edge
#' width can reflect rank scores, and edge color distinguishes intra- vs.
#' inter-cell type interactions.
#'
#' @param lr_interactions A result from the function \code{infer_lr_interactions}, which is a data frame containing ligand-receptor interaction information.
#'   Required columns include \code{sender}, \code{receiver}, \code{ligand}, \code{receptor}.
#'   Optional columns include \code{rank_score}, \code{ligand_score}, \code{receptor_score},
#'   \code{mean_diff}, and \code{pval_adj}.
#' @param topn Integer. Number of top interactions to retain based on \code{rank_score}.
#'   If \code{NULL}, all interactions are kept.
#' @param mean_diff_cutoff Numeric. Minimum value of \code{mean_diff} to retain interactions.
#'   Interactions with smaller values are filtered out. Default is \code{NULL}.
#' @param pval_cutoff Numeric. Maximum adjusted p-value (\code{pval_adj}) to retain interactions.
#'   Interactions with larger values are filtered out. Default is \code{NULL}.
#' @param palette Optional character vector of colors to use for cell types (views). If not
#'   provided, colors are automatically assigned based on the number of cell types(views).
#' @param global_views Optional character vector specifying all cell types (views). This ensures
#'   consistent coloring of nodes across different networks.
#'
#' @return A \code{ggraph} plot object visualizing the ligand-receptor network.
#'   - Nodes are colored by view and sized by ligand/receptor scores.
#'   - Directed edges indicate interactions from ligand to receptor.
#'   - Edge color differentiates inter-cell type (red) vs intra-cell type (gray40) interactions.
#'   - Edge width reflects \code{rank_score} if present.
#'


#' @export


plot_lr_interactions <- function(
    lr_interactions,
    topn = NULL,
    mean_diff_cutoff = NULL,
    pval_cutoff = NULL,
    palette = NULL,
    global_views = NULL
) {
  # --- Determine available columns ---
  has_rank <- "rank_score" %in% colnames(lr_interactions)
  has_scores <- all(c("ligand_score", "receptor_score") %in% colnames(lr_interactions))

  # --- Filtering ---
  if (!is.null(mean_diff_cutoff) && "mean_diff" %in% colnames(lr_interactions)) {
    lr_interactions <- lr_interactions %>% filter(mean_diff >= mean_diff_cutoff)
  }
  if (!is.null(pval_cutoff) && "pval_adj" %in% colnames(lr_interactions)) {
    lr_interactions <- lr_interactions %>% filter(pval_adj <= pval_cutoff)
  }
  if (nrow(lr_interactions) == 0) stop("No interactions left after filtering!")

  # --- Top-N selection ---
  if (!is.null(topn) && nrow(lr_interactions) > topn) {
    if (has_rank) {
      lr_interactions <- lr_interactions %>% arrange(desc(rank_score)) %>% slice(1:topn)
    } else {
      lr_interactions <- lr_interactions[1:topn, ]
    }
  }

  # --- Node identifiers ---
  lr_interactions <- lr_interactions %>%
    mutate(
      ligand_node = paste0(sender, "_", ligand),
      receptor_node = paste0(receiver, "_", receptor)
    )

  # --- Build igraph object ---
  g <- igraph::graph_from_data_frame(
    d = lr_interactions %>% select(from = ligand_node, to = receptor_node),
    directed = TRUE
  )

  # --- Edge width ---
  if (has_rank) {
    E(g)$edge_width <- lr_interactions$rank_score
    edge_width_range <- c(0.2, 2)
    show_edge_legend <- TRUE
  } else {
    E(g)$edge_width <- 0.25
    edge_width_range <- NULL
    show_edge_legend <- FALSE
  }

  # --- Edge interaction type ---
  edge_ends <- igraph::ends(g, E(g))
  E(g)$interaction_type <- ifelse(
    sub("_.*$", "", edge_ends[,1]) != sub("_.*$", "", edge_ends[,2]),
    "Inter-cell type", "Intra-cell type"
  )

  # --- Vertex attributes ---
  all_nodes <- unique(c(edge_ends[,1], edge_ends[,2]))
  vertices <- data.frame(
    name = all_nodes,
    view = sub("_.*$", "", all_nodes),
    gene = sub("^.*?_", "", all_nodes),
    stringsAsFactors = FALSE
  )

  if (has_scores) {
    ligand_scores <- lr_interactions %>%
      select(node = ligand_node, ligand_score) %>%
      group_by(node) %>%
      summarise(ligand_score = mean(ligand_score, na.rm = TRUE), .groups = "drop")
    receptor_scores <- lr_interactions %>%
      select(node = receptor_node, receptor_score) %>%
      group_by(node) %>%
      summarise(receptor_score = mean(receptor_score, na.rm = TRUE), .groups = "drop")
    vertices <- vertices %>%
      left_join(ligand_scores, by = c("name" = "node")) %>%
      left_join(receptor_scores, by = c("name" = "node"))
    vertices$score <- dplyr::coalesce(vertices$ligand_score, vertices$receptor_score, 0)
  } else {
    vertices$score <- 1
  }

  V(g)$gene <- vertices$gene
  V(g)$score <- vertices$score

  # --- Safe factor assignment for view with global_views ---
  if (!is.null(global_views)) {
    # Only keep views actually present in the graph
    present_views <- intersect(global_views, vertices$view)

    # Assign factor to vertex attribute
    V(g)$view <- factor(vertices$view, levels = present_views)

    # Colors for all global_views
    n_views <- length(global_views)
    if (!is.null(palette)) {
      if (length(palette) < n_views) palette <- rep(palette, length.out = n_views)
      view_colors <- setNames(palette[seq_len(n_views)], global_views)
    } else {
      if (n_views <= 8) view_colors <- setNames(RColorBrewer::brewer.pal(n_views, "Dark2"), global_views)
      else if (n_views <= 12) view_colors <- setNames(RColorBrewer::brewer.pal(n_views, "Set3"), global_views)
      else view_colors <- setNames(viridis::viridis(n_views, option = "D"), global_views)
    }

    # Use only colors for present views
    used_colors <- view_colors[present_views]

  } else {
    present_views <- unique(vertices$view)
    V(g)$view <- factor(vertices$view, levels = present_views)
    n_views <- length(present_views)
    if (!is.null(palette)) {
      if (length(palette) < n_views) palette <- rep(palette, length.out = n_views)
      used_colors <- setNames(palette[seq_len(n_views)], present_views)
    } else {
      if (n_views <= 8) used_colors <- setNames(RColorBrewer::brewer.pal(n_views, "Dark2"), present_views)
      else if (n_views <= 12) used_colors <- setNames(RColorBrewer::brewer.pal(n_views, "Set3"), present_views)
      else used_colors <- setNames(viridis::viridis(n_views, option = "D"), present_views)
    }
  }

  # --- Edge colors ---
  edge_colors <- c("Inter-cell type" = "red", "Intra-cell type" = "gray40")

  # --- Build plot ---
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(edge_width = edge_width, edge_colour = interaction_type),
                   arrow = grid::arrow(length = grid::unit(3, 'mm')),
                   end_cap = circle(3, 'mm'),
                   alpha = 0.8) +
    geom_node_point(aes(size = score, color = V(g)$view), alpha = 0.9) +
    geom_node_text(aes(label = gene), repel = TRUE, size = 3) +
    scale_edge_colour_manual(values = edge_colors, name = "Interaction Type") +
    {if (show_edge_legend)
      scale_edge_width(range = edge_width_range, name = "Rank Score")
      else
        scale_edge_width(range = c(0.25, 0.25), guide = "none")} +
    scale_color_manual(values = used_colors, name = "Cell Types", drop = FALSE) +
    scale_size_continuous(name = "Ligand/Receptor Score", trans = "identity") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    ) +
    ggtitle(paste0("Ligand-Receptor Interaction Network (", nrow(lr_interactions), " interactions)"))

  return(p)
}





#######
#' Plot top genes for a latent factor
#'
#' Generate a heatmap of the top genes for a specific latent factor (LF) across
#' multiple cell types/views. The top genes are selected based on their latent factor loadings.
#'
#' @param coPathRes A result object from \code{mcCoordPath}
#' @param lf_index Integer. Index of the latent factor to visualize. Default is 1.
#' @param topn Integer. Number of top genes to select per view. Default is 20.
#' @param scale_rows Logical. If TRUE, rows are z-score scaled before plotting. Default is FALSE.
#' @param view_names Optional character vector specifying the names of the views. If NULL,
#'   If NULL, defaults to "View1", "View2", etc.
#' @param cluster_cols Logical. If TRUE, columns (views) will be clustered in the heatmap. Default is FALSE.
#'
#' @return A heatmap displaying the top genes (rows) across views (columns) for the selected LF.
#'   Missing values (genes not present in a view) are shown as blank cells.
#'


#' @export


plot_top_genes_LF <- function(coPathRes, lf_index = 1, topn = 20, scale_rows = F, view_names = NULL,cluster_cols=F) {
  if (is.null(coPathRes$Z_list)) stop("res must contain Z_list.")
  Z_list <- coPathRes$Z_list
  n_views <- length(Z_list)

  AUC_list <- coPathRes$AUC_list

  # Assign view names
  if (is.null(names(AUC_list)) & is.null(view_names)) {
    names(AUC_list) <- paste0("View", seq_along(AUC_list))
  }else if(!is.null(view_names) & length(AUC_list)==length(view_names)){
    names(AUC_list)<-view_names}else {warining("the length of view names does not match the number of views in the res")}


  # Step 1: Get top genes per view (Z >= 0, no abs)
  top_genes_list <- lapply(Z_list, function(Z) {
    if (lf_index > ncol(Z)) stop("lf_index exceeds number of LFs in Z.")
    Z_lf <- Z[, lf_index]
    names(sort(Z_lf, decreasing = TRUE))[1:min(topn, length(Z_lf))]
  })

  # Step 2: Union of all top genes
  top_genes_all <- unique(unlist(top_genes_list))

  # Step 3: Create gene × view matrix, fill missing genes with NA
  heatmap_mat <- do.call(cbind, lapply(Z_list, function(Z) {
    genes_in_Z <- intersect(top_genes_all, rownames(Z))
    Z_lf <- rep(NA, length(top_genes_all))
    names(Z_lf) <- top_genes_all
    Z_lf[genes_in_Z] <- Z[genes_in_Z, lf_index]
    return(Z_lf)
  }))

  rownames(heatmap_mat) <- top_genes_all
  colnames(heatmap_mat) <- names(AUC_list)

  # Step 4: Optional row scaling, skip rows with all NA
  if (scale_rows) {
    heatmap_mat <- t(apply(heatmap_mat, 1, function(x) {
      if (all(is.na(x))) return(x)
      scale(x)
    }))
  }

  # Step 5: Plot
   pheatmap(
    heatmap_mat,
    cluster_rows = TRUE,
    cluster_cols = cluster_cols,
    fontsize_row = 6,
    main = paste0("Top ", topn, " Genes per View for LF", lf_index),
    na_col = "white"  # optional: color for NA
  )
}

######
#' Plot gene expression across multiple cell types/views
#'
#' Visualize expression of selected genes across multiple cell types/views
#' using line plots or grouped boxplots depending on whether a grouping variable is provided.
#'
#' @param data_list A list of expression matrices (genes × samples) for each view.
#' @param genes Character vector of gene names to plot.
#' @param view_names Optional character vector specifying names of the views. Default is "View1", "View2", etc.
#' @param sample_group Optional vector specifying sample groups. Length must match the number of columns in each matrix.
#'   If provided, boxplots are plotted; otherwise, line plots across samples are used.
#' @param colors Optional vector of colors for views. If NULL, a default palette is used.
#' @param sample_names Optional vector of sample names for x-axis labels.
#' @param test_method Statistical test to use for pairwise comparisons when `group` is provided.
#'   Default is "t.test".
#'
#' @return A \code{ggplot} object displaying gene expression across views and optionally groups.
#'

#' @export



plot_genes_in_multiViews <- function(data_list, genes, view_names = NULL, sample_group = NULL, colors = NULL,
                                     sample_names = NULL, test_method = "t.test") {
  n_views <- length(data_list)

  if (!is.null(sample_group)) sample_group <- as.factor(sample_group)
  if (is.null(view_names)) view_names <- paste0("View", seq_len(n_views))

  # Prepare data frame
  plot_df <- lapply(seq_len(n_views), function(i) {
    data_mat <- data_list[[i]]
    common_genes <- intersect(genes, rownames(data_mat))
    if (length(common_genes) == 0) return(NULL)

    df <- data_mat[common_genes, , drop = FALSE] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene") %>%
      tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expression")

    df$view <- view_names[i]

    # Group info
    if (!is.null(sample_group)) {
      if (length(sample_group) != ncol(data_mat)) stop("Length of sample_group vector must match number of samples in each view.")
      df$group <- rep(sample_group, each = length(common_genes))
    } else {
      df$group <- "All"
    }

    # Replace sample names if provided
    if (!is.null(sample_names)) {
      df$sample <- rep(sample_names, each = length(common_genes))
    }

    # Combine view and group for boxplot x-axis
    df$view_group <- factor(
      paste0(df$view, "_", df$group),
      levels = as.vector(t(outer(view_names, unique(df$group), paste, sep = "_")))
    )

    df
  }) %>% dplyr::bind_rows()

  if (nrow(plot_df) == 0) stop("No matching genes found in data_list.")

  # Ensure view factor ordered by view_names
  plot_df$view <- factor(plot_df$view, levels = view_names)

  # Colors for views
  unique_views <- levels(plot_df$view)
  if (is.null(colors)) {
    n_colors <- length(unique_views)
    if (n_colors <= 8) colors <- RColorBrewer::brewer.pal(n_colors, "Dark2")
    else colors <- viridis::viridis(n_colors)
  }
  view_colors <- setNames(colors[seq_along(unique_views)], unique_views)

  # Map group to alpha
  if (!is.null(sample_group)) {
    group_levels <- unique(plot_df$group)
    alpha_map <- setNames(c(1, 0.4)[seq_along(group_levels)], group_levels)
  } else {
    alpha_map <- 1
  }

  # --- Plot ---
  if (is.null(sample_group) || all(plot_df$group == "All")) {
    # Line plot across samples
    plot_df <- plot_df %>%
      dplyr::group_by(gene, view) %>%
      dplyr::mutate(sample_index = seq_along(sample))

    p <- ggplot(plot_df, aes(x = sample_index, y = expression, color = view, group = view)) +
      geom_line(size = 1) +
      geom_point(size = 2) +
      facet_wrap(~gene, scales = "free_y") +
      scale_color_manual(values = view_colors) +
      theme_bw() +
      labs(x = "Sample", y = "Expression", color = "View") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right")

    # Add custom sample names if provided
    if (!is.null(sample_names)) {
      p <- p + scale_x_continuous(breaks = seq_along(sample_names),
                                  labels = sample_names)
    }

  } else {
    # Boxplot + jitter + pairwise tests
    comparisons <- combn(unique(plot_df$group), 2, simplify = FALSE)

    p <- suppressWarnings(
      ggplot(plot_df, aes(x = view_group, y = expression, fill = view, alpha = group)) +
        geom_boxplot(outlier.shape = NA, color = "black") +
        geom_jitter(aes(color = group), width = 0.15, size = 1.5, alpha = 0.7) +
        facet_wrap(~gene, scales = "free_y") +
        scale_fill_manual(values = view_colors) +
        scale_alpha_manual(values = alpha_map) +
        theme_bw() +
        labs(x = "Cell Type and Group", y = "Expression", fill = "View", color = "Group", alpha = "Group") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right")
    )

    if (length(comparisons) > 0) {
      p <- p + ggpubr::stat_compare_means(aes(x = view_group, y = expression),
                                          method = test_method,
                                          comparisons = comparisons,
                                          label = "p.signif")
    }
  }

  return(p)
}
