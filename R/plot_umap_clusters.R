#' Plot a UMAP embedding colored by gene cluster assignment
#'
#' @description
#' Creates a publication-ready UMAP scatter plot where each point represents a
#' gene and points are colored by cluster assignment. The function joins the
#' UMAP embedding coordinates to \code{cluster_table} via gene identifiers and
#' returns a \pkg{ggplot2} object suitable for saving with
#' \code{ggplot2::ggsave()} or a graphics device such as \pkg{svglite}.
#'
#' This function performs plotting only. It does not compute the UMAP embedding
#' or gene clusters. Input validation is delegated to an internal helper and
#' errors are emitted via \code{rlang::abort()}.
#'
#' @param cluster_table
#' Data frame or tibble with one row per gene and at least the columns
#' \code{gene} and \code{cluster}. The \code{gene} column must match
#' \code{rownames(umap_embedding)}.
#'
#' @param umap_embedding
#' Numeric matrix of UMAP coordinates with genes in rows and embedding
#' dimensions in columns. Must have at least two columns. Row names must be
#' gene identifiers used to join against \code{cluster_table$gene}.
#'
#' @param point_size
#' Numeric scalar giving the point size used in \code{ggplot2::geom_point()}.
#' Must be strictly positive.
#'
#' @return
#' A \pkg{ggplot2} object representing the UMAP scatter plot colored by cluster.
#'
#' @examples
#' set.seed(1)
#' genes <- paste0("gene", 1:20)
#'
#' cluster_table <- tibble::tibble(
#'   gene = genes,
#'   cluster = sample(1:3, length(genes), replace = TRUE)
#' )
#'
#' umap_embedding <- matrix(
#'   rnorm(length(genes) * 2),
#'   ncol = 2,
#'   dimnames = list(genes, c("UMAP1", "UMAP2"))
#' )
#'
#' p <- plot_umap_clusters(
#'   cluster_table  = cluster_table,
#'   umap_embedding = umap_embedding,
#'   point_size     = 1.2
#' )
#' p
#' 
#' @importFrom rlang .data
#'
#' @export
#' 
plot_umap_clusters <- function(
        cluster_table,
        umap_embedding,
        point_size = 1.2
) {
    .check_plot_umap_clusters_input(
        cluster_table  = cluster_table,
        umap_embedding = umap_embedding,
        point_size     = point_size
    )
    
    df_emb <- tibble::tibble(
        gene  = rownames(umap_embedding),
        UMAP1 = umap_embedding[, 1],
        UMAP2 = umap_embedding[, 2]
    )
    
    df <- dplyr::left_join(
        df_emb,
        cluster_table,
        by = "gene"
    )
    
    if (all(is.na(df$cluster))) {
        rlang::abort(c(
            "All clusters are NA after join.",
            "'cluster_table$gene' does not match embedding rownames."
        ))
    }
    
    df$cluster <- factor(df$cluster)
    
    ggplot2::ggplot(
        df,
        ggplot2::aes(
            x     = .data$UMAP1,
            y     = .data$UMAP2,
            color = .data$cluster
        )
    ) +
        ggplot2::geom_point(
            size = point_size,
            alpha = 0.9
        ) +
        ggplot2::labs(
            x = "UMAP 1",
            y = "UMAP 2",
            color = "Cluster"
        ) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(
            legend.position = "right",
            legend.title = ggplot2::element_text(size = 11),
            legend.text  = ggplot2::element_text(size = 10)
        )
}


# Level 1 function definitions -------------------------------------------------


#' Validate inputs for UMAP cluster plotting
#'
#' @noRd
#'
#' @description
#' Internal helper that validates the structure and contents of the inputs
#' passed to \code{plot_umap_clusters()}.
#'
#' This function verifies that:
#'
#' \itemize{
#'   \item{\code{cluster_table} is a data frame (or tibble) containing the
#'   required columns \code{gene} and \code{cluster}.}
#'   \item{\code{umap_embedding} is a numeric matrix with at least two columns
#'   and non-\code{NULL} row names.}
#'   \item{\code{point_size} is a single positive numeric value.}
#'   \item{At least one gene identifier overlaps between
#'   \code{cluster_table$gene} and \code{rownames(umap_embedding)}.}
#' }
#'
#' If only a subset of embedding genes is present in \code{cluster_table},
#' a warning is emitted via \code{rlang::warn()} to indicate that some points
#' will be plotted without a cluster assignment.
#'
#' All checks are performed at the interface level so that downstream plotting
#' code may assume well-formed inputs.
#'
#' @param cluster_table
#' Data frame with one row per gene containing at least the columns
#' \code{gene} and \code{cluster}. The \code{gene} column is used to align
#' cluster assignments to \code{umap_embedding}.
#'
#' @param umap_embedding
#' Numeric matrix of UMAP coordinates with genes in rows and embedding
#' dimensions in columns. Must have at least two columns and non-\code{NULL}
#' row names corresponding to gene identifiers.
#'
#' @param point_size
#' Numeric scalar specifying the point size used for plotting. Must be
#' strictly positive.
#'
#' @return
#' Invisibly returns \code{TRUE} if all checks pass. Otherwise, raises an
#' error via \code{rlang::abort()}.
#' 
.check_plot_umap_clusters_input <- function(
        cluster_table,
        umap_embedding,
        point_size
) {
    # cluster_table
    if (!inherits(cluster_table, "data.frame")) {
        rlang::abort(
            "`cluster_table` must be a data.frame or tibble."
        )
    }
    
    required_cols <- c("gene", "cluster")
    missing_cols <- setdiff(required_cols, colnames(cluster_table))
    if (length(missing_cols) > 0) {
        rlang::abort(c(
            "`cluster_table` is missing required columns.",
            paste0("Missing: ", paste(missing_cols, collapse = ", "))
        ))
    }
    
    if (anyNA(cluster_table$gene)) {
        rlang::abort("`cluster_table$gene` contains NA values.")
    }
    
    # umap_embedding
    if (!is.matrix(umap_embedding)) {
        rlang::abort("`umap_embedding` must be a numeric matrix.")
    }
    
    if (!is.numeric(umap_embedding)) {
        rlang::abort("`umap_embedding` must be numeric.")
    }
    
    if (ncol(umap_embedding) < 2) {
        rlang::abort("`umap_embedding` must have at least 2 columns.")
    }
    
    if (is.null(rownames(umap_embedding))) {
        rlang::abort(
            "`umap_embedding` must have rownames corresponding to gene IDs."
        )
    }
    
    # point_size 
    if (!is.numeric(point_size) 
        || length(point_size) != 1 
        || point_size <= 0) {
        rlang::abort("`point_size` must be a single positive numeric value.")
    }
    
    # alignment diagnostics 
    overlap <- sum(rownames(umap_embedding) %in% cluster_table$gene)
    
    if (overlap == 0) {
        rlang::abort(c(
            "No overlap between `umap_embedding` rownames and ",
            "`cluster_table$gene`.",
            "Check that both use the same gene identifiers."
        ))
    }
    
    if (overlap < nrow(umap_embedding)) {
        rlang::warn(c(
            "Some UMAP points have no cluster assignment.",
            paste0(
                "Matched genes: ", overlap, " / ", nrow(umap_embedding)
            )
        ))
    }

    invisible(TRUE)
}
