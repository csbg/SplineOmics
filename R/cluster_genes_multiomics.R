#' Gene-Centric Multi-Omics Clustering Across Blocks and Layers
#'
#' @description
#' Performs gene-centric clustering of multi-omics time-series data across
#' multiple blocks (e.g., time effect, interaction effect) and data modalities
#' (e.g., transcript, protein, feature-level modalities with many-to-one gene
#' mapping). The function integrates multiple modalities within each block by
#' computing modality-specific gene–gene distances, combining them via 
#' user-defined weights, and clustering genes based on the resulting unified
#'  distance matrix.
#'
#' The function is flexible with respect to how spline trajectories,
#' interaction representations, or feature-level signatures are constructed:
#' these are precomputed outside the function, and supplied as matrices inside
#' the \code{blocks} structure. Internally, the function simply performs
#' harmonization of genes, distance computation, weighted integration, and
#' clustering.
#'
#' @param blocks
#' A named nested list specifying all data used for clustering.
#' The outer list corresponds to analytical \emph{blocks} (e.g.,
#' \code{time_Ctrl}, \code{interaction_Ctrl_vs_Treat}).
#' Each element of the outer list is itself a named list whose elements are
#' \emph{modalities} (e.g., \code{rna}, \code{protein}, \code{phospho}), each
#'  being a numeric matrix of dimension \code{features x spline_points}.
#'
#' For one-to-one (gene-level) modalities, rows represent genes directly.
#' In this case, row names must be the gene identifiers themselves and must
#' follow the pattern \code{<gene_id>}. The angle brackets are shown for
#' illustration only and must not be included in the actual row names.
#' Gene identifiers must be consistent across all one-to-one layers; otherwise,
#' genes cannot be matched across omics layers during distance computation and
#' clustering.
#'
#' For many-to-one modalities (e.g., phospho sites, probes), rows represent
#' features that map to genes and are summarized into gene-level pattern
#' signatures
#' based on the metadata tables. For these modalities, row names must follow the
#' pattern \code{<gene_id>_<feature_id>}, where the gene identifier precedes
#' the first underscore. Again, the angle brackets are for illustration only
#' and must not be included in the actual row names.
#'
#' This row-naming convention is critical, as it defines how features are
#' associated with genes and how genes are aligned across modalities prior to
#' signature construction and downstream clustering.
#'
#' @param block_meta
#' A data frame containing \emph{block-level metadata}. One row per block.  
#' Must include:
#' \describe{
#'   \item{\code{block}}{Block identifier (must match names in \code{blocks}).}
#'   \item{\code{block_k}}{Number of gene clusters for this block.}
#'   \item{\code{result_category}}{Numeric result category label (e.g., 1 for
#'   time effect, 3 for interaction).}
#'   \item{\code{cond1}}{Primary condition associated with this block (e.g.,
#'   condition for time-effect blocks, or first condition in a contrast).}
#'   \item{\code{cond2}}{Secondary condition for contrast-based blocks. Set to
#'   \code{NA} for single-condition blocks.}
#' }
#'
#' @param modality_meta
#' A data frame containing \emph{modality-level metadata}. One row per (block ×
#' modality). Must include:
#' \describe{
#'   \item{\code{block}}{Block identifier linking to \code{block_meta}.}
#'   \item{\code{layer}}{Modality name within the block.}
#'   \item{\code{layer_k}}{Number of pattern clusters to use for building
#'   pattern signatures for many-to-one modalities. \code{NA} for modalities
#'   that are already gene-level.}
#'   \item{\code{layer_w}}{
#'     Relative weight of this modality when combining modality-specific
#'     gene–gene distances within the block. Values are treated as
#'     *relative* weights and are normalized internally, so they do not
#'     need to sum to 1 (e.g. 1, 1, 2 means the third modality has twice
#'     the weight of the others).
#'   }
#' }
#'
#' @param gene_mode
#' Character string specifying how genes should be harmonized across modalities
#' within each block prior to clustering.  
#' \describe{
#'   \item{\code{"intersection"}}{Retain only genes present in \emph{all}
#'   modalities of the block. Produces the most interpretable multi-omics
#'    clusters.}
#'
#'   \item{\code{"union"}}{Retain genes present in \emph{any} modality of the
#'   block. Gene–gene distances are computed using only shared modalities per
#'   gene pair, with weights renormalized accordingly. Increases coverage but
#'   results in heterogeneous information across genes.}
#' }
#' 
#' @param verbose Boolean flag indicating if info messages are be shown.
#'
#' @return
#' A named list with two tibbles:
#'
#' \describe{
#'   \item{\code{cluster_table}}{A tibble with one row per gene containing
#'   block-specific cluster assignments suitable for downstream enrichment
#'   analyses. Columns include the gene identifier and one clustering column
#'   per analytical block (e.g., \code{cluster_<cond>} or
#'   \code{cluster_cat3_<cond1>_vs_<cond2>}). Genes not included in a given
#'   block are assigned \code{NA}.}
#'
#'   \item{\code{centroid_info}}{A tibble with one row per block, modality,
#'   and cluster, summarizing modality-specific cluster centroid
#'   trajectories and within-cluster coherence. Columns include the block
#'   and layer identifiers, cluster label, gene coverage statistics, mean
#'   and standard deviation of per-gene R^2 values, optional per-gene
#'   R^2 vectors, and the centroid trajectory stored as a list-column.}
#' }
#' 
#' @examples
#' set.seed(1)
#'
#' genes <- paste0("gene", 1:6)
#'
#' rna_time_ctrl <- matrix(
#'     rnorm(6 * 5),
#'     nrow = 6,
#'     ncol = 5,
#'     dimnames = list(genes, NULL)
#' )
#' rna_time_treat <- matrix(
#'     rnorm(6 * 5),
#'     nrow = 6,
#'     ncol = 5,
#'     dimnames = list(genes, NULL)
#' )
#'
#' blocks <- list(
#'     time_Ctrl = list(rna = rna_time_ctrl),
#'     time_Treat = list(rna = rna_time_treat)
#' )
#'
#' block_meta <- data.frame(
#'     block           = c("time_Ctrl", "time_Treat"),
#'     block_k         = c(2L, 2L),
#'     result_category = c(1, 1),
#'     cond1           = c("Ctrl", "Treat"),
#'     cond2           = c(NA_character_, NA_character_),
#'     stringsAsFactors = FALSE
#' )
#'
#' modality_meta <- data.frame(
#'     block   = c("time_Ctrl", "time_Treat"),
#'     layer   = c("rna", "rna"),
#'     layer_k = c(NA_real_, NA_real_),
#'     layer_w = c(1, 1),
#'     stringsAsFactors = FALSE
#' )
#'
#' cluster_table <- cluster_genes_multiomics(
#'     blocks     = blocks,
#'     block_meta = block_meta,
#'     modality_meta = modality_meta,
#'     gene_mode  = "intersection"
#' )
#'
#' cluster_table
#'
#' @export
#' 
cluster_genes_multiomics <- function(
        blocks,
        block_meta,
        modality_meta,
        gene_mode = "intersection",
        verbose = TRUE
) {
    start_time <- Sys.time()
    .check_cluster_genes_multiomics_input(
        blocks     = blocks,
        block_meta = block_meta,
        modality_meta = modality_meta,
        gene_mode  = gene_mode,
        verbose = verbose
    )

    block_ids <- unique(block_meta$block)
    block_clusters <- vector("list", length(block_ids))
    names(block_clusters) <- block_ids
    
    centroid_list <- vector("list", length(block_ids))
    names(centroid_list) <- block_ids
    
    for (i in seq_along(block_ids)) {
        b <- block_ids[i]
        if (verbose) {
            message(
                "[cluster_genes_multiomics] Block ",
                i,
                "/",
                length(block_ids),
                ": '",
                b,
                "'"
            )
        }
        res_b <- .cluster_genes_multiomics_block(
            block_id   = b,
            blocks     = blocks,
            block_meta = block_meta,
            modality_meta = modality_meta,
            gene_mode  = gene_mode,
            verbose = verbose
        )
        
        block_clusters[[b]] <- res_b$cl_b
        centroid_list[[b]]  <- res_b$centroid_info
    }
    
    genes_all <- sort(
        unique(unlist(lapply(block_clusters, names)))
    )
    cluster_table <- data.frame(
        gene = genes_all,
        stringsAsFactors = FALSE
    )
    
    for (i in seq_along(block_ids)) {
        b <- block_ids[i]
        cl_b <- block_clusters[[b]]
        if (is.null(cl_b)) next
        
        mb <- block_meta[block_meta$block == b, , drop = FALSE][1L, ]
        
        rc <- mb$result_category
        c1 <- mb$cond1
        c2 <- mb$cond2
        
        if (rc == 1) {
            col_name <- paste0("cluster_", c1)
        } else if (rc == 3) {
            col_name <- paste0(
                "cluster_cat3_",
                c1, "_vs_", c2
            )
        } else {
            next
        }
        
        if (!col_name %in% colnames(cluster_table)) {
            cluster_table[[col_name]] <- NA_integer_
        }
        
        genes_b <- names(cl_b)
        idx_tbl <- match(genes_b, cluster_table$gene)
        cluster_table[[col_name]][idx_tbl] <- as.integer(cl_b)
    }
    
    centroid_info <- do.call(
        rbind,
        centroid_list
    )
    
    if (verbose) {
        end_time <- Sys.time()
        elapsed  <- end_time - start_time
        formatted <- format(elapsed, digits = 2)
        message("[cluster_genes_multiomics] total runtime: ", formatted)
    }
    list(
        cluster_table = tibble::as_tibble(cluster_table),
        centroid_info = tibble::as_tibble(centroid_info)
    )
}


# Level 1 function definitions -------------------------------------------------


#' Validate inputs for multi-omics gene clustering
#'
#' @noRd
#'
#' @description
#' Internal helper that validates all inputs passed to
#' \code{cluster_genes_multiomics()}.  
#'
#' It performs:
#' \itemize{
#'   \item argument-wise checks for \code{blocks}, \code{block_meta},
#'         \code{modality_meta}, and \code{gene_mode}, and
#'   \item cross-argument consistency checks linking blocks, layers,
#'         and metadata.
#' }
#'
#' This function is a thin orchestrator that delegates to specialized
#' internal checkers:
#' \code{.check_cluster_genes_multiomics_gene_mode()},
#' \code{.check_cluster_genes_multiomics_blocks()},
#' \code{.check_cluster_genes_multiomics_block_meta()},
#' \code{.check_cluster_genes_multiomics_modality_meta()}, and
#' \code{.check_cluster_genes_multiomics_cross_args()}.
#'
#' @param blocks
#' Named list of blocks, one element per block.  
#' Each block is itself a named list of layers, where each layer is a
#' numeric matrix (\code{features x spline_points}) with non-empty
#' row names (feature or gene IDs).
#'
#' @param block_meta
#' Data frame with block-level metadata.  
#' Must contain the columns \code{block}, \code{block_k},
#' \code{result_category}, \code{cond1}, and \code{cond2}.  
#' The column \code{block_k} gives the number of clusters per block and
#' must be positive and constant within each block.
#'
#' @param modality_meta
#' Data frame with layer-level metadata.  
#' Must contain the columns \code{block}, \code{layer}, \code{layer_k},
#' and \code{layer_w}.  
#' The column \code{layer_k} encodes whether a layer is many-to-one
#' (positive integer) or one-to-one (NA).  
#' The column \code{layer_w} provides non-negative layer weights whose
#' sum must be positive within each block.
#'
#' @param gene_mode
#' Character scalar specifying how genes should be combined across
#' layers within a block.  
#' Must be either \code{"intersection"} or \code{"union"} and is
#' validated via \code{\link[base]{match.arg}}.
#' 
#' @param verbose Boolean flag indicating if info messages are be shown.
#'
#' @return
#' Invisibly returns \code{TRUE} if all checks pass.  
#' Otherwise, it raises an informative error (or warning) describing the
#' first detected issue.
#' 
.check_cluster_genes_multiomics_input <- function(
        blocks,
        block_meta,
        modality_meta,
        gene_mode,
        verbose
) {
    .check_cluster_genes_multiomics_gene_mode(gene_mode)
    .check_cluster_genes_multiomics_blocks(blocks)
    .check_cluster_genes_multiomics_block_meta(block_meta)
    .check_cluster_genes_multiomics_modality_meta(modality_meta)
    .check_cluster_genes_multiomics_cross_args(
        blocks,
        block_meta,
        modality_meta
        )
    
    # verbose must be TRUE or FALSE (length 1 logical)
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop_call_false(
            "`verbose` must be either TRUE or FALSE (a single logical value)."
        )
    }
    
    invisible(TRUE)
}


#' Cluster genes within a single block and compute layer-wise centroids
#'
#' @noRd
#'
#' @description
#' Internal helper that performs all clustering steps for a single block
#' in the multi-omics workflow.  
#'
#' For the specified block, the function:
#' \itemize{
#'   \item extracts the relevant block-level and layer-level metadata,
#'   \item converts many-to-one layers into gene-level signatures,
#'   \item selects genes according to the requested \code{gene_mode}
#'         (\code{"intersection"} or \code{"union"}),
#'   \item constructs a block-level distance matrix by combining
#'         layer-level distances using per-layer weights, and
#'   \item runs the clustering algorithm to obtain block-specific
#'         gene–cluster assignments.
#' }
#'
#' After clustering, it computes per-layer cluster centroids and quality
#' metrics through \code{.compute_block_centroids()}.
#'
#' @param block_id
#' Identifier of the block being processed. Must match a value in
#' \code{block_meta$block} and \code{names(blocks)}.
#'
#' @param blocks
#' Named list of blocks as supplied to
#' \code{cluster_genes_multiomics()}.  
#' Each block is a named list of raw layer matrices.
#'
#' @param block_meta
#' Block-level metadata table.  
#' Used to obtain \code{block_k} and to filter metadata for the current
#' block.
#'
#' @param modality_meta
#' Layer-level metadata table.  
#' Used to retrieve layer names, layer types (\code{layer_k}), and
#' layer weights for the current block.
#'
#' @param gene_mode
#' Character scalar specifying how to combine genes across layers:
#' either \code{"intersection"} or \code{"union"}.  
#' Determines whether only shared genes or all genes across layers
#' contribute to the block-level distance.
#'
#' @return
#' A list with two components:
#' \describe{
#'   \item{\code{cl_b}}{
#'     Named vector of cluster assignments for all genes in the block.
#'   }
#'
#'   \item{\code{centroid_info}}{
#'     Data frame with per-layer, per-cluster centroid trajectories and
#'     associated summary metrics, as produced by
#'     \code{.compute_block_centroids()}.
#'   }
#' }
#'
.cluster_genes_multiomics_block <- function(
        block_id,
        blocks,
        block_meta,
        modality_meta,
        gene_mode,
        verbose
) {
    meta_b_block <- block_meta[block_meta$block == block_id, , drop = FALSE]
    meta_b_layer <- modality_meta[modality_meta$block == block_id, , drop = FALSE]
    
    block_k <- meta_b_block$block_k[1L]
    
    layer_names   <- meta_b_layer$layer
    layer_weights <- meta_b_layer$layer_w
    
    layer_mats_raw <- lapply(
        layer_names,
        function(ln) blocks[[block_id]][[ln]]
    )
    names(layer_mats_raw) <- layer_names
    
    layer_mats <- vector("list", length(layer_names))
    names(layer_mats) <- layer_names
    
    if (verbose) {
        message(
            "  [block '",
            block_id,
            "'] building gene-level layer matrices..."
            )
    }
    
    for (i in seq_along(layer_names)) {
        ln      <- layer_names[i]
        lk      <- meta_b_layer$layer_k[i]
        mat_raw <- layer_mats_raw[[ln]]
        
        if (is.na(lk)) {
            layer_mats[[ln]] <- mat_raw
        } else {
            rn      <- rownames(mat_raw)
            gene_ids <- sub("_.*$", "", rn)
            feature_to_gene <- gene_ids
            
            sig <- .build_site_signatures(
                layer_mat       = mat_raw,
                feature_to_gene = feature_to_gene,
                layer_k         = lk
            )$signatures
            
            sig <- sig[order(rownames(sig)), , drop = FALSE]
            layer_mats[[ln]] <- sig
        }
    }
    
    genes_per_layer <- lapply(layer_mats, rownames)
    
    if (verbose) {
        message(
            "  [block '",
            block_id,
            "'] computing and combining layer-wise distance matrices..."
            )
    }
    
    if (gene_mode == "intersection") {
        D_block <- .compute_block_distance_intersection(
            layer_mats    = layer_mats,
            layer_names   = layer_names,
            layer_weights = layer_weights
        )
    } else {    # gene_mode == "union"
        D_block <- .compute_block_distance_union(
            layer_mats    = layer_mats,
            layer_names   = layer_names,
            layer_weights = layer_weights
        )
    }
    if (verbose) {
        message(
            "  [block '",
            block_id,
            "'] clustering genes (k = ",
            block_k,
            ")..."
            )
    }
    cl_b <- .cluster_with_kmeans(
        dist_mat = D_block,
        k        = block_k
    )
    
    centroid_info_b <- .compute_block_centroids(
        block_id     = block_id,
        layer_mats   = layer_mats,
        cl_b         = cl_b,
        meta_b_layer = meta_b_layer
    )
    
    list(
        cl_b          = cl_b,
        centroid_info = centroid_info_b
    )
}


# Level 2 function definitions -------------------------------------------------


#' Build gene-level pattern signatures from many-to-one features
#'
#' @noRd
#'
#' @description
#' Internal helper that converts a many-to-one omics layer (e.g. sites,
#' probes) into **gene-level pattern signatures**.  
#'
#' All features in the layer are first clustered *jointly* into
#' \code{layer_k} global dynamic archetypes based on their temporal
#' trajectories. These global archetypes represent recurring temporal
#' shapes observed across the entire layer.
#'
#' Each gene is then represented by a fixed-length signature vector that
#' quantifies how its own features distribute across these shared
#' archetypes—effectively a mixture over global temporal patterns.
#'
#' @param layer_mat
#' Numeric matrix (\code{features x spline_points}) where each row is a
#' feature-level trajectory (e.g. phosphosite, probe).  
#' Row names must uniquely identify features; no specific format is
#' required.
#'
#' @param feature_to_gene
#' Character vector of length \code{nrow(layer_mat)} mapping each feature
#' (row of \code{layer_mat}) to its corresponding gene.  
#' The mapping is **strictly positional**:  
#' the \(i\)-th entry of \code{feature_to_gene} corresponds to the
#' \(i\)-th row of \code{layer_mat}.  
#' Names are ignored.
#'
#' @param layer_k
#' Integer scalar giving the number of global archetype clusters to
#' compute across all features.  
#' Must satisfy \code{2 <= layer_k <= nrow(layer_mat)}.
#'
#' @param seed
#' Optional integer seed passed to \code{\link[base]{set.seed}} to make
#' the feature-level clustering reproducible.
#'
#' @return
#' A list with two components:
#'
#' \describe{
#'
#'   \item{\code{signatures}}{
#'     A numeric matrix of dimension \code{genes x layer_k}.  
#'     Each row represents a gene; each column corresponds to one global
#'     temporal archetype.  
#'
#'     Values are the **proportion of that gene's features** assigned to
#'     each archetype (rows sum to 1 for genes that have at least one
#'     mapped feature).
#'   }
#'
#'   \item{\code{feature_clusters}}{
#'     Integer vector of length \code{nrow(layer_mat)} giving the
#'     archetype assignment for each feature.  
#'     Names correspond to feature IDs (rownames of
#'     \code{layer_mat}).
#'   }
#'
#' }
#'
#' @importFrom stats dist
#' 
.build_site_signatures <- function(
        layer_mat,
        feature_to_gene,
        layer_k,
        seed = NULL
) {
    layer_k <- as.integer(layer_k)
    
    # positional mapping: i-th row of layer_mat -> feature_to_gene[i]
    feat_map <- as.character(feature_to_gene)
    
    if (length(feat_map) != nrow(layer_mat)) {
        stop(
            "`feature_to_gene` must have length equal to nrow(layer_mat). ",
            "Mapping is positional: element i is the gene ID for row i.",
            call. = FALSE
        )
    }
    
    # 1) Global feature-by-feature distance matrix (all sites together)
    D_feat <- .compute_layer_distance(layer_mat)
    
    # 2) Cluster all features (sites) into layer_k global shape archetypes
    feature_clusters <- .cluster_with_kmeans(
        dist_mat = D_feat,
        k        = layer_k,
        seed     = seed
    )
    names(feature_clusters) <- rownames(layer_mat)
    
    # 3) Build gene x K signature matrix (fractions per global archetype)
    genes <- sort(unique(feat_map))
    sig_mat <- matrix(
        0,
        nrow = length(genes),
        ncol = layer_k
    )
    rownames(sig_mat) <- genes
    colnames(sig_mat) <- paste0("cluster_", seq_len(layer_k))
    
    for (g in genes) {
        idx <- which(feat_map == g)
        if (length(idx) == 0L) {
            next
        }
        
        clust_g <- feature_clusters[idx]
        tab <- table(clust_g)
        sig_mat[g, as.integer(names(tab))] <- as.numeric(tab)
        
        row_sum <- sum(sig_mat[g, ])
        if (row_sum > 0) {
            sig_mat[g, ] <- sig_mat[g, ] / row_sum
        }
    }
    
    list(
        signatures       = sig_mat,        # genes x layer_k, fractions
        feature_clusters = feature_clusters # site -> global pattern
    )
}


#' Compute block-level distance matrix (intersection mode)
#'
#' @noRd
#'
#' @description
#' Internal helper that computes a block-level gene-gene distance matrix
#' in \code{gene_mode = "intersection"}.  
#'
#' Only genes present in \emph{all} layers of the block are retained.
#' For this common gene set, the function:
#' \itemize{
#'   \item computes a layer-specific distance matrix for each layer, and
#'   \item combines these matrices via a weighted sum using
#'         \code{layer_weights}.
#' }
#'
#' The result is a single symmetric distance matrix whose row and column
#' names correspond to the intersected gene set.
#'
#' @param layer_mats
#' Named list of gene-level matrices, one per layer, with matching layer
#' names in \code{layer_names}. Row names are gene identifiers.
#'
#' @param layer_names
#' Character vector of layer names to use when extracting from
#' \code{layer_mats}.
#'
#' @param layer_weights
#' Numeric vector of layer weights, same length and order as
#' \code{layer_names}. Used as relative weights when combining
#' layer-specific distance matrices.
#'
#' @return
#' A numeric symmetric matrix \code{D_block} containing gene-gene
#' distances for the intersection of genes across all layers.  
#' Row and column names are the retained gene identifiers.
#'
.compute_block_distance_intersection <- function(
        layer_mats,
        layer_names,
        layer_weights
) {
    genes_per_layer <- lapply(layer_mats, rownames)
    genes_block <- Reduce(intersect, genes_per_layer)
    genes_block <- sort(genes_block)
    
    dist_list <- lapply(
        layer_names,
        function(ln) {
            mat <- layer_mats[[ln]][genes_block, , drop = FALSE]
            .compute_layer_distance(mat)
        }
    )
    names(dist_list) <- layer_names
    
    D_block <- .combine_layer_distances(
        dist_list = dist_list,
        weights   = layer_weights
    )
    
    D_block
}


#' Compute block-level distance matrix (union mode)
#'
#' @noRd
#'
#' @description
#' Internal helper that computes a block-level gene-gene distance matrix
#' in \code{gene_mode = "union"}.  
#'
#' All genes present in \emph{any} layer of the block are retained.
#' Distances are combined in a pairwise, coverage-aware fashion:
#' for each gene pair \code{(gi, gj)}, the function:
#' \itemize{
#'   \item considers only layers where \emph{both} genes are present,
#'   \item computes the layer-specific distance for that pair,
#'   \item accumulates \code{weight * distance} in a numerator matrix,
#'   \item accumulates \code{weight} in a denominator matrix, and
#'   \item forms a weighted average over all contributing layers.
#' }
#'
#' Pairs of genes that never co-occur in any layer effectively receive no
#' distance information and remain unused in practice.
#'
#' @param layer_mats
#' Named list of gene-level matrices, one per layer, with matching layer
#' names in \code{layer_names}. Row names are gene identifiers.
#'
#' @param layer_names
#' Character vector of layer names to use when extracting from
#' \code{layer_mats}.
#'
#' @param layer_weights
#' Numeric vector of layer weights, same length and order as
#' \code{layer_names}. Used as relative weights when combining
#' layer-specific distances.
#'
#' @return
#' A numeric symmetric matrix \code{D_block} containing gene-gene
#' distances for the union of genes across all layers.  
#' Each entry corresponds to a pairwise weighted average of distances
#' over the layers where both genes are observed.
#'
.compute_block_distance_union <- function(
        layer_mats,
        layer_names,
        layer_weights
) {
    genes_per_layer <- lapply(layer_mats, rownames)
    genes_block <- Reduce(union, genes_per_layer)
    genes_block <- sort(genes_block)
    
    n_g <- length(genes_block)
    D_block <- matrix(0, nrow = n_g, ncol = n_g)
    W_block <- matrix(0, nrow = n_g, ncol = n_g)
    rownames(D_block) <- genes_block
    colnames(D_block) <- genes_block
    rownames(W_block) <- genes_block
    colnames(W_block) <- genes_block
    
    for (i in seq_along(layer_names)) {
        ln <- layer_names[i]
        w  <- layer_weights[i]
        
        mat_l   <- layer_mats[[ln]]
        genes_l <- rownames(mat_l)
        
        idx <- match(genes_l, genes_block)
        
        D_l <- .compute_layer_distance(mat_l)
        
        D_block[idx, idx] <- D_block[idx, idx] + w * D_l
        W_block[idx, idx] <- W_block[idx, idx] + w
    }
    
    valid <- W_block > 0
    D_block[valid] <- D_block[valid] / W_block[valid]
    diag(D_block) <- 0
    
    D_block
}


#' Cluster genes using k-means based on a distance matrix
#'
#' @noRd
#'
#' @description
#' Internal helper that performs k-means clustering of genes based on a
#' precomputed pairwise distance matrix.  
#' The distance matrix is first embedded into a low-dimensional Euclidean
#' space using classical multidimensional scaling (MDS), after which
#' k-means (or MiniBatch k-means for large datasets) is applied.  
#' If multiple values of \code{k} are supplied, the best solution is
#' selected using the Bayesian Information Criterion (BIC).
#'
#' @param dist_mat
#' Square numeric matrix of pairwise distances (\code{genes x genes}).
#' Row and column names must be gene identifiers, which are propagated to
#' the returned cluster vector.
#'
#' @param k
#' Integer or integer vector specifying the number(s) of clusters to
#' evaluate.  
#' - If a single value is given, clustering is performed for this \code{k}.  
#' - If multiple values are supplied, each \code{k} is evaluated and the
#'   best solution is chosen via BIC.  
#' All values must satisfy \code{2 <= k < nrow(dist_mat)}.
#'
#' @param seed
#' Optional integer seed passed to \code{\link[base]{set.seed}} before
#' clustering. If \code{NULL}, no seed is set.
#'
#' @return
#' An integer vector of cluster assignments, one entry per gene.
#' Names correspond to the row names of \code{dist_mat}.
#'
#' @importFrom stats kmeans cmdscale dist
#' @importFrom ClusterR MiniBatchKmeans predict_KMeans
#' 
.cluster_with_kmeans <- function(
        dist_mat,
        k,
        seed = NULL
) {
    n_obs <- nrow(dist_mat)
    if (n_obs < 2L) {
        stop("Need at least 2 observations for clustering.", call. = FALSE)
    }
    
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    k_vec <- sort(unique(as.integer(k)))
    k_vec <- k_vec[k_vec > 1L & k_vec < n_obs]
    if (length(k_vec) == 0L) {
        stop(
            "`k` must contain values between 2 and n_obs - 1.",
            call. = FALSE
        )
    }
    
    # embed distance matrix into Euclidean space (MDS)
    k_mds <- min(10L, n_obs - 1L)
    curve_values <- stats::cmdscale(
        stats::as.dist(dist_mat),
        k = k_mds
    )
    
    use_mb <- n_obs > 1000L
    
    run_one_k <- function(kk) {
        if (!use_mb) {
            fit <- stats::kmeans(
                curve_values,
                centers  = kk,
                nstart   = 10L,
                iter.max = 300L
            )
            list(
                tot_within = fit$tot.withinss,
                cluster    = fit$cluster
            )
        } else {
            batch_size <- min(
                n_obs,
                max(20L, 2L * kk, floor(0.05 * n_obs))
            )
            fit <- ClusterR::MiniBatchKmeans(
                data            = curve_values,
                clusters        = kk,
                batch_size      = batch_size,
                num_init        = 10L,
                max_iters       = 300L,
                init_fraction   = 1.0,
                early_stop_iter = 10L,
                tol             = 1e-4,
                verbose         = FALSE,
                seed            = if (is.null(seed)) 42L else seed
            )
            list(
                tot_within = sum(fit$WCSS_per_cluster),
                cluster    = ClusterR::predict_KMeans(
                    data      = curve_values,
                    CENTROIDS = fit$centroids
                )
            )
        }
    }
    
    if (length(k_vec) == 1L) {
        res <- run_one_k(k_vec)
        cl <- res$cluster
        names(cl) <- rownames(dist_mat)
        return(cl)
    }
    
    fits <- lapply(k_vec, run_one_k)
    tot_within <- vapply(
        fits,
        function(x) x$tot_within,
        numeric(1)
    )
    
    p <- ncol(curve_values)
    bic <- n_obs * log(tot_within / n_obs) +
        k_vec * log(n_obs) * p
    best_idx <- which.min(bic)
    
    cl <- fits[[best_idx]]$cluster
    names(cl) <- rownames(dist_mat)
    
    cl
}


#' Compute per-layer cluster centroids within a block
#'
#' @noRd
#'
#' @description
#' Internal helper that summarizes gene-level clusters for a single
#' block by computing layer-specific centroid trajectories and quality
#' metrics.  
#'
#' For each layer and each cluster, the function:
#' \itemize{
#'   \item identifies the genes assigned to that cluster,
#'   \item restricts to those genes present in the layer matrix,
#'   \item z-scores each gene's trajectory (row-wise) to obtain
#'         shape-centric profiles,
#'   \item computes the centroid as the mean z-scored trajectory across
#'         genes, and
#'   \item evaluates per-gene agreement with the centroid via
#'         squared correlation (\eqn{R^2}), summarizing with the mean and
#'         standard deviation.
#' }
#'
#' Layer types are inferred from \code{meta_b_layer$layer_k}:
#' \code{NA} indicates a one-to-one gene-level layer; positive values
#' indicate many-to-one layers.
#'
#' @param block_id
#' Character (or factor) scalar giving the identifier of the block for
#' which centroids are computed.
#'
#' @param layer_mats
#' Named list of numeric matrices, one per layer in the block.  
#' Each matrix has dimensions \code{genes x spline_points}, with row
#' names corresponding to gene identifiers.
#'
#' @param cl_b
#' Vector of cluster assignments for genes in the block.  
#' Names must be gene IDs; values are cluster labels (e.g. integers).
#'
#' @param meta_b_layer
#' Data frame with layer-level metadata for this block.  
#' Must contain at least the columns \code{layer} (matching
#' \code{names(layer_mats)}) and \code{layer_k} for layer-type
#' inference.
#'
#' @return
#' A data frame with one row per \code{(layer, cluster)} combination in
#' the block and the following columns:
#' \describe{
#'   \item{\code{block}}{Block identifier (same as \code{block_id}).}
#'   \item{\code{layer}}{Layer name.}
#'   \item{\code{layer_type}}{Character label, either
#'         \code{"one_to_one"} or \code{"many_to_one"}.}
#'   \item{\code{cluster}}{Cluster label.}
#'   \item{\code{n_genes_cluster}}{Number of genes in the cluster
#'         (according to \code{cl_b}).}
#'   \item{\code{n_genes_used}}{Number of cluster genes present in the
#'         layer matrix and used to compute the centroid.}
#'   \item{\code{coverage}}{Fraction of cluster genes represented in the
#'         layer (\code{n_genes_used / n_genes_cluster}).}
#'   \item{\code{mean_R2}}{Mean squared correlation between each used
#'         gene's z-scored trajectory and the centroid.}
#'   \item{\code{sd_R2}}{Standard deviation of the \eqn{R^2} values
#'         across genes.}
#'   \item{\code{r2_member}}{List-column; each entry is a named numeric
#'         vector of per-gene \eqn{R^2} values (names are gene IDs) for
#'         genes used to compute the centroid in this layer and cluster.}
#'   \item{\code{centroid}}{List-column; each entry is a numeric vector
#'         giving the centroid trajectory (length
#'         \code{ncol(layer_mats[[layer]])}).}
#' }
#'
.compute_block_centroids <- function(
        block_id,
        layer_mats,
        cl_b,
        meta_b_layer
) {
    centroid_rows <- list()
    row_i <- 0L
    
    clusters <- sort(unique(cl_b))
    
    for (ln in names(layer_mats)) {
        mat_l <- layer_mats[[ln]]
        
        # determine layer type from meta_b_layer (NA = one-to-one)
        lk <- meta_b_layer$layer_k[meta_b_layer$layer == ln][1L]
        layer_type <- if (is.na(lk)) "one_to_one" else "many_to_one"
        
        for (c in clusters) {
            genes_cluster <- names(cl_b)[cl_b == c]
            n_cluster <- length(genes_cluster)
            
            genes_layer <- intersect(genes_cluster, rownames(mat_l))
            n_used <- length(genes_layer)
            
            coverage <- if (n_cluster > 0L) {
                n_used / n_cluster
            } else {
                NA_real_
            }
            
            if (n_used == 0L) {
                centroid <- rep(NA_real_, ncol(mat_l))
                r2_member  <- setNames(numeric(0), character(0))
                mean_R2 <- NA_real_
                sd_R2 <- NA_real_
            } else {
                X <- mat_l[genes_layer, , drop = FALSE]
                
                # row-wise z-score (shape-centric)
                X_z <- t(scale(t(X)))
                
                if (n_used == 1L) {
                    centroid <- as.numeric(X_z[1L, ])
                    r2_member <- setNames(1, genes_layer)
                    mean_R2 <- 1
                    sd_R2 <- 0
                } else {
                    centroid <- colMeans(X_z, na.rm = TRUE)
                    
                    r2_vec <- apply(
                        X_z,
                        1L,
                        function(x) {
                            r <- suppressWarnings(
                                stats::cor(x, centroid)
                            )
                            r^2
                        }
                    )
                    r2_member <- r2_vec
                    names(r2_member) <- rownames(X_z)
                    mean_R2 <- mean(r2_vec, na.rm = TRUE)
                    sd_R2 <- stats::sd(r2_vec, na.rm = TRUE)
                }
            }
            
            row_i <- row_i + 1L
            centroid_rows[[row_i]] <- data.frame(
                block           = block_id,
                modality        = ln,
                modality_type   = layer_type,
                cluster         = c,
                n_genes_cluster = n_cluster,
                n_genes_used    = n_used,
                coverage        = coverage,
                mean_R2         = mean_R2,
                sd_R2           = sd_R2,
                r2_member       = I(list(r2_member)),
                centroid        = I(list(centroid)),
                stringsAsFactors = FALSE
            )
        }
    }
    
    if (length(centroid_rows) == 0L) {
        centroid_info_b <- data.frame(
            block           = character(0),
            modality        = character(0),
            modality_type   = character(0),
            cluster         = integer(0),
            n_genes_cluster = integer(0),
            n_genes_used    = integer(0),
            coverage        = numeric(0),
            mean_R2         = numeric(0),
            sd_R2           = numeric(0),
            r2_member       = I(list()),
            centroid        = I(list()),
            stringsAsFactors = FALSE
        )
    } else {
        centroid_info_b <- do.call(
            rbind,
            centroid_rows
        )
    }
    
    centroid_info_b
}


#' Validate the `gene_mode` argument
#'
#' @noRd
#'
#' @description
#' Internal helper that checks the validity of the \code{gene_mode}
#' argument used by \code{cluster_genes_multiomics()}.  
#'
#' The mode determines how genes across layers within a block should be
#' combined and must be one of the two supported options.
#'
#' @param gene_mode
#' Character scalar specifying how genes are aggregated across layers.  
#' Must be either \code{"intersection"} or \code{"union"}.  
#' The value is validated using \code{\link[base]{match.arg}}.
#'
#' @return
#' Invisibly returns the validated \code{gene_mode} value.
#' 
.check_cluster_genes_multiomics_gene_mode <- function(gene_mode) {
    if (!is.character(gene_mode) || length(gene_mode) != 1L) {
        stop(
            "`gene_mode` must be a single character string ",
            "('intersection' or 'union').",
            call. = FALSE
        )
    }
    gene_mode <- match.arg(gene_mode, c("intersection", "union"))
    invisible(gene_mode)
}


#' Validate the `blocks` argument
#'
#' @noRd
#'
#' @description
#' Internal helper that validates the structure and content of the
#' \code{blocks} argument used by \code{cluster_genes_multiomics()}.  
#'
#' The function checks:
#' \itemize{
#'   \item that \code{blocks} is a non-empty named list of blocks,
#'   \item that each block is a non-empty named list of layers, and
#'   \item that each layer is a numeric matrix with non-empty row names
#'         and non-zero dimensions.
#' }
#'
#' @param blocks
#' Named list of blocks, one element per block.  
#' Each block contains a named list of layers, where each layer is a
#' numeric matrix (\code{features x spline_points}) with valid row names.
#'
#' @return
#' Invisibly returns \code{blocks} if all checks pass.
#' 
.check_cluster_genes_multiomics_blocks <- function(blocks) {
    # blocks: outer structure
    if (!is.list(blocks) || length(blocks) == 0L) {
        stop(
            "`blocks` must be a non-empty named list ",
            "(one element per block).",
            call. = FALSE
        )
    }
    if (is.null(names(blocks)) || any(names(blocks) == "")) {
        stop(
            "`blocks` must be a named list; all blocks must have ",
            "non-empty names.",
            call. = FALSE
        )
    }
    
    # blocks: inner structure (type/shape only)
    for (block_name in names(blocks)) {
        block_obj <- blocks[[block_name]]
        
        if (!is.list(block_obj) || length(block_obj) == 0L) {
            stop(
                "Each element of `blocks` must be a non-empty list of ",
                "layers. Block '", block_name, "' is not valid.",
                call. = FALSE
            )
        }
        if (is.null(names(block_obj)) || any(names(block_obj) == "")) {
            stop(
                "Each block in `blocks` must be a named list of layers. ",
                "Block '", block_name, "' has unnamed layers.",
                call. = FALSE
            )
        }
        
        for (layer_name in names(block_obj)) {
            mat <- block_obj[[layer_name]]
            
            if (!is.matrix(mat)) {
                stop(
                    "Layer '", layer_name, "' in block '", block_name,
                    "' must be a numeric matrix (features x spline_points).",
                    call. = FALSE
                )
            }
            if (!is.numeric(mat)) {
                stop(
                    "Layer '", layer_name, "' in block '", block_name,
                    "' must be a numeric matrix.",
                    call. = FALSE
                )
            }
            if (nrow(mat) == 0L || ncol(mat) == 0L) {
                stop(
                    "Layer '", layer_name, "' in block '", block_name,
                    "' has zero rows or columns.",
                    call. = FALSE
                )
            }
            if (is.null(rownames(mat)) || any(rownames(mat) == "")) {
                stop(
                    "Layer '", layer_name, "' in block '", block_name,
                    "' must have rownames (feature or gene IDs).",
                    call. = FALSE
                )
            }
        }
    }
    
    invisible(blocks)
}


#' Validate the `block_meta` argument
#'
#' @noRd
#'
#' @description
#' Internal helper that validates the block-level metadata supplied in
#' \code{block_meta}.  
#'
#' The function checks:
#' \itemize{
#'   \item that \code{block_meta} is a data frame with required columns
#'         \code{block}, \code{block_k}, \code{result_category},
#'         \code{cond1}, and \code{cond2},
#'   \item that \code{block_k} values are numeric, positive, non-missing,
#'         and constant within each block,
#'   \item that \code{result_category} is numeric, and
#'   \item that \code{cond1} and \code{cond2} are character vectors.
#' }
#'
#' @param block_meta
#' Data frame containing block-level metadata used by
#' \code{cluster_genes_multiomics()}.  
#' Must include the required columns and satisfy the structural and
#' consistency constraints listed above.
#'
#' @return
#' Invisibly returns \code{block_meta} if all checks pass.
#'
.check_cluster_genes_multiomics_block_meta <- function(block_meta) {
    # block_meta: structure
    if (!is.data.frame(block_meta)) {
        stop(
            "`block_meta` must be a data.frame (or tibble) with ",
            "block-level metadata.",
            call. = FALSE
        )
    }
    required_block_cols <- c(
        "block", "block_k", "result_category", "cond1", "cond2"
    )
    missing_block_cols <- setdiff(required_block_cols, colnames(block_meta))
    if (length(missing_block_cols) > 0L) {
        stop(
            "`block_meta` is missing required columns: ",
            paste(missing_block_cols, collapse = ", "),
            call. = FALSE
        )
    }
    
    # block_k checks
    if (!is.numeric(block_meta$block_k)) {
        stop(
            "`block_meta$block_k` must be numeric (positive integers).",
            call. = FALSE
        )
    }
    if (any(is.na(block_meta$block_k))) {
        stop(
            "`block_meta$block_k` contains NA. Each block must have a ",
            "defined number of clusters.",
            call. = FALSE
        )
    }
    if (any(block_meta$block_k <= 0)) {
        stop(
            "`block_meta$block_k` must contain positive values.",
            call. = FALSE
        )
    }
    
    # block_k constant within block
    by_block_k <- tapply(
        block_meta$block_k,
        block_meta$block,
        function(x) length(unique(x))
    )
    if (any(by_block_k > 1L)) {
        bad_blocks <- names(by_block_k)[by_block_k > 1L]
        stop(
            "`block_meta$block_k` must be identical for all rows of a ",
            "block. Blocks with inconsistent `block_k`: ",
            paste(bad_blocks, collapse = ", "),
            call. = FALSE
        )
    }
    
    # result_category numeric
    if (!is.numeric(block_meta$result_category)) {
        stop(
            "`block_meta$result_category` must be numeric (e.g. 1, 3).",
            call. = FALSE
        )
    }
    
    # cond1/cond2 character
    if (!is.character(block_meta$cond1)) {
        stop(
            "`block_meta$cond1` must be a character vector ",
            "(condition names).",
            call. = FALSE
        )
    }
    if (!is.character(block_meta$cond2)) {
        stop(
            "`block_meta$cond2` must be a character vector ",
            "(condition names or NA).",
            call. = FALSE
        )
    }
    
    invisible(block_meta)
}


#' Validate the `modality_meta` argument
#'
#' @noRd
#'
#' @description
#' Internal helper that validates the layer-level metadata supplied in
#' \code{modality_meta}.  
#'
#' The function checks:
#' \itemize{
#'   \item that \code{modality_meta} is a data frame with required columns
#'         \code{block}, \code{layer}, \code{layer_k}, and \code{layer_w},
#'   \item that \code{layer_k} is numeric, positive where non-NA, and may
#'         be \code{NA} to indicate one-to-one gene-level layers,
#'   \item that \code{layer_w} is numeric, non-negative, non-missing,
#'         and sums to a strictly positive value within each block, and
#'   \item that each \code{(block, layer)} combination appears exactly
#'         once.
#' }
#'
#' @param modality_meta
#' Data frame containing layer-level metadata used by
#' \code{cluster_genes_multiomics()}.  
#' Must include valid entries for \code{layer_k} and \code{layer_w},
#' and uniquely identify each layer within each block.
#'
#' @return
#' Invisibly returns \code{modality_meta} if all checks pass.
#'
.check_cluster_genes_multiomics_modality_meta <- function(modality_meta) {
    # modality_meta: structure
    if (!is.data.frame(modality_meta)) {
        stop(
            "`modality_meta` must be a data.frame (or tibble) with ",
            "layer-level metadata.",
            call. = FALSE
        )
    }
    required_layer_cols <- c("block", "layer", "layer_k", "layer_w")
    missing_layer_cols <- setdiff(required_layer_cols, colnames(modality_meta))
    if (length(missing_layer_cols) > 0L) {
        stop(
            "`modality_meta` is missing required columns: ",
            paste(missing_layer_cols, collapse = ", "),
            call. = FALSE
        )
    }
    
    # layer_k: NA or positive numeric
    if (!is.numeric(modality_meta$layer_k)) {
        stop(
            "`modality_meta$layer_k` must be numeric (positive integers) ",
            "or NA.",
            call. = FALSE
        )
    }
    if (any(modality_meta$layer_k[!is.na(modality_meta$layer_k)] <= 0)) {
        stop(
            "`modality_meta$layer_k` must be positive where it is not NA.",
            call. = FALSE
        )
    }
    
    # layer_w: numeric, non-negative; per-block sum > 0
    if (!is.numeric(modality_meta$layer_w)) {
        stop("`modality_meta$layer_w` must be numeric.", call. = FALSE)
    }
    if (any(is.na(modality_meta$layer_w))) {
        stop(
            "`modality_meta$layer_w` contains NA. Please provide weights ",
            "for all rows.",
            call. = FALSE
        )
    }
    if (any(modality_meta$layer_w < 0)) {
        stop(
            "`modality_meta$layer_w` must be non-negative.",
            call. = FALSE
        )
    }
    
    w_by_block <- tapply(modality_meta$layer_w, modality_meta$block, sum)
    if (any(w_by_block <= 0)) {
        bad_blocks <- names(w_by_block)[w_by_block <= 0]
        stop(
            "For each block, the sum of `modality_meta$layer_w` must be > 0. ",
            "Blocks with invalid weights: ",
            paste(bad_blocks, collapse = ", "),
            call. = FALSE
        )
    }
    
    # ensure (block, layer) pairs in modality_meta are unique
    bl_pairs <- paste(modality_meta$block, modality_meta$layer, sep = "||")
    if (any(duplicated(bl_pairs))) {
        dup <- unique(bl_pairs[duplicated(bl_pairs)])
        stop(
            "Duplicate (block, layer) combinations in `modality_meta`: ",
            paste(dup, collapse = ", "),
            ". Each (block, layer) must appear only once.",
            call. = FALSE
        )
    }
    
    invisible(modality_meta)
}


#' Cross-validate `blocks`, `block_meta`, and `modality_meta`
#'
#' @noRd
#'
#' @description
#' Internal helper that performs cross-argument consistency checks across
#' \code{blocks}, \code{block_meta}, and \code{modality_meta}.  
#'
#' The function verifies:
#' \itemize{
#'   \item that all blocks referenced in \code{block_meta} and
#'         \code{modality_meta} exist in \code{blocks},
#'   \item that every block in \code{block_meta} has at least one
#'         corresponding row in \code{modality_meta},
#'   \item that all blocks mentioned in \code{modality_meta} also appear in
#'         \code{block_meta},
#'   \item and issues a warning for blocks present in \code{blocks} but
#'         not referenced in \code{block_meta}.
#' }
#'
#' In addition, it performs per-layer checks using metadata:
#' \itemize{
#'   \item verifies that no layer matrices contain \code{NA} values,
#'   \item enforces that one-to-one gene-level layers (with
#'         \code{layer_k = NA}) have unique gene IDs (no duplicated
#'         rownames), and
#'   \item enforces that many-to-one layers (with positive
#'         \code{layer_k}) use rownames of the form \code{"<gene>_<site>"}
#'         and have at least \code{max(2, layer_k)} features.
#' }
#'
#' Finally, it ensures that within each block, all gene-level layers
#' (one-to-one layers) share at least one common gene ID.
#'
#' @param blocks
#' Named list of blocks as used by \code{cluster_genes_multiomics()},
#' where each block contains a named list of layer matrices.
#'
#' @param block_meta
#' Data frame with block-level metadata, including at least a
#' \code{block} column that must match the block identifiers in
#' \code{blocks}.
#'
#' @param modality_meta
#' Data frame with layer-level metadata, including \code{block},
#' \code{layer}, and \code{layer_k} columns that must be consistent with
#' \code{blocks} and \code{block_meta}.
#'
#' @return
#' Invisibly returns \code{TRUE} if all cross-argument checks pass.
#'
.check_cluster_genes_multiomics_cross_args <- function(
        blocks,
        block_meta,
        modality_meta
) {
    # block_meta: blocks exist in blocks
    unknown_blocks <- setdiff(block_meta$block, names(blocks))
    if (length(unknown_blocks) > 0L) {
        stop(
            "Blocks in `block_meta$block` not present in `blocks`: ",
            paste(unique(unknown_blocks), collapse = ", "),
            call. = FALSE
        )
    }
    
    # each (block, layer) in modality_meta exists in blocks
    for (i in seq_len(nrow(modality_meta))) {
        b <- modality_meta$block[i]
        l <- modality_meta$layer[i]
        
        if (!b %in% names(blocks)) {
            stop(
                "Row ", i, " of `modality_meta` refers to block '", b,
                "', which is not present in `blocks`.",
                call. = FALSE
            )
        }
        if (!l %in% names(blocks[[b]])) {
            stop(
                "Row ", i, " of `modality_meta` refers to layer '", l,
                "' in block '", b, "', which is not present in ",
                "`blocks[[\"",
                b, "\"]]`.",
                call. = FALSE
            )
        }
    }
    
    # each block in block_meta must appear in modality_meta
    blocks_without_layers <- setdiff(
        block_meta$block,
        unique(modality_meta$block)
    )
    if (length(blocks_without_layers) > 0L) {
        stop(
            "Blocks in `block_meta` with no corresponding rows in ",
            "`modality_meta`: ",
            paste(blocks_without_layers, collapse = ", "),
            call. = FALSE
        )
    }
    
    # consistency: blocks in modality_meta appear in block_meta
    unknown_in_block_meta <- setdiff(
        unique(modality_meta$block),
        block_meta$block
    )
    if (length(unknown_in_block_meta) > 0L) {
        stop(
            "Blocks in `modality_meta$block` not present in ",
            "`block_meta$block`: ",
            paste(unknown_in_block_meta, collapse = ", "),
            call. = FALSE
        )
    }
    
    # warn on blocks defined in blocks but not in block_meta
    unused_blocks <- setdiff(names(blocks), block_meta$block)
    if (length(unused_blocks) > 0L) {
        warning(
            "Blocks present in `blocks` but not referenced in ",
            "`block_meta`: ",
            paste(unused_blocks, collapse = ", "),
            call. = FALSE
        )
    }
    
    # per-layer checks using meta: NA, one-to-one vs many-to-one
    for (i in seq_len(nrow(modality_meta))) {
        b  <- modality_meta$block[i]
        l  <- modality_meta$layer[i]
        lk <- modality_meta$layer_k[i]
        mat <- blocks[[b]][[l]]
        
        if (any(is.na(mat))) {
            stop(
                "Layer '", l, "' in block '", b,
                "' contains NA values. Please remove or impute ",
                "missing values before calling `cluster_genes_multiomics()`.",
                call. = FALSE
            )
        }
        
        rn <- rownames(mat)
        
        if (is.na(lk)) {
            # one-to-one gene-level layer: unique gene IDs
            if (any(duplicated(rn))) {
                stop(
                    "Layer '", l, "' in block '", b,
                    "' is marked as one-to-one (layer_k is NA) but has ",
                    "duplicated rownames. Gene-level layers must have ",
                    "one row per gene.",
                    call. = FALSE
                )
            }
        } else {
            # many-to-one layer: "<gene>_<site>" and enough features
            if (!all(grepl("_", rn, fixed = TRUE))) {
                stop(
                    "Layer '", l, "' in block '", b,
                    "' is marked as many-to-one (layer_k = ", lk,
                    ") but some rownames do not contain an underscore. ",
                    "Many-to-one layers must use rownames of the form ",
                    "'<gene>_<site>'.",
                    call. = FALSE
                )
            }
            if (nrow(mat) < 2L) {
                stop(
                    "Layer '", l, "' in block '", b,
                    "' is marked as many-to-one (layer_k = ", lk,
                    ") but has fewer than 2 features (rows). At least 2 ",
                    "mapped features are required for clustering.",
                    call. = FALSE
                )
            }
            if (nrow(mat) < lk) {
                stop(
                    "Layer '", l, "' in block '", b,
                    "' is marked as many-to-one with layer_k = ", lk,
                    " but only ", nrow(mat), " features are available. ",
                    "`layer_k` cannot exceed the number of features.",
                    call. = FALSE
                )
            }
        }
    }
    
    # ensure gene-level layers within each block share some genes
    for (b in unique(modality_meta$block)) {
        rows_b <- modality_meta$block == b &
            is.na(modality_meta$layer_k)
        layers_gene <- modality_meta$layer[rows_b]
        
        if (length(layers_gene) < 2L) {
            next
        }
        
        mats_b <- blocks[[b]][layers_gene]
        genes_lists <- lapply(mats_b, rownames)
        inter_b <- Reduce(intersect, genes_lists)
        
        if (length(inter_b) == 0L) {
            stop(
                "No common genes across gene-level layers for block '",
                b, "'. Check that rownames of gene-level matrices in ",
                "`blocks[[\"",
                b, "\"]]` refer to the same gene IDs.",
                call. = FALSE
            )
        }
    }
    
    invisible(TRUE)
}


# Level 3 function definitions -------------------------------------------------


#' Compute pairwise distances for an omics layer
#'
#' @noRd
#'
#' @description
#' Internal helper that computes a pairwise distance matrix between all
#' rows (features or genes) of a layer matrix.  
#' The matrix is first z-scored across rows, and Euclidean distances are
#' then computed using \code{\link[stats]{dist}}.
#'
#' @param layer_mat
#' Numeric matrix (\code{n_features x n_points}) where each row is a
#' feature- or gene-level trajectory (e.g. spline-evaluated time course).
#' Row names identify the features/genes and are propagated to the output
#' distance matrix.
#'
#' @return
#' A square numeric matrix \code{D} of size
#' \code{n_features x n_features}, where  
#' \code{D[i, j]} is the Euclidean distance between the scaled rows
#' \code{layer_mat[i, ]} and \code{layer_mat[j, ]}.  
#' Row and column names correspond to the row names of
#' \code{layer_mat}.
#'
#' @importFrom stats dist
#' 
.compute_layer_distance <- function(layer_mat) {
    mat_scaled <- scale(layer_mat)
    d <- stats::dist(
        mat_scaled,
        method = "euclidean"
    )
    D <- as.matrix(d)
}


#' Combine layer-specific distance matrices into a unified distance matrix
#'
#' @noRd
#'
#' @description
#' Internal helper that merges multiple layer-specific gene–gene distance
#' matrices into a single unified distance matrix for a block.  
#' Each layer contributes proportionally to its user-defined weight, and
#' all weights are internally normalized to sum to 1.
#'
#' @param dist_list
#' A list of square numeric matrices.  
#' Each matrix must represent pairwise distances between the same set of
#' genes (identical row and column names, identical ordering).
#'
#' @param weights
#' Numeric vector of the same length as \code{dist_list}, giving the
#' relative weight of each layer when building the combined distance
#' matrix.  
#' Values may be any non-negative numbers; they are internally
#' normalized so that the sum of weights equals 1.
#'
#' @return
#' A numeric matrix of the same dimensions as the matrices in
#' \code{dist_list}.  
#' The \code{(i, j)} entry is the weighted sum of the corresponding
#' distances from each layer:
#' \deqn{
#'   D_{ij} = \sum_{l=1}^{L} w_l \cdot D^{(l)}_{ij}
#' }
#' where \(w_l\) are the normalized weights.
#'
#' Row and column names are inherited from the first element of
#' \code{dist_list}.
#' 
.combine_layer_distances <- function(
        dist_list,
        weights
) {
    w <- weights / sum(weights)
    
    D_unified <- matrix(
        0,
        nrow = nrow(dist_list[[1L]]),
        ncol = ncol(dist_list[[1L]])
    )
    rownames(D_unified) <- rownames(dist_list[[1L]])
    colnames(D_unified) <- colnames(dist_list[[1L]])
    
    for (i in seq_along(dist_list)) {
        D_unified <- D_unified + w[i] * dist_list[[i]]
    }
    
    D_unified
}