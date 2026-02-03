#' Gene-Centric Multi-Omics Clustering Across Blocks and Data Modalities
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
#' Gene identifiers must be consistent across all one-to-one modalities;
#' otherwise, genes cannot be matched across omics modalities during distance
#' computation and clustering.
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
#' @param block_clusters
#' A named list specifying the amount of clusters per block. The list names
#' must match the names of \code{blocks}. Each element value specifies the
#' number of gene clusters (\code{k}) to compute for the corresponding
#' block.
#'
#' @param modality_meta
#' A data frame containing \emph{modality-level metadata}. One row per (block ×
#' modality). Must include:
#' \describe{
#'   \item{\code{block}}{Block identifier linking to \code{block_meta}.}
#'   \item{\code{modality}}{Modality name within the block.}
#'   \item{\code{many_to_one_k}}{Number of pattern clusters to use for building
#'   pattern signatures for many-to-one modalities. \code{NA} for modalities
#'   that are already gene-level.}
#'   \item{\code{modality_w}}{
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
#' A named list with three tibbles:
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
#'   and cluster, summarizing modality-specific cluster centroids and
#'   within-cluster coherence. Columns include the block and modality
#'   identifiers, cluster label, gene coverage statistics, the QC method
#'   used (\code{qc_method}; \code{"Pearson R2"} for one-to-one modalities
#'   and \code{"BC(HD)"} for many-to-one modalities), mean and standard
#'   deviation of per-gene QC values (\code{mean_qc}, \code{sd_qc}),
#'   optional per-gene QC vectors (\code{qc_member}) as a list-column, and
#'   the centroid representation stored as a list-column (\code{centroid}).}
#'
#'   \item{\code{many_to_one_clustering_qc}}{A tibble (or \code{NULL} if no
#'   many-to-one modalities are present) providing clustering quality
#'   diagnostics for many-to-one feature clustering steps. One row per
#'   block, modality, and feature-cluster, including the many-to-one
#'   \code{k} value used, the number of features used, the QC method
#'   (\code{"Pearson R2"}), mean and standard deviation of per-feature QC
#'   values (\code{mean_qc}, \code{sd_qc}), optional per-feature QC vectors
#'   (\code{qc_member}) as a list-column, and the centroid representation
#'   stored as a list-column (\code{centroid}).}
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
#' block_clusters <- list(
#'     time_Ctrl = 2L,
#'     time_Treat = 2L
#' )
#'
#' modality_meta <- data.frame(
#'     block   = c("time_Ctrl", "time_Treat"),
#'     modality   = c("rna", "rna"),
#'     many_to_one_k = c(NA_real_, NA_real_),
#'     modality_w = c(1, 1),
#'     stringsAsFactors = FALSE
#' )
#'
#' cluster_table <- cluster_genes_multiomics(
#'     blocks     = blocks,
#'     block_clusters = block_clusters,
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
        block_clusters,
        modality_meta,
        gene_mode = "intersection",
        verbose = TRUE
) {
    start_time <- Sys.time()
    
    .check_cluster_genes_multiomics_input(
        blocks = blocks,
        block_clusters = block_clusters,
        modality_meta = modality_meta,
        gene_mode = gene_mode,
        verbose = verbose
    )
    
    block_ids <- names(blocks)
    
    block_assignments <- vector("list", length(block_ids))
    names(block_assignments) <- block_ids
    
    centroid_list <- vector("list", length(block_ids))
    names(centroid_list) <- block_ids
    many_to_one_clustering_qc_list <- vector("list", length(block_ids))
    names(many_to_one_clustering_qc_list) <- block_ids
    
    for (i in seq_along(block_ids)) {
        b <- block_ids[i]
        if (isTRUE(verbose)) {
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
        
        k <- block_clusters[[b]]

        res_b <- .cluster_genes_multiomics_block(
            block_id = b,
            blocks = blocks,
            block_k = k,
            modality_meta = modality_meta,
            gene_mode = gene_mode,
            verbose = verbose
        )
        
        block_assignments[[b]] <- res_b$cl_b
        centroid_list[[b]] <- res_b$centroid_info
        many_to_one_clustering_qc_list[[b]] <- res_b$many_to_one_clustering_qc
    }
    
    genes_all <- sort(unique(unlist(lapply(block_assignments, names))))
    cluster_table <- data.frame(gene = genes_all, stringsAsFactors = FALSE)
    
    for (i in seq_along(block_ids)) {
        b <- block_ids[i]
        cl_b <- block_assignments[[b]]
        if (is.null(cl_b)) next
        
        col_name <- b
        if (!col_name %in% colnames(cluster_table)) {
            cluster_table[[col_name]] <- NA_integer_
        }
        
        genes_b <- names(cl_b)
        idx_tbl <- match(genes_b, cluster_table$gene)
        cluster_table[[col_name]][idx_tbl] <- as.integer(cl_b)
    }
    
    centroid_info <- do.call(rbind, centroid_list)
    many_to_one_clustering_qc_list <- Filter(
        Negate(is.null),
        many_to_one_clustering_qc_list
        )
    
    many_to_one_clustering_qc <- if (
        length(many_to_one_clustering_qc_list) > 0L
        ) {
        do.call(rbind, many_to_one_clustering_qc_list)
    } else {
        NULL
    }
    
    if (isTRUE(verbose)) {
        end_time <- Sys.time()
        elapsed <- end_time - start_time
        formatted <- format(elapsed, digits = 2)
        message("[cluster_genes_multiomics] total runtime: ", formatted)
    }
    
    list(
        cluster_table = tibble::as_tibble(cluster_table),
        centroid_info = tibble::as_tibble(centroid_info),
        many_to_one_clustering_qc =
            if (is.null(many_to_one_clustering_qc)) {
                NULL
            } else {
                tibble::as_tibble(many_to_one_clustering_qc)
            }
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
#'   \item cross-argument consistency checks linking blocks, modalities,
#'         and metadata.
#' }
#'
#' This function is a thin orchestrator that delegates to specialized
#' internal checkers:
#' \code{.check_cluster_genes_multiomics_gene_mode()},
#' \code{.check_cluster_genes_multiomics_blocks()},
#' \code{.check_cluster_genes_multiomics_block_clusters()},
#' \code{.check_cluster_genes_multiomics_modality_meta()}, and
#' \code{.check_cluster_genes_multiomics_cross_args()}.
#'
#' @param blocks
#' Named list of blocks, one element per block.  
#' Each block is itself a named list of modalities, where each modality is a
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
#' Data frame with modality-level metadata.  
#' Must contain the columns \code{block}, \code{modality}, \code{many_to_one_k},
#' and \code{modality_w}.  
#' The column \code{many_to_one_k} encodes whether a modality is many-to-one
#' (positive integer) or one-to-one (NA).  
#' The column \code{modality_w} provides non-negative modality weights whose
#' sum must be positive within each block.
#'
#' @param gene_mode
#' Character scalar specifying how genes should be combined across
#' modalities within a block.  
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
        block_clusters,
        modality_meta,
        gene_mode,
        verbose
) {
    .check_cluster_genes_multiomics_gene_mode(gene_mode)
    .check_cluster_genes_multiomics_blocks(blocks)
    .check_cluster_genes_multiomics_block_clusters(
        block_clusters = block_clusters,
        blocks = blocks
        )
    .check_cluster_genes_multiomics_modality_meta(modality_meta)
    .check_cluster_genes_multiomics_cross_args(
        blocks = blocks,
        modality_meta = modality_meta
        )
    
    # verbose must be TRUE or FALSE (length 1 logical)
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop_call_false(
            "`verbose` must be either TRUE or FALSE (a single logical value)."
        )
    }
    
    invisible(TRUE)
}


#' Cluster genes within a single block and compute modality-wise centroids
#'
#' @noRd
#'
#' @description
#' Internal helper that performs all clustering steps for a single block
#' in the multi-omics workflow.  
#'
#' For the specified block, the function:
#' \itemize{
#'   \item extracts the relevant block-level and modality-level metadata,
#'   \item converts many-to-one modalities into gene-level signatures,
#'   \item selects genes according to the requested \code{gene_mode}
#'         (\code{"intersection"} or \code{"union"}),
#'   \item constructs a block-level distance matrix by combining
#'         modality-level distances using per-modality weights, and
#'   \item runs the clustering algorithm to obtain block-specific
#'         gene–cluster assignments.
#' }
#'
#' After clustering, it computes per-modality cluster centroids and quality
#' metrics through \code{.compute_block_centroids()}.
#'
#' @param block_id
#' Identifier of the block being processed. Must match a value in
#' \code{block_meta$block} and \code{names(blocks)}.
#'
#' @param blocks
#' Named list of blocks as supplied to
#' \code{cluster_genes_multiomics()}.  
#' Each block is a named list of raw modality matrices.
#'
#' @param block_meta
#' Block-level metadata table.  
#' Used to obtain \code{block_k} and to filter metadata for the current
#' block.
#'
#' @param modality_meta
#' Modality-level metadata table.  
#' Used to retrieve modality names, modality types (\code{many_to_one_k}), and
#' modality weights for the current block.
#'
#' @param gene_mode
#' Character scalar specifying how to combine genes across modalities:
#' either \code{"intersection"} or \code{"union"}.  
#' Determines whether only shared genes or all genes across modalities
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
#'     Data frame with per-modality, per-cluster centroid trajectories and
#'     associated summary metrics, as produced by
#'     \code{.compute_block_centroids()}.
#'   }
#' }
#'
.cluster_genes_multiomics_block <- function(
        block_id,
        blocks,
        block_k,
        modality_meta,
        gene_mode,
        verbose
) {
    meta_b_modality <- 
        modality_meta[modality_meta$block == block_id, , drop = FALSE]
    
    modality_names <- meta_b_modality$modality
    modality_weights <- meta_b_modality$modality_w
    
    modality_mats_raw <- lapply(
        modality_names,
        function(ln) blocks[[block_id]][[ln]]
    )
    names(modality_mats_raw) <- modality_names
    
    modality_mats <- vector("list", length(modality_names))
    names(modality_mats) <- modality_names
    feature_cluster_qc_by_modality <- list()
    
    if (isTRUE(verbose)) {
        message(
            "  [block '",
            block_id,
            "'] building gene-level modality matrices..."
        )
    }
    
    many_to_one <- setNames(
        rep(FALSE, length(modality_names)),
        modality_names
    )
    
    for (i in seq_along(modality_names)) {
        ln <- modality_names[i]
        if (isTRUE(verbose)) message("    Modality: ", ln)
        
        lk <- meta_b_modality$many_to_one_k[i]
        mat_raw <- modality_mats_raw[[ln]]
        
        is_m2o <- !is.na(lk)
        many_to_one[[ln]] <- is_m2o
        
        if (!is_m2o) {
            # one-to-one: normalize trajectories (row-wise z-score)
            modality_mats[[ln]] <- .normalize_modality_mat(
                mat_raw,
                is_many_to_one = FALSE
                )
        } else {
            rn <- rownames(mat_raw)
            gene_ids <- sub("_.*$", "", rn)
            feature_to_gene <- gene_ids
            
            sig_obj <- .build_site_signatures(
                modality_mat     = mat_raw,
                feature_to_gene  = feature_to_gene,
                many_to_one_k    = lk
            )
            
            sig <- sig_obj$signatures
            sig <- sig[order(rownames(sig)), , drop = FALSE]
            # many-to-one: Hellinger normalize signatures
            modality_mats[[ln]] <- .normalize_modality_mat(
                sig,
                is_many_to_one = TRUE
            )
            
            # collect QC if present/non-empty
            qc <- sig_obj$many_to_one_clustering_qc
            if (is.data.frame(qc) && nrow(qc) > 0L) {
                qc$block <- block_id
                qc$modality <- ln
                qc$many_to_one_k <- lk
                feature_cluster_qc_by_modality[[ln]] <- qc
            }
        }
    }
    
    many_to_one_clustering_qc <- if (
        length(feature_cluster_qc_by_modality
               ) > 0L) {
        do.call(rbind, feature_cluster_qc_by_modality)
    } else {
        NULL
    }
    
    if (isTRUE(verbose)) {
        message(
            "  [block '",
            block_id,
            "'] building joint feature matrix for multi-omics clustering..."
        )
    }

    if (gene_mode == "intersection") {
        X_block <- .build_block_feature_matrix_intersection(
            modality_mats     = modality_mats,
            modality_names    = modality_names,
            modality_weights  = modality_weights,
            many_to_one       = many_to_one
        )
    } else {
        X_block <- .build_block_feature_matrix_union(
            modality_mats     = modality_mats,
            modality_names    = modality_names,
            modality_weights  = modality_weights,
            many_to_one       = many_to_one
        ) 
    }
    
    if (isTRUE(verbose)) {
        message(
            "  [block '",
            block_id,
            "'] clustering genes (k = ",
            block_k,
            ")..."
        )
    }

    cl_b <- .cluster_feature_matrix(
        feature_mat = X_block,
        k           = block_k
    )

    centroid_info_b <- .compute_block_centroids(
        block_id = block_id,
        modality_mats = modality_mats,
        cl_b = cl_b,
        meta_b_modality = meta_b_modality
    )
    centroid_info_b <- .rank_clusters(centroid_info_b)
    
    list(
        cl_b = cl_b,
        centroid_info = centroid_info_b,
        many_to_one_clustering_qc = many_to_one_clustering_qc
    )
}


# Level 2 function definitions -------------------------------------------------


#' Build gene-level pattern signatures from many-to-one features
#'
#' @noRd
#'
#' @description
#' Internal helper that converts a many-to-one omics modality (e.g. sites,
#' probes) into **gene-level pattern signatures**.
#'
#' All feature-level trajectories are first **clustered jointly** into
#' \code{many_to_one_k} global temporal archetypes based directly on their
#' numeric trajectories (after z-scoring). These archetypes represent
#' recurring dynamic patterns shared across the entire modality.
#'
#' Each gene is then represented by a fixed-length signature vector that
#' quantifies how its own features distribute across these shared
#' archetypes—effectively a mixture over global temporal patterns.
#'
#' @param modality_mat
#' Numeric matrix (\code{features x spline_points}) where each row is a
#' feature-level trajectory (e.g. phosphosite, probe).
#' Row names must uniquely identify features; no specific format is
#' required.
#'
#' @param feature_to_gene
#' Character vector of length \code{nrow(modality_mat)} mapping each feature
#' (row of \code{modality_mat}) to its corresponding gene.
#' The mapping is **strictly positional**: the \(i\)-th element of
#' \code{feature_to_gene} corresponds to the \(i\)-th row of
#' \code{modality_mat}. Names are ignored.
#'
#' @param many_to_one_k
#' Integer scalar giving the number of global temporal archetype clusters to
#' compute across all features.
#' Must satisfy \code{2 <= many_to_one_k < nrow(modality_mat)}.
#'
#' @return
#' A list with two components:
#'
#' \describe{
#'   \item{\code{signatures}}{
#'     Numeric matrix of dimension \code{genes x many_to_one_k}.
#'     Each row represents a gene; each column corresponds to one global
#'     temporal archetype.
#'
#'     Entries give the **fraction of that gene's features** assigned to
#'     each archetype. Rows sum to 1 for genes with at least one mapped
#'     feature.
#'   }
#'
#'   \item{\code{feature_clusters}}{
#'     Integer vector of length \code{nrow(modality_mat)} giving the
#'     archetype assignment for each feature.
#'     Names correspond to feature identifiers (rownames of
#'     \code{modality_mat}).
#'   }
#' }
#' 
.build_site_signatures <- function(
        modality_mat,
        feature_to_gene,
        many_to_one_k
) {
    many_to_one_k <- as.integer(many_to_one_k)
    
    feat_map <- as.character(feature_to_gene)
    if (length(feat_map) != nrow(modality_mat)) {
        stop(
            "`feature_to_gene` must have length equal to nrow(modality_mat). ",
            "Mapping is positional: element i is the gene ID for row i.",
            call. = FALSE
        )
    }

    feature_clusters <- .cluster_feature_matrix(
        feature_mat = scale(modality_mat),
        k           = many_to_one_k
    )

    many_to_one_clustering_qc <- .compute_cluster_centroids_qc(
        X  = scale(modality_mat),
        cl = feature_clusters,
        qc_method = "Pearson R2",
        center_scale_rows = TRUE,
        require = "intersection"
    )
    
    # 2) Build gene x K signature matrix (fractions per global archetype)
    genes <- sort(unique(feat_map))
    sig_mat <- matrix(
        0,
        nrow = length(genes),
        ncol = many_to_one_k
    )
    rownames(sig_mat) <- genes
    colnames(sig_mat) <- paste0("cluster_", seq_len(many_to_one_k))
    
    for (g in genes) {
        idx <- which(feat_map == g)
        if (length(idx) == 0L) next
        
        tab <- table(feature_clusters[idx])
        sig_mat[g, as.integer(names(tab))] <- as.numeric(tab)
        
        rs <- sum(sig_mat[g, ])
        if (rs > 0) sig_mat[g, ] <- sig_mat[g, ] / rs
    }
    
    list(
        signatures       = sig_mat,
        feature_clusters = feature_clusters,
        many_to_one_clustering_qc = many_to_one_clustering_qc
    )
}


#' Build joint block-level feature matrix (intersection mode)
#'
#' @noRd
#'
#' @description
#' Internal helper that constructs a joint gene-level feature matrix for
#' multi-omics clustering in \code{gene_mode = "intersection"}.
#'
#' Only genes present in \emph{all} modalities of the block are retained.
#' For this common gene set, the function:
#' \itemize{
#'   \item extracts gene-level representations per modality,
#'   \item applies modality-specific preprocessing and transformations,
#'   \item standardizes features within each modality,
#'   \item rescales modality blocks to enforce user-defined weights while
#'         correcting for differing feature dimensionalities, and
#'   \item concatenates all modality blocks into a single joint feature
#'         matrix suitable for Euclidean clustering.
#' }
#'
#' Many-to-one modalities (e.g. site or probe signatures) are automatically
#' transformed using a Hellinger (square-root) transform prior to
#' standardization, yielding a geometry appropriate for Euclidean methods.
#'
#' The resulting matrix can be passed directly to
#' \code{.cluster_feature_matrix()} for joint gene-centric clustering across
#' all modalities, without constructing gene--gene distance matrices.
#'
#' @param modality_mats
#' Named list of gene-level numeric matrices, one per modality. Row names
#' are gene identifiers; columns are modality-specific features (e.g.
#' spline points or signature components).
#'
#' @param modality_names
#' Character vector specifying which modalities to include and in what
#' order. Must correspond to names in \code{modality_mats}.
#'
#' @param modality_weights
#' Numeric vector of non-negative modality weights, same length and order
#' as \code{modality_names}. Weights control the relative contribution of
#' each modality to the joint feature space.
#'
#' @param many_to_one
#' Logical or character vector indicating which modalities represent
#' many-to-one gene summaries (e.g. site signatures). If logical, it must
#' be named by modality. If character, it is interpreted as the set of
#' many-to-one modality names.
#'
#' @return
#' A numeric matrix of dimension \code{genes x features}, where each row
#' corresponds to a gene present in all modalities and columns represent
#' the weighted, standardized, and concatenated feature blocks across
#' modalities. Row names are gene identifiers.
#' 
#' @importFrom purrr map reduce
#' @importFrom dplyr case_when
#' @importFrom rlang abort
#' @importFrom stats setNames
#'
.build_block_feature_matrix_intersection <- function(
        modality_mats,
        modality_names,
        modality_weights,
        many_to_one = NULL
) {
    if (is.null(many_to_one)) {
        many_to_one <- rep(FALSE, length(modality_names))
        names(many_to_one) <- modality_names
    } else if (is.logical(many_to_one)) {
        if (is.null(names(many_to_one))) {
            names(many_to_one) <- modality_names
        }
        many_to_one <- many_to_one[modality_names]
    } else {
        many_to_one <- stats::setNames(
            modality_names %in% many_to_one,
            modality_names
        )
    }

    genes_block <- modality_mats[modality_names] |>
        purrr::map(rownames) |>
        purrr::reduce(intersect) |>
        sort()
    
    w <- modality_weights / sum(modality_weights)
    names(w) <- modality_names
    
    blocks <- purrr::map(
        modality_names,
        function(m) {
            X <- modality_mats[[m]][genes_block, , drop = FALSE] |>
                as.matrix()
            p <- ncol(X)
            X * sqrt(w[[m]] / p)
        }
    )
    
    X_block <- purrr::reduce(blocks, cbind)
    rownames(X_block) <- genes_block
    
    X_block
}


#' Build joint block-level feature matrix (union mode)
#'
#' @noRd
#'
#' @description
#' Internal helper that constructs a joint gene-level feature matrix for
#' multi-omics clustering in \code{gene_mode = "union"}.
#'
#' All genes present in \emph{any} modality of the block are retained. For
#' each modality, genes not observed in that modality are represented by
#' zero-valued feature blocks after standardization, such that the modality
#' contributes no signal for those genes.
#'
#' For the union gene set, the function:
#' \itemize{
#'   \item aligns modality-specific gene-level feature matrices to the
#'         union of genes,
#'   \item applies modality-specific preprocessing and transformations,
#'   \item standardizes features within each modality,
#'   \item rescales modality blocks to enforce user-defined weights while
#'         correcting for differing feature dimensionalities, and
#'   \item concatenates all modality blocks into a single joint feature
#'         matrix suitable for Euclidean clustering.
#' }
#'
#' Many-to-one modalities (e.g. site or probe signatures) are automatically
#' transformed using a Hellinger (square-root) transform prior to
#' standardization, yielding a geometry appropriate for Euclidean methods.
#'
#' The resulting matrix can be passed directly to
#' \code{.cluster_feature_matrix()} for joint gene-centric clustering across
#' all modalities, without constructing gene--gene distance matrices.
#'
#' @param modality_mats
#' Named list of gene-level numeric matrices, one per modality. Row names
#' are gene identifiers; columns are modality-specific features (e.g.
#' spline points or signature components).
#'
#' @param modality_names
#' Character vector specifying which modalities to include and in what
#' order. Must correspond to names in \code{modality_mats}.
#'
#' @param modality_weights
#' Numeric vector of non-negative modality weights, same length and order
#' as \code{modality_names}. Weights control the relative contribution of
#' each modality to the joint feature space.
#'
#' @param many_to_one
#' Logical or character vector indicating which modalities represent
#' many-to-one gene summaries (e.g. site signatures). If logical, it must
#' be named by modality. If character, it is interpreted as the set of
#' many-to-one modality names.
#'
#' @return
#' A numeric matrix of dimension \code{genes x features}, where each row
#' corresponds to a gene present in at least one modality and columns
#' represent the weighted, standardized, and concatenated feature blocks
#' across modalities. Row names are gene identifiers.
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr case_when
#' @importFrom rlang abort
#' @importFrom stats setNames
#'
.build_block_feature_matrix_union <- function(
        modality_mats,
        modality_names,
        modality_weights,
        many_to_one = NULL
) {
    if (is.null(many_to_one)) {
        many_to_one <- rep(FALSE, length(modality_names))
        names(many_to_one) <- modality_names
    } else if (is.logical(many_to_one)) {
        if (is.null(names(many_to_one))) {
            names(many_to_one) <- modality_names
        }
        many_to_one <- many_to_one[modality_names]
    } else {
        many_to_one <- stats::setNames(
            modality_names %in% many_to_one,
            modality_names
        )
    }
    
    genes_block <- modality_mats[modality_names] |>
        purrr::map(rownames) |>
        purrr::reduce(union) |>
        sort()
    
    w <- modality_weights / sum(modality_weights)
    names(w) <- modality_names
    
    blocks <- purrr::map(
        modality_names,
        function(m) {
            mat <- modality_mats[[m]]
            X <- matrix(
                0,
                nrow = length(genes_block),
                ncol = ncol(mat),
                dimnames = list(genes_block, colnames(mat))
            )
            idx <- match(rownames(mat), genes_block)
            X[idx, ] <- as.matrix(mat)
            p <- ncol(X)
            X * sqrt(w[[m]] / p)
        }
    )
    
    X_block <- purrr::reduce(blocks, cbind)
    rownames(X_block) <- genes_block
    
    X_block
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
#' @return
#' An integer vector of cluster assignments, one entry per gene.
#' Names correspond to the row names of \code{dist_mat}.
#'
#' @importFrom stats kmeans cmdscale dist
#' @importFrom ClusterR MiniBatchKmeans predict_KMeans
#' 
.cluster_feature_matrix <- function(
        feature_mat,
        k,
        verbose = FALSE
) {
    n_obs <- nrow(feature_mat)
    if (n_obs < 2L) {
        stop("Need at least 2 observations for clustering.", call. = FALSE)
    }

    k_vec <- sort(unique(as.integer(k)))
    k_vec <- k_vec[k_vec > 1L & k_vec < n_obs]
    if (length(k_vec) == 0L) {
        stop(
            "`k` must contain values between 2 and n_obs - 1.",
            call. = FALSE
        )
    }
    
    x <- feature_mat
    if (anyNA(x) || any(!is.finite(x))) {
        stop(
            "`feature_mat` must not contain NA/NaN/Inf values.",
            call. = FALSE
        )
    }
    
    rn <- rownames(x)
    if (is.null(rn)) {
        rn <- paste0("row_", seq_len(nrow(x)))
        rownames(x) <- rn
    }
    
    use_mb <- n_obs > 1000L
    
    run_one_k <- function(kk) {
        if (isTRUE(verbose)) {
            message("Clustering with k = ", kk, " (n = ", n_obs, ")")
        }
        
        if (!use_mb) {
            fit <- stats::kmeans(
                x,
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
                data            = x,
                clusters        = kk,
                batch_size      = batch_size,
                num_init        = 10L,
                max_iters       = 300L,
                init_fraction   = 1.0,
                early_stop_iter = 10L,
                tol             = 1e-4,
                verbose         = FALSE,
                seed            = 42
            )
            list(
                tot_within = sum(fit$WCSS_per_cluster),
                cluster    = ClusterR::predict_KMeans(
                    data      = x,
                    CENTROIDS = fit$centroids
                )
            )
        }
    }

    if (length(k_vec) == 1L) {
        cl <- run_one_k(k_vec)$cluster
    } else {
        fits <- lapply(k_vec, run_one_k)
        
        tot_within <- vapply(
            fits,
            function(z) z$tot_within,
            numeric(1)
        )
        
        p <- ncol(x)
        bic <- n_obs * log(tot_within / n_obs) +
            k_vec * log(n_obs) * p
        
        best_idx <- which.min(bic)
        cl <- fits[[best_idx]]$cluster
    }
    
    names(cl) <- rownames(x)
    cl
}


#' Compute per-modality cluster centroids and QC metrics within a block
#'
#' @noRd
#'
#' @description
#' Internal wrapper that summarizes gene-level clusters for a single block
#' by computing modality-specific centroid representations and associated
#' quality-control (QC) metrics, then adding block- and modality-level
#' bookkeeping such as coverage and modality type.
#'
#' Core centroid and QC calculations are delegated to
#' \code{.compute_cluster_centroids_qc()} and are performed independently
#' for each modality using a modality-appropriate similarity measure:
#' \itemize{
#'   \item \strong{one-to-one modalities}: Pearson correlation
#'         coefficient squared (\eqn{R^2}),
#'   \item \strong{many-to-one modalities}: Bhattacharyya coefficient in
#'         Hellinger space (\code{"BC(HD)"}).
#' }
#'
#' It is assumed that modality matrices have already been normalized
#' upstream (row-wise z-scoring for one-to-one modalities; Hellinger
#' transformation for many-to-one modalities). No additional scaling is
#' applied here.
#'
#' For each modality and each cluster, the function:
#' \itemize{
#'   \item computes centroid vectors and per-gene QC values using
#'         \code{.compute_cluster_centroids_qc()},
#'   \item counts cluster membership (\code{n_genes_cluster}) from
#'         \code{cl_b},
#'   \item counts how many cluster genes are present in the modality
#'         matrix (\code{n_genes_used}),
#'   \item computes coverage as
#'         \code{n_genes_used / n_genes_cluster}.
#' }
#'
#' Modality types are inferred from \code{meta_b_modality$many_to_one_k}:
#' \code{NA} indicates a one-to-one gene-level modality; non-\code{NA}
#' values indicate many-to-one modalities.
#'
#' @param block_id
#' Character (or factor) scalar giving the identifier of the block for which
#' centroids and QC metrics are computed.
#'
#' @param modality_mats
#' Named list of numeric matrices, one per modality in the block. Each
#' matrix has dimensions \code{genes x features}, with row names
#' corresponding to gene identifiers. Matrices are assumed to be already
#' normalized according to modality type.
#'
#' @param cl_b
#' Named integer vector of cluster assignments for genes in the block.
#' Names must be gene identifiers used to match rows in
#' \code{modality_mats}; values are cluster labels.
#'
#' @param meta_b_modality
#' Data frame with modality-level metadata for this block. Must contain at
#' least the columns \code{modality} (matching
#' \code{names(modality_mats)}) and \code{many_to_one_k} for modality-type
#' inference.
#'
#' @return
#' A data frame with one row per \code{(modality, cluster)} combination in
#' the block and the following columns:
#' \describe{
#'   \item{\code{block}}{Block identifier (same as \code{block_id}).}
#'   \item{\code{modality}}{Modality name.}
#'   \item{\code{modality_type}}{Character label, either
#'         \code{"one_to_one"} or \code{"many_to_one"}.}
#'   \item{\code{cluster}}{Cluster label.}
#'   \item{\code{n_genes_cluster}}{Number of genes in the cluster according
#'         to \code{cl_b}.}
#'   \item{\code{n_genes_used}}{Number of cluster genes present in the
#'         modality matrix and used to compute centroids and QC metrics.}
#'   \item{\code{coverage}}{Fraction of cluster genes represented in the
#'         modality (\code{n_genes_used / n_genes_cluster}).}
#'   \item{\code{qc_method}}{QC metric used for this modality
#'         (\code{"Pearson R2"} or \code{"BC(HD)"}).}
#'   \item{\code{mean_qc}}{Mean QC value across genes used in this
#'         modality and cluster.}
#'   \item{\code{sd_qc}}{Standard deviation of the QC values across genes
#'         used in this modality and cluster.}
#'   \item{\code{qc_member}}{List-column; each entry is a named numeric
#'         vector of per-gene QC values (names are gene IDs) for genes used
#'         in this modality and cluster.}
#'   \item{\code{centroid}}{List-column; each entry is a numeric vector
#'         giving the centroid representation for this modality and
#'         cluster.}
#' }
#' 
.compute_block_centroids <- function(
        block_id,
        modality_mats,
        cl_b,
        meta_b_modality
) {
    out <- list()
    row_i <- 0L
    
    clusters <- sort(unique(as.integer(cl_b)))
    
    for (ln in names(modality_mats)) {
        mat_l <- modality_mats[[ln]]
        
        lk <- meta_b_modality$many_to_one_k[meta_b_modality$modality == ln][1L]
        modality_type <- if (is.na(lk)) "one_to_one" else "many_to_one"
        
        qc_method <- if (modality_type == "one_to_one") {
            "Pearson R2"
        } else {
            "BC(HD)"
        }
        
        stats_l <- .compute_cluster_centroids_qc(
            X  = mat_l,
            cl = cl_b,
            qc_method = qc_method,
            center_scale_rows = FALSE,
            require = "intersection"
        )
        
        for (c in clusters) {
            genes_cluster <- names(cl_b)[cl_b == c]
            n_cluster <- length(genes_cluster)
            
            genes_modality <- intersect(genes_cluster, rownames(mat_l))
            n_used <- length(genes_modality)
            
            coverage <- if (n_cluster > 0L) n_used / n_cluster else NA_real_
            
            st <- stats_l[stats_l$cluster == c, , drop = FALSE]
            if (nrow(st) == 0L) {
                st <- data.frame(
                    cluster   = c,
                    n_used    = 0L,
                    qc_method = qc_method,
                    mean_qc   = NA_real_,
                    sd_qc     = NA_real_,
                    qc_member = I(list(setNames(numeric(0), character(0)))),
                    centroid  = I(list(rep(NA_real_, ncol(mat_l)))),
                    stringsAsFactors = FALSE
                )
            }
            
            row_i <- row_i + 1L
            out[[row_i]] <- data.frame(
                block           = block_id,
                modality        = ln,
                modality_type   = modality_type,
                cluster         = c,
                n_genes_cluster = n_cluster,
                n_genes_used    = n_used,
                coverage        = coverage,
                qc_method       = st$qc_method,
                mean_qc         = st$mean_qc,
                sd_qc           = st$sd_qc,
                qc_member       = st$qc_member,
                centroid        = st$centroid,
                stringsAsFactors = FALSE
            )
        }
    }
    
    if (length(out) == 0L) {
        return(data.frame(
            block           = character(0),
            modality        = character(0),
            modality_type   = character(0),
            cluster         = integer(0),
            n_genes_cluster = integer(0),
            n_genes_used    = integer(0),
            coverage        = numeric(0),
            qc_method       = character(0),
            mean_qc         = numeric(0),
            sd_qc           = numeric(0),
            qc_member       = I(list()),
            centroid        = I(list()),
            stringsAsFactors = FALSE
        ))
    }
    
    do.call(rbind, out)
}


#' Validate the `gene_mode` argument
#'
#' @noRd
#'
#' @description
#' Internal helper that checks the validity of the \code{gene_mode}
#' argument used by \code{cluster_genes_multiomics()}.  
#'
#' The mode determines how genes across modalities within a block should be
#' combined and must be one of the two supported options.
#'
#' @param gene_mode
#' Character scalar specifying how genes are aggregated across modalities.  
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
#'   \item that each block is a non-empty named list of modalities, and
#'   \item that each modality is a numeric matrix with non-empty row names
#'         and non-zero dimensions.
#' }
#'
#' @param blocks
#' Named list of blocks, one element per block.  
#' Each block contains a named list of modalities, where each modality is a
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
                "modalities. Block '", block_name, "' is not valid.",
                call. = FALSE
            )
        }
        if (is.null(names(block_obj)) || any(names(block_obj) == "")) {
            stop(
                "Each block in `blocks` must be a named list of modalities. ",
                "Block '", block_name, "' has unnamed modalities.",
                call. = FALSE
            )
        }
        
        for (modality_name in names(block_obj)) {
            mat <- block_obj[[modality_name]]
            
            if (!is.matrix(mat)) {
                stop(
                    "Modality '", modality_name, "' in block '", block_name,
                    "' must be a numeric matrix (features x spline_points).",
                    call. = FALSE
                )
            }
            if (!is.numeric(mat)) {
                stop(
                    "Modality '", modality_name, "' in block '", block_name,
                    "' must be a numeric matrix.",
                    call. = FALSE
                )
            }
            if (nrow(mat) == 0L || ncol(mat) == 0L) {
                stop(
                    "Modality '", modality_name, "' in block '", block_name,
                    "' has zero rows or columns.",
                    call. = FALSE
                )
            }
            if (is.null(rownames(mat)) || any(rownames(mat) == "")) {
                stop(
                    "Modality '", modality_name, "' in block '", block_name,
                    "' must have rownames (feature or gene IDs).",
                    call. = FALSE
                )
            }
        }
    }
    
    invisible(blocks)
}


#' Validate the `block_clusters` argument
#'
#' @noRd
#'
#' @description
#' Internal helper that validates the block-level clustering parameters
#' supplied in \code{block_clusters}. The object must be a named list whose
#' names match the block names in \code{blocks}. Each list element specifies
#' the number of clusters \code{k} to compute for that block.
#'
#' The function checks:
#' \itemize{
#'   \item that \code{block_clusters} is a named list,
#'   \item that \code{blocks} is a named list (outer list = blocks),
#'   \item that block names match in both directions (no missing/excess),
#'   \item that each \code{k} is a scalar, integer-like, and > 0.
#' }
#'
#' @param block_clusters A named list. Names are block identifiers and must
#'   match \code{names(blocks)}. Values are the number of clusters \code{k}
#'   (integer-like scalars > 0) for the corresponding block.
#'
#' @param blocks A named nested list specifying all data used for clustering.
#'   The outer list corresponds to analytical \emph{blocks}. Each block is a
#'   named list of modalities, each a numeric matrix of dimension
#'   \code{features x spline_points}. Row naming conventions define how
#'   features map to genes and how genes are aligned across modalities.
#'
#' @return Invisibly returns \code{block_clusters} if all checks pass.
#' 
.check_cluster_genes_multiomics_block_clusters <- function(
        block_clusters,
        blocks
) {
    if (!is.list(block_clusters)) {
        stop("`block_clusters` must be a named list.")
    }
    if (is.null(names(block_clusters)) || any(names(block_clusters) == "")) {
        stop("`block_clusters` must have non-empty names (block identifiers).")
    }
    
    if (!is.list(blocks)) {
        stop("`blocks` must be a named list (outer list = blocks).")
    }
    if (is.null(names(blocks)) || any(names(blocks) == "")) {
        stop("`blocks` must have non-empty names (block identifiers).")
    }
    
    bc_names <- names(block_clusters)
    b_names <- names(blocks)
    
    missing_in_bc <- setdiff(b_names, bc_names)
    extra_in_bc <- setdiff(bc_names, b_names)
    
    if (length(missing_in_bc) > 0L) {
        stop(
            "Missing entries in `block_clusters` for blocks: ",
            paste(missing_in_bc, collapse = ", ")
        )
    }
    if (length(extra_in_bc) > 0L) {
        stop(
            "Unknown block names in `block_clusters` (not in `blocks`): ",
            paste(extra_in_bc, collapse = ", ")
        )
    }
    
    is_int_like1 <- function(x) {
        if (length(x) != 1L || is.na(x)) return(FALSE)
        if (is.integer(x)) return(TRUE)
        if (!is.numeric(x)) return(FALSE)
        identical(x, as.numeric(as.integer(x)))
    }
    
    for (bn in b_names) {
        k <- block_clusters[[bn]]
        
        if (!is_int_like1(k)) {
            stop(
                "Invalid `block_clusters[[", bn, "]]`: must be an int-like ",
                "scalar (e.g., 3L) and not NA."
            )
        }
        if (as.integer(k) <= 0L) {
            stop(
                "Invalid `block_clusters[[", bn, "]]`: must be > 0."
            )
        }
    }
    
    invisible(block_clusters)
}


#' Validate the `modality_meta` argument
#'
#' @noRd
#'
#' @description
#' Internal helper that validates the modality-level metadata supplied in
#' \code{modality_meta}.  
#'
#' The function checks:
#' \itemize{
#'   \item that \code{modality_meta} is a data frame with required columns
#'         \code{block}, \code{modality}, \code{many_to_one_k}, 
#'         and \code{modality_w},
#'   \item that \code{many_to_one_k} is numeric, positive where non-NA, and may
#'         be \code{NA} to indicate one-to-one gene-level modalities,
#'   \item that \code{modality_w} is numeric, non-negative, non-missing,
#'         and sums to a strictly positive value within each block, and
#'   \item that each \code{(block, modality)} combination appears exactly
#'         once.
#' }
#'
#' @param modality_meta
#' Data frame containing modality-level metadata used by
#' \code{cluster_genes_multiomics()}.  
#' Must include valid entries for \code{many_to_one_k} and \code{modality_w},
#' and uniquely identify each modality within each block.
#'
#' @return
#' Invisibly returns \code{modality_meta} if all checks pass.
#'
.check_cluster_genes_multiomics_modality_meta <- function(modality_meta) {
    # modality_meta: structure
    if (!is.data.frame(modality_meta)) {
        stop(
            "`modality_meta` must be a data.frame (or tibble) with ",
            "modality-level metadata.",
            call. = FALSE
        )
    }
    required_modality_cols <- c(
        "block", 
        "modality",
        "many_to_one_k",
        "modality_w"
        )
    missing_modality_cols <- setdiff(
        required_modality_cols,
        colnames(modality_meta)
        )
    if (length(missing_modality_cols) > 0L) {
        stop(
            "`modality_meta` is missing required columns: ",
            paste(missing_modality_cols, collapse = ", "),
            call. = FALSE
        )
    }
    
    # many_to_one_k: NA or positive numeric
    if (!is.numeric(modality_meta$many_to_one_k)) {
        stop(
            "`modality_meta$many_to_one_k` must be numeric ",
            "(positive integers) ",
            "or NA.",
            call. = FALSE
        )
    }
    if (any(modality_meta$many_to_one_k[!is.na(modality_meta$many_to_one_k)]
            <= 0)) {
        stop(
            "`modality_meta$many_to_one_k` must be positive when not NA.",
            call. = FALSE
        )
    }
    
    # modality_w: numeric, non-negative; per-block sum > 0
    if (!is.numeric(modality_meta$modality_w)) {
        stop("`modality_meta$modality_w` must be numeric.", call. = FALSE)
    }
    if (any(is.na(modality_meta$modality_w))) {
        stop(
            "`modality_meta$modality_w` contains NA. Please provide weights ",
            "for all rows.",
            call. = FALSE
        )
    }
    if (any(modality_meta$modality_w < 0)) {
        stop(
            "`modality_meta$modality_w` must be non-negative.",
            call. = FALSE
        )
    }
    
    w_by_block <- tapply(modality_meta$modality_w, modality_meta$block, sum)
    if (any(w_by_block <= 0)) {
        bad_blocks <- names(w_by_block)[w_by_block <= 0]
        stop(
            "For each block, the sum of `modality_meta$modality_w` ",
            "must be > 0. ",
            "Blocks with invalid weights: ",
            paste(bad_blocks, collapse = ", "),
            call. = FALSE
        )
    }
    
    # ensure (block, modality) pairs in modality_meta are unique
    bl_pairs <- paste(modality_meta$block, modality_meta$modality, sep = "||")
    if (any(duplicated(bl_pairs))) {
        dup <- unique(bl_pairs[duplicated(bl_pairs)])
        stop(
            "Duplicate (block, modality) combinations in `modality_meta`: ",
            paste(dup, collapse = ", "),
            ". Each (block, modality) must appear only once.",
            call. = FALSE
        )
    }
    
    invisible(modality_meta)
}


#' Cross-validate `blocks` and `modality_meta`
#'
#' @noRd
#'
#' @description
#' Internal helper that performs cross-argument consistency checks across
#' \code{blocks} and \code{modality_meta}.
#'
#' The function verifies:
#' \itemize{
#'   \item that \code{blocks} is a named list of blocks and each block is a
#'   named list of numeric matrices,
#'   \item that every (block, modality) in \code{modality_meta} exists in
#'   \code{blocks},
#'   \item that each block has metadata rows for all its modalities (and no
#'   extra rows for unknown modalities),
#'   \item and performs per-modality checks driven by \code{modality_meta}.
#' }
#'
#' Per-modality checks:
#' \itemize{
#'   \item no \code{NA} values in modality matrices,
#'   \item one-to-one modalities (\code{many_to_one_k = NA}) have
#'    unique rownames,
#'   \item many-to-one modalities (\code{many_to_one_k > 0}) have rownames
#'   containing
#'   an underscore and have at least \code{max(2, many_to_one_k)} rows.
#' }
#'
#' Finally, within each block, all one-to-one modalities share at least one
#' common gene ID.
#'
#' @param blocks Named list of blocks as used by
#'   \code{cluster_genes_multiomics()}, where each block contains a named
#'   list of modality matrices.
#'
#' @param modality_meta Data frame with modality-level metadata, including at
#'   least \code{block}, \code{modality}, and \code{many_to_one_k}. Entries
#'    must be consistent with \code{blocks}.
#'
#' @return Invisibly returns \code{TRUE} if all checks pass.
#' 
.check_cluster_genes_multiomics_cross_args <- function(
        blocks,
        modality_meta
) {
    if (!is.list(blocks) ||
        is.null(names(blocks)) ||
        any(names(blocks) == "")) {
        stop("`blocks` must be a named list (outer list = blocks).")
    }
    
    if (!is.data.frame(modality_meta)) {
        stop("`modality_meta` must be a data.frame/tibble.")
    }
    
    req <- c("block", "modality", "many_to_one_k")
    miss <- setdiff(req, names(modality_meta))
    if (length(miss) > 0L) {
        stop(
            "`modality_meta` is missing required columns: ",
            paste(miss, collapse = ", ")
        )
    }
    
    b_names <- names(blocks)
    
    for (b in b_names) {
        blk <- blocks[[b]]
        if (!is.list(blk) ||
            is.null(names(blk)) ||
            any(names(blk) == "")) {
            stop(
                "`blocks[[", b,
                "]]` must be a named list of modality matrices."
            )
        }
    }
    
    for (i in seq_len(nrow(modality_meta))) {
        b <- modality_meta$block[i]
        l <- modality_meta$modality[i]
        
        if (!is.character(b) ||
            length(b) != 1L ||
            is.na(b) ||
            b == "") {
            stop("`modality_meta$block` must contain non-empty strings.")
        }
        
        if (!is.character(l) ||
            length(l) != 1L ||
            is.na(l) ||
            l == "") {
            stop("`modality_meta$modality` must contain non-empty strings.")
        }
        
        if (!b %in% b_names) {
            stop(
                "Row ", i,
                " of `modality_meta` refers to unknown block '",
                b, "'."
            )
        }
        
        if (!l %in% names(blocks[[b]])) {
            stop(
                "Row ", i,
                " of `modality_meta` refers to unknown modality '",
                l, "' in block '", b, "'."
            )
        }
    }
    
    for (b in b_names) {
        rows_b <- modality_meta$block == b
        if (!any(rows_b)) {
            stop(
                "Block '", b,
                "' has no corresponding rows in `modality_meta`."
            )
        }
        
        modalitys_meta <- modality_meta$modality[rows_b]
        modalitys_blk <- names(blocks[[b]])
        
        miss_modalitys <- setdiff(modalitys_blk, modalitys_meta)
        extra_modalitys <- setdiff(modalitys_meta, modalitys_blk)
        
        if (length(miss_modalitys) > 0L) {
            stop(
                "Block '", b,
                "' missing modality rows in `modality_meta`: ",
                paste(miss_modalitys, collapse = ", ")
            )
        }
        
        if (length(extra_modalitys) > 0L) {
            stop(
                "Block '", b,
                "' has unknown modalities in `modality_meta`: ",
                paste(extra_modalitys, collapse = ", ")
            )
        }
    }
    
    for (i in seq_len(nrow(modality_meta))) {
        b <- modality_meta$block[i]
        l <- modality_meta$modality[i]
        lk <- modality_meta$many_to_one_k[i]
        mat <- blocks[[b]][[l]]
        
        if (!is.matrix(mat) || !is.numeric(mat)) {
            stop(
                "Modality '", l,
                "' in block '", b,
                "' must be a numeric matrix."
            )
        }
        
        if (any(is.na(mat))) {
            stop(
                "Modality '", l,
                "' in block '", b,
                "' contains NA values."
            )
        }
        
        rn <- rownames(mat)
        if (is.null(rn) ||
            any(is.na(rn)) ||
            any(rn == "")) {
            stop(
                "Modality '", l,
                "' in block '", b,
                "' must have non-empty rownames."
            )
        }
        
        if (is.na(lk)) {
            if (any(duplicated(rn))) {
                stop(
                    "Modality '", l,
                    "' in block '", b,
                    "' is gene-level (many_to_one_k is NA) but has ",
                    "duplicated ",
                    "rownames."
                )
            }
        } else {
            if (!is.numeric(lk) ||
                length(lk) != 1L ||
                is.na(lk) ||
                lk <= 0) {
                stop(
                    "`many_to_one_k` must be NA or a positive numeric scalar ",
                    "for ",
                    "modality '", l,
                    "' in block '", b, "'."
                )
            }
            
            if (!all(grepl("_", rn, fixed = TRUE))) {
                stop(
                    "Many-to-one modality '", l,
                    "' in block '", b,
                    "' must use rownames '<gene>_<feature>'."
                )
            }
            
            min_n <- max(2L, as.integer(lk))
            if (nrow(mat) < min_n) {
                stop(
                    "Many-to-one modality '", l,
                    "' in block '", b,
                    "' has too few rows. Need at least ",
                    min_n, "."
                )
            }
            
            if (nrow(mat) < lk) {
                stop(
                    "Many-to-one modality '", l,
                    "' in block '", b,
                    "' has nrow < many_to_one_k. Reduce `many_to_one_k` ",
                    "or add ",
                    "features."
                )
            }
        }
    }
    
    for (b in b_names) {
        rows_b <- modality_meta$block == b &
            is.na(modality_meta$many_to_one_k)
        modalitys_gene <- modality_meta$modality[rows_b]
        
        if (length(modalitys_gene) < 2L) next
        
        mats_b <- blocks[[b]][modalitys_gene]
        genes_lists <- lapply(mats_b, rownames)
        inter_b <- Reduce(intersect, genes_lists)
        
        if (length(inter_b) == 0L) {
            stop(
                "No common genes across gene-level modalities for block '",
                b, "'."
            )
        }
    }
    
    invisible(TRUE)
}


#' Rank cluster centroid summaries by overall within-cluster coherence
#'
#' @noRd
#'
#' @description
#' Internal helper that ranks cluster-level centroid summaries by their
#' overall coherence across modalities.
#'
#' For each cluster, the function computes the average of the per-modality
#' \code{mean_R2} values (ignoring missing values), yielding a single
#' summary statistic that reflects how consistently cluster members align
#' with their centroids across all available omics layers. The input table
#' is then augmented with this average metric and sorted such that clusters
#' with the highest overall coherence appear first.
#'
#' @param centroid_info_b
#' Data frame as returned by \code{.compute_block_centroids()}, containing
#' one row per \code{(block, modality, cluster)} combination and at least
#' the columns \code{cluster} and \code{mean_R2}.
#'
#' @return
#' A data frame with the same rows as \code{centroid_info_b}, augmented by an
#' additional column \code{avg_mean_R2} giving the mean of \code{mean_R2}
#' across modalities for each cluster. Rows are ordered in decreasing order
#' of \code{avg_mean_R2}, with ties broken by cluster label and modality.
#' 
#' @importFrom stats aggregate
#'
.rank_clusters <- function(centroid_info_b) {
    score <- centroid_info_b$mean_qc
    is_hel <- centroid_info_b$qc_method == "BC(HD)"
    score[is_hel] <- 1 - score[is_hel]
    
    tmp <- centroid_info_b
    tmp$qc_score <- score
    
    avg <- stats::aggregate(
        qc_score ~ cluster,
        data = tmp,
        FUN = function(x) mean(x, na.rm = TRUE)
    )
    names(avg)[names(avg) == "qc_score"] <- "avg_qc_score"
    
    out <- merge(
        tmp,
        avg,
        by = "cluster",
        all.x = TRUE,
        sort = FALSE
    )
    
    out <- out[order(-out$avg_qc_score, out$cluster, out$modality), ,
               drop = FALSE]
    
    out$qc_score <- NULL
    out
}


# Level 3 function definitions -------------------------------------------------


#' Compute cluster centroids and per-member QC values
#'
#' @noRd
#'
#' @description
#' Internal helper that computes centroid representations and per-feature
#' quality-control (QC) metrics for clusters defined on the rows of a
#' numeric matrix. The function operates on a single matrix and a named
#' cluster assignment vector and is agnostic to any higher-level biological,
#' block, or modality-specific context.
#'
#' For each cluster, the function:
#' \itemize{
#'   \item restricts the input matrix to rows present in both \code{X} and
#'         \code{cl},
#'   \item optionally applies row-wise z-scoring to obtain shape-centric
#'         profiles,
#'   \item computes the centroid as the mean vector across rows, and
#'   \item evaluates per-row agreement with the centroid using a
#'         user-specified QC metric.
#' }
#'
#' Two QC metrics are supported:
#' \itemize{
#'   \item \code{"Pearson R2"}: squared Pearson correlation between each
#'         row and the centroid (higher values indicate stronger agreement),
#'   \item \code{"BC(HD)"}: Bhattacharyya coefficient computed in Hellinger
#'         space (assumes rows represent square-root–transformed
#'         compositions; values lie in \eqn{[0,1]} with higher values
#'         indicating greater similarity).
#' }
#'
#' The function is designed to be reusable across contexts (e.g. block-level
#' centroid summaries, site-level clustering quality control) and performs
#' no bookkeeping beyond cluster-local statistics.
#'
#' @param X
#' Numeric matrix with rows corresponding to features and columns to ordered
#' measurements (e.g. spline points, time points, or signature components).
#' Row names must uniquely identify features.
#'
#' @param cl
#' Named integer vector of cluster assignments. Names must correspond to
#' row names of \code{X} (or a subset thereof); values are cluster labels.
#'
#' @param qc_method
#' Character scalar specifying the QC metric to use. Either
#' \code{"Pearson R2"} or \code{"BC(HD)"}.
#'
#' @param center_scale_rows
#' Logical scalar; if \code{TRUE}, each row of \code{X} is z-scored
#' (centered and scaled) prior to centroid computation and QC evaluation.
#' If \code{FALSE}, centroids and QC metrics are computed on the original
#' scale of \code{X}.
#'
#' @param require
#' Character scalar controlling behavior when there is no overlap between
#' \code{names(cl)} and \code{rownames(X)}. If \code{"error"} (default), an
#' error is raised; if \code{"intersection"}, an empty result is returned.
#'
#' @return
#' A data frame with one row per cluster and the following columns:
#' \describe{
#'   \item{\code{cluster}}{Cluster label.}
#'   \item{\code{n_used}}{Number of rows (features) used to compute the
#'         centroid for this cluster.}
#'   \item{\code{qc_method}}{QC metric used
#'         (\code{"Pearson R2"} or \code{"BC(HD)"}).}
#'   \item{\code{mean_qc}}{Mean QC value across rows in the cluster.}
#'   \item{\code{sd_qc}}{Standard deviation of the QC values across rows.}
#'   \item{\code{qc_member}}{List-column; each entry is a named numeric
#'         vector of per-row QC values (names are feature IDs).}
#'   \item{\code{centroid}}{List-column; each entry is a numeric vector
#'         giving the centroid representation (length \code{ncol(X)}).}
#' }
#'
.compute_cluster_centroids_qc <- function(
        X,
        cl,
        qc_method = c("Pearson R2", "BC(HD)"),
        center_scale_rows = FALSE,
        require = c("intersection", "error")
) {
    qc_method <- match.arg(qc_method)
    require <- match.arg(require)
    
    common <- intersect(names(cl), rownames(X))
    if (length(common) == 0L) {
        if (require == "error") {
            stop_call_false("No overlap between names(cl) and rownames(X).")
        }
        return(data.frame(
            cluster   = integer(0),
            n_used    = integer(0),
            qc_method = character(0),
            mean_qc   = numeric(0),
            sd_qc     = numeric(0),
            qc_member = I(list()),
            centroid  = I(list()),
            stringsAsFactors = FALSE
        ))
    }
    
    X2  <- X[common, , drop = FALSE]
    cl2 <- cl[common]
    clusters <- sort(unique(as.integer(cl2)))
    
    rows <- vector("list", length(clusters))
    for (i in seq_along(clusters)) {
        c <- clusters[i]
        idx <- which(cl2 == c)
        feats <- names(cl2)[idx]
        Xc <- X2[feats, , drop = FALSE]
        n_used <- nrow(Xc)
        
        if (n_used == 0L) {
            centroid  <- rep(NA_real_, ncol(X2))
            qc_member <- setNames(numeric(0), character(0))
            mean_qc <- NA_real_
            sd_qc <- NA_real_
        } else {
            if (isTRUE(center_scale_rows)) {
                X_use <- t(scale(t(Xc)))
                X_use[is.na(X_use)] <- 0
            } else {
                X_use <- Xc
            }
            
            if (n_used == 1L) {
                centroid <- as.numeric(X_use[1L, ])
                
                if (qc_method == "Pearson R2") {
                    qc_member <- setNames(1, rownames(X_use))
                    mean_qc <- 1
                    sd_qc <- 0
                } else {
                    qc_member <- setNames(1, rownames(X_use))
                    mean_qc <- 1
                    sd_qc <- 0
                }
            } else {
                centroid <- colMeans(X_use, na.rm = TRUE)
                
                if (qc_method == "Pearson R2") {
                    qc_vec <- apply(
                        X_use,
                        1L,
                        function(x) {
                            r <- stats::cor(
                                x,
                                centroid
                                )
                            r^2
                        }
                    )
                } else {
                    # Bhattacharyya coefficient in [0, 1] for Hellinger space
                    # If X_use rows are sqrt(p) and centroid is their mean,
                    # BC can be computed as sum_i x_i * centroid_i.
                    qc_vec <- apply(
                        X_use,
                        1L,
                        function(x) {
                            sum(x * centroid)
                        }
                    )
                    qc_vec <- pmax(0, pmin(1, qc_vec))
                }
                
                names(qc_vec) <- rownames(X_use)
                qc_member <- qc_vec
                mean_qc <- mean(qc_vec, na.rm = TRUE)
                sd_qc   <- stats::sd(qc_vec, na.rm = TRUE)
            }
        }
        
        rows[[i]] <- data.frame(
            cluster   = c,
            n_used    = n_used,
            qc_method = qc_method,
            mean_qc   = mean_qc,
            sd_qc     = sd_qc,
            qc_member = I(list(qc_member)),
            centroid  = I(list(centroid)),
            stringsAsFactors = FALSE
        )
    }
    
    do.call(rbind, rows)
}


#' Normalize a modality-specific gene feature matrix
#'
#' @noRd
#'
#' @description
#' Internal helper that applies modality-appropriate normalization to a
#' gene-level feature matrix prior to multi-omics integration.
#'
#' For one-to-one modalities (e.g. transcriptomics, proteomics), features
#' are interpreted as temporal trajectories and normalization is performed
#' by row-wise z-scoring (shape-centric normalization), ensuring that
#' clustering is driven by relative temporal patterns rather than absolute
#' magnitudes.
#'
#' For many-to-one modalities (e.g. site- or probe-level signatures),
#' normalization is performed using a Hellinger transform (square root of
#' non-negative values), yielding a geometry appropriate for Euclidean
#' methods when working with compositional or fractional signatures.
#'
#' @param mat
#' Numeric matrix of dimension \code{genes x features}. Row names are gene
#' identifiers; columns represent modality-specific features (e.g. spline
#' points or signature components).
#'
#' @param is_many_to_one
#' Logical scalar indicating whether the modality represents a many-to-one
#' gene summary. If \code{TRUE}, Hellinger normalization is applied; if
#' \code{FALSE}, row-wise z-scoring is used.
#'
#' @return
#' A numeric matrix of the same dimension as \code{mat}, normalized according
#' to the modality type. For one-to-one modalities, rows have mean zero and
#' unit variance (with constant rows mapped to zero). For many-to-one
#' modalities, entries are non-negative and square-root transformed.
#'
.normalize_modality_mat <- function(
        mat,
        is_many_to_one
        ) {
    mat <- as.matrix(mat)
    
    if (isTRUE(is_many_to_one)) {
        # Hellinger (assumes non-negative; safe guard)
        return(sqrt(pmax(mat, 0)))
    }
    
    # one-to-one: row-wise z-score (shape-centric)
    z <- t(scale(t(mat)))
    z[is.na(z)] <- 0
    z
}