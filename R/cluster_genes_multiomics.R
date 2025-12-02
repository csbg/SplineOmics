#' Gene-Centric Multi-Omics Clustering Across Blocks and Layers
#'
#' @description
#' Performs gene-centric clustering of multi-omics time-series data across
#' multiple _blocks_ (e.g., time effect, interaction effect) and _layers_
#' (e.g., transcript, protein, feature-level layers with many-to-one gene
#' mapping). The function integrates multiple omics layers within each block by
#' computing layer-specific gene–gene distances, combining them via user-defined
#' weights, and clustering genes based on the resulting unified distance matrix.
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
#' \emph{layers} (e.g., \code{rna}, \code{protein}, \code{phospho}), each being
#' a numeric matrix of dimension \code{features x spline_points}.
#'
#' For gene-level layers, rows represent genes.  
#' For many-to-one layers (e.g., phospho sites, probes), rows represent features
#' mapped externally to genes and summarized into pattern signatures based on
#' the metadata tables.
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
#' @param layer_meta
#' A data frame containing \emph{layer-level metadata}. One row per (block ×
#' layer). Must include:
#' \describe{
#'   \item{\code{block}}{Block identifier linking to \code{block_meta}.}
#'   \item{\code{layer}}{Layer name within the block.}
#'   \item{\code{layer_k}}{Number of pattern clusters to use for building
#'   pattern signatures for many-to-one layers. \code{NA} for layers that are
#'   already gene-level.}
#'   \item{\code{layer_w}}{Relative weight of this layer in the block-specific
#'   distance integration. Weights are normalized to sum to 1 within each
#'   block.}
#' }
#'
#' @param gene_mode
#' Character string specifying how genes should be harmonized across layers
#' within each block prior to clustering.  
#' \describe{
#'   \item{\code{"intersection"}}{Retain only genes present in \emph{all}
#'   layers of the block. Produces the most interpretable multi-omics clusters.}
#'
#'   \item{\code{"union"}}{Retain genes present in \emph{any} layer of the
#'   block. Gene–gene distances are computed using only shared layers per gene
#'   pair, with weights renormalized accordingly. Increases coverage but results
#'   in heterogeneous information across genes.}
#' }
#'
#' @return
#' A tibble \code{cluster_table} with one row per gene and columns containing
#' all block-specific clustering results in a format suitable for downstream
#' enrichment analyses.  
#' The table includes:
#'
#' \describe{
#'   \item{\code{gene}}{Gene identifier used for clustering.}
#'
#'   \item{\code{feature_nr}}{Optional numeric identifier if available in the
#'   input.}
#'
#'   \item{\code{feature_name}}{Optional human-readable feature/gene name.}
#'
#'   \item{\code{cluster_<cond>}}{Integer cluster assignment for each
#'    time-effect block (result category 1). One column per condition, derived
#'    from \code{block_meta$cond1}.}
#'
#'   \item{\code{cluster_cat3_<cond1>_vs_<cond2>}}{Character cluster label for
#'   each interaction block (result category 3). Column names follow the pattern
#'   \code{"cluster_cat3_<cond1>_vs_<cond2>"}. Values are cluster assignments
#'   for genes clustered in the corresponding block, and \code{NA} for genes
#'   not included in that block.}
#' }
#'
#' All clustering columns contain \code{NA} for genes that were not included in
#' the given block (e.g., due to \code{gene_mode = "intersection"} or missing
#' layers in a block).
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
#' layer_meta <- data.frame(
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
#'     layer_meta = layer_meta,
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
        layer_meta,
        gene_mode = "intersection"
) {
    .check_cluster_genes_multiomics_input(
        blocks     = blocks,
        block_meta = block_meta,
        layer_meta = layer_meta,
        gene_mode  = gene_mode
    )
    
    gene_mode <- match.arg(gene_mode, c("intersection", "union"))
    
    block_ids <- unique(block_meta$block)
    block_clusters <- vector("list", length(block_ids))
    names(block_clusters) <- block_ids
    
    centroid_list <- vector("list", length(block_ids))
    names(centroid_list) <- block_ids
    
    for (b in block_ids) {
        meta_b_block <- block_meta[block_meta$block == b, , drop = FALSE]
        meta_b_layer <- layer_meta[layer_meta$block == b, , drop = FALSE]
        
        block_k <- meta_b_block$block_k[1L]
        
        layer_names <- meta_b_layer$layer
        layer_weights <- meta_b_layer$layer_w
        
        layer_mats_raw <- lapply(
            layer_names,
            function(ln) blocks[[b]][[ln]]
        )
        names(layer_mats_raw) <- layer_names
        
        layer_mats <- vector("list", length(layer_names))
        names(layer_mats) <- layer_names
        
        for (i in seq_along(layer_names)) {
            ln <- layer_names[i]
            lk <- meta_b_layer$layer_k[i]
            mat_raw <- layer_mats_raw[[ln]]
            
            if (is.na(lk)) {
                layer_mats[[ln]] <- mat_raw
            } else {
                rn <- rownames(mat_raw)
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
        
        if (gene_mode == "intersection") {
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
        } else {
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
                
                mat_l <- layer_mats[[ln]]
                genes_l <- rownames(mat_l)
                
                idx <- match(genes_l, genes_block)
                
                D_l <- .compute_layer_distance(mat_l)
                D_block[idx, idx] <- D_block[idx, idx] + w * D_l
                W_block[idx, idx] <- W_block[idx, idx] + w
            }
            
            valid <- W_block > 0
            D_block[valid] <- D_block[valid] / W_block[valid]
            diag(D_block) <- 0
        }
        
        cl_b <- .cluster_with_kmeans(
            dist_mat = D_block,
            k        = block_k
        )
        block_clusters[[b]] <- cl_b
        
        # centroids for this block, based on gene-level layer_mats
        centroid_info_b <- .compute_block_centroids(
            block_id    = b,
            layer_mats  = layer_mats,
            cl_b        = cl_b,
            meta_b_layer = meta_b_layer
        )
        centroid_list[[b]] <- centroid_info_b
    }
    
    genes_all <- sort(
        unique(unlist(lapply(block_clusters, names)))
    )
    cluster_table <- data.frame(
        gene = genes_all,
        stringsAsFactors = FALSE
    )
    
    for (b in block_ids) {
        cl_b <- block_clusters[[b]]
        if (is.null(cl_b)) next
        
        mb <- block_meta[block_meta$block == b, , drop = FALSE][1L, ]
        
        rc   <- mb$result_category
        c1   <- mb$cond1
        c2   <- mb$cond2
        
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
    
    list(
        cluster_table = tibble::as_tibble(cluster_table),
        centroid_info = tibble::as_tibble(centroid_info)
    )
}


# Level 1 function definitions -------------------------------------------------


.check_cluster_genes_multiomics_input <- function(
        blocks,
        block_meta,
        layer_meta,
        gene_mode
) {
    # gene_mode
    if (!is.character(gene_mode) || length(gene_mode) != 1L) {
        stop(
            "`gene_mode` must be a single character string ",
            "('intersection' or 'union').",
            call. = FALSE
        )
    }
    gene_mode <- match.arg(gene_mode, c("intersection", "union"))
    
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
    
    # block_meta: blocks exist in blocks
    unknown_blocks <- setdiff(block_meta$block, names(blocks))
    if (length(unknown_blocks) > 0L) {
        stop(
            "Blocks in `block_meta$block` not present in `blocks`: ",
            paste(unique(unknown_blocks), collapse = ", "),
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
    
    # layer_meta: structure
    if (!is.data.frame(layer_meta)) {
        stop(
            "`layer_meta` must be a data.frame (or tibble) with ",
            "layer-level metadata.",
            call. = FALSE
        )
    }
    required_layer_cols <- c("block", "layer", "layer_k", "layer_w")
    missing_layer_cols <- setdiff(required_layer_cols, colnames(layer_meta))
    if (length(missing_layer_cols) > 0L) {
        stop(
            "`layer_meta` is missing required columns: ",
            paste(missing_layer_cols, collapse = ", "),
            call. = FALSE
        )
    }
    
    # each (block, layer) in layer_meta exists in blocks
    for (i in seq_len(nrow(layer_meta))) {
        b <- layer_meta$block[i]
        l <- layer_meta$layer[i]
        
        if (!b %in% names(blocks)) {
            stop(
                "Row ", i, " of `layer_meta` refers to block '", b,
                "', which is not present in `blocks`.",
                call. = FALSE
            )
        }
        if (!l %in% names(blocks[[b]])) {
            stop(
                "Row ", i, " of `layer_meta` refers to layer '", l,
                "' in block '", b, "', which is not present in ",
                "`blocks[[\"",
                b, "\"]]`.",
                call. = FALSE
            )
        }
    }
    
    # each block in block_meta must appear in layer_meta
    blocks_without_layers <- setdiff(
        block_meta$block,
        unique(layer_meta$block)
    )
    if (length(blocks_without_layers) > 0L) {
        stop(
            "Blocks in `block_meta` with no corresponding rows in ",
            "`layer_meta`: ",
            paste(blocks_without_layers, collapse = ", "),
            call. = FALSE
        )
    }
    
    # layer_k: NA or positive numeric
    if (!is.numeric(layer_meta$layer_k)) {
        stop(
            "`layer_meta$layer_k` must be numeric (positive integers) ",
            "or NA.",
            call. = FALSE
        )
    }
    if (any(layer_meta$layer_k[!is.na(layer_meta$layer_k)] <= 0)) {
        stop(
            "`layer_meta$layer_k` must be positive where it is not NA.",
            call. = FALSE
        )
    }
    
    # layer_w: numeric, non-negative; per-block sum > 0
    if (!is.numeric(layer_meta$layer_w)) {
        stop("`layer_meta$layer_w` must be numeric.", call. = FALSE)
    }
    if (any(is.na(layer_meta$layer_w))) {
        stop(
            "`layer_meta$layer_w` contains NA. Please provide weights ",
            "for all rows.",
            call. = FALSE
        )
    }
    if (any(layer_meta$layer_w < 0)) {
        stop(
            "`layer_meta$layer_w` must be non-negative.",
            call. = FALSE
        )
    }
    
    w_by_block <- tapply(layer_meta$layer_w, layer_meta$block, sum)
    if (any(w_by_block <= 0)) {
        bad_blocks <- names(w_by_block)[w_by_block <= 0]
        stop(
            "For each block, the sum of `layer_meta$layer_w` must be > 0. ",
            "Blocks with invalid weights: ",
            paste(bad_blocks, collapse = ", "),
            call. = FALSE
        )
    }
    
    # consistency: blocks in layer_meta appear in block_meta
    unknown_in_block_meta <- setdiff(
        unique(layer_meta$block),
        block_meta$block
    )
    if (length(unknown_in_block_meta) > 0L) {
        stop(
            "Blocks in `layer_meta$block` not present in ",
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
    
    # ensure (block, layer) pairs in layer_meta are unique
    bl_pairs <- paste(layer_meta$block, layer_meta$layer, sep = "||")
    if (any(duplicated(bl_pairs))) {
        dup <- unique(bl_pairs[duplicated(bl_pairs)])
        stop(
            "Duplicate (block, layer) combinations in `layer_meta`: ",
            paste(dup, collapse = ", "),
            ". Each (block, layer) must appear only once.",
            call. = FALSE
        )
    }

    # per-layer checks using meta: NA, one-to-one vs many-to-one
    for (i in seq_len(nrow(layer_meta))) {
        b  <- layer_meta$block[i]
        l  <- layer_meta$layer[i]
        lk <- layer_meta$layer_k[i]
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
    for (b in unique(layer_meta$block)) {
        rows_b <- layer_meta$block == b &
            is.na(layer_meta$layer_k)
        layers_gene <- layer_meta$layer[rows_b]
        
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


#' Compute centroids and quality metrics for one block
#'
#' @noRd
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
                mean_R2 <- NA_real_
                sd_R2 <- NA_real_
            } else {
                X <- mat_l[genes_layer, , drop = FALSE]
                
                # row-wise z-score (shape-centric)
                X_z <- t(scale(t(X)))
                
                if (n_used == 1L) {
                    centroid <- as.numeric(X_z[1L, ])
                    r2_vec <- 1
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
                    
                    mean_R2 <- mean(r2_vec, na.rm = TRUE)
                    sd_R2 <- stats::sd(r2_vec, na.rm = TRUE)
                }
            }
            
            row_i <- row_i + 1L
            centroid_rows[[row_i]] <- data.frame(
                block           = block_id,
                layer           = ln,
                layer_type      = layer_type,
                cluster         = c,
                n_genes_cluster = n_cluster,
                n_genes_used    = n_used,
                coverage        = coverage,
                mean_R2         = mean_R2,
                sd_R2           = sd_R2,
                centroid        = I(list(centroid)),
                stringsAsFactors = FALSE
            )
        }
    }
    
    if (length(centroid_rows) == 0L) {
        centroid_info_b <- data.frame(
            block           = character(0),
            layer           = character(0),
            layer_type      = character(0),
            cluster         = integer(0),
            n_genes_cluster = integer(0),
            n_genes_used    = integer(0),
            coverage        = numeric(0),
            mean_R2         = numeric(0),
            sd_R2           = numeric(0),
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
