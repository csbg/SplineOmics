#' Gene-centric multi-omics clustering across data modalities
#'
#' @description
#' Performs gene-centric clustering of multi-omics feature representations that
#' are supplied as modality-specific matrices. The function harmonizes genes
#' across modalities, constructs a gene-centric joint feature matrix, computes a
#' UMAP neighborhood graph using \pkg{uwot}, and performs spectral clustering on
#' the resulting UMAP fuzzy graph (\code{fgraph}).
#'
#' The input matrices are assumed to be precomputed gene-level trajectories
#' (e.g., spline values or coefficients) or many-to-one feature-level matrices
#' (e.g., phospho sites) that are summarized into gene-level pattern signature
#' vectors inside this function.
#'
#' This function performs statistical computation and returns the UMAP fit for
#' downstream visualization. Clustering is performed on the UMAP graph rather
#' than on the two-dimensional embedding coordinates.
#'
#' @param data
#' Named list defining the multi-condition, multi-modality input data.
#' The outer list corresponds to experimental conditions (e.g. control, and
#' treatment). Each element of the outer list must itself be a
#' named list of numeric matrices, one per data modality.
#'
#' For each condition, the inner list entries represent modalities and must
#' share the same modality names across conditions. Each modality matrix must
#' have row names and can be of two types:
#'
#' \describe{
#'   \item{One-to-one (gene-level) modality}{Rows represent genes directly.
#'   Row names must be gene identifiers in the format: \code{<gene_id>}.
#'   Columns represent the
#'   modality-specific representation used for clustering, such as spline
#'   values at fixed time points, spline coefficients, or other numeric
#'   features.}
#'
#'   \item{Many-to-one (feature-level) modality}{Rows represent features that
#'   map to genes (e.g., phospho sites, probes). Row names must follow the
#'   pattern \code{<gene_id>_<feature_id>} where the gene identifier precedes
#'   the first underscore. Such modalities are aggregated into gene-level
#'   pattern signature vectors using an internal feature clustering step prior
#'   to integration.}
#' }
#' 
#' The <gene_id> parts of the rownames of all modalities across all conditions
#' should match, otherwise, the gene-centric clustering is not possible!
#' Modality matrices may differ in their number of columns both within and
#' across conditions. After aggregation (for many-to-one modalities) and
#' normalization,
#' all condition- and modality-specific gene-level matrices are aligned
#' according to \code{gene_mode} and concatenated to form the gene-centric
#' feature matrix used to construct the UMAP graph and clustering.
#'
#' @param meta
#' Data frame with one row per modality, providing modality-level parameters
#' required for gene-centric representation building. The data frame must
#' contain at least the following columns:
#'
#' \describe{
#'   \item{\code{modality}}{
#'     Character scalar giving the modality identifier. Each value must match
#'     a modality name present in \code{data[[condition]]}. The order of rows
#'     defines the modality ordering used throughout downstream processing.}
#'
#'   \item{\code{many_to_one_k}}{
#'     Integer or \code{NA}. If \code{NA}, the modality is treated as one-to-one
#'     (gene-level) and passed through unchanged. If an integer, the modality
#'     is treated as many-to-one (feature-level) and collapsed to gene-level
#'     signatures using \code{many_to_one_k} global archetypes.}
#'
#'   \item{\code{modality_w}}{
#'     Numeric, non-negative scalar giving the relative weight of the modality
#'     in the joint feature space. Weights are normalized internally to sum to
#'     one across modalities before being applied.}
#' }
#'
#' Additional columns may be present but are ignored.
#'
#' @param k
#' Integer. Number of gene clusters for spectral clustering.
#'
#' @param gene_mode
#' Character string specifying how to harmonize genes across modalities prior
#' to constructing the joint feature matrix:
#' \describe{
#'   \item{\code{"intersection"}}{Retain only genes present in all modalities.}
#'   \item{\code{"union"}}{Retain genes present in any modality. Gene vectors
#'   are constructed using available modalities; missing modality blocks are
#'   handled internally.}
#' }
#'
#' @param n_neighbors
#' Integer. Size of the local neighborhood used by UMAP to construct the
#' k-nearest-neighbor graph. Larger values emphasize broader structure, smaller
#' values emphasize local structure.
#'
#' @param verbose
#' Logical scalar indicating whether progress messages should be emitted via
#' \code{rlang::inform()} by internal helpers.
#'
#' @return
#' A named list with the following elements:
#' \describe{
#'   \item{\code{cluster_table}}{A tibble with one row per gene and columns
#'   \code{gene} and \code{cluster}.}
#'
#'   \item{\code{centroid_info}}{A tibble summarizing modality-specific cluster
#'   centroids and within-cluster QC. The centroid representation is stored as a
#'   list-column.}
#'
#'   \item{\code{many_to_one_clustering_qc}}{A tibble (or \code{NULL}) with QC
#'   diagnostics for many-to-one feature clustering steps, if any are present.}
#'
#'   \item{\code{umap_fit}}{The object returned by \code{uwot::umap()},
#'   including \code{$embedding} (genes \eqn{\times} \code{n_components}) and,
#'   if requested via \code{ret_extra}, \code{$fgraph} (the UMAP fuzzy graph
#'   used for spectral clustering).}
#' }
#'
#' @examples
#' set.seed(1)
#' genes <- paste0("gene", 1:8)
#'
#' rna <- matrix(
#'   rnorm(length(genes) * 5),
#'   nrow = length(genes),
#'   dimnames = list(genes, NULL)
#' )
#'
#' prot <- matrix(
#'   rnorm(length(genes) * 3),
#'   nrow = length(genes),
#'   dimnames = list(genes, NULL)
#' )
#'
#' data <- list(
#'   condition1 = list(rna = rna, protein = prot)
#' )
#'
#' meta <- data.frame(
#'   modality = c("rna", "protein"),
#'   many_to_one_k = c(NA_real_, NA_real_),
#'   modality_w = c(1, 1),
#'   stringsAsFactors = FALSE
#' )
#'
#' res <- cluster_genes_multiomics(
#'   data = data,
#'   meta = meta,
#'   k = 3L,
#'   gene_mode = "intersection",
#'   n_neighbors = 5L
#' )
#'
#' res$cluster_table
#'
#' @export
#' 
cluster_genes_multiomics <- function(
        data,
        meta,
        k,
        gene_mode = c("intersection", "union"),
        n_neighbors = 15L,
        verbose = FALSE
) {
    start_time <- Sys.time()
    .check_cluster_genes_multiomics_dependencies()
    gene_mode <- match.arg(
        gene_mode,
        c("intersection", "union")
        )
    
    .check_cluster_genes_multiomics_input(
        data        = data,
        meta        = meta,
        k           = k,
        n_neighbors = n_neighbors,
        verbose     = verbose
    )

    rep <- .build_gene_centric_representation(
        data      = data,
        meta      = meta,
        gene_mode = gene_mode,
        verbose   = verbose
    )

    fit <- uwot::umap(
        X            = rep$X,
        n_neighbors  = n_neighbors,
        n_components = 2,
        ret_model    = TRUE,
        ret_extra    = "fgraph"
    )

    cl <- .spectral_clustering(
        W = fit$fgraph,
        k = k
    )
    names(cl) <- rep$gene_ids

    centroid_info <- .compute_centroids(
        mats_norm = rep$mats_norm,
        cl        = cl,
        meta      = rep$aligned_meta
    )

    centroid_info <- .rank_clusters(centroid_info)

    cluster_table <- tibble::tibble(
        gene    = rep$gene_ids,
        cluster = as.integer(cl)
    )

    if (isTRUE(verbose)) {
        end_time <- Sys.time()
        elapsed <- end_time - start_time
        rlang::inform(c(
            "[cluster_genes_multiomics] total runtime: ",
            format(elapsed, digits = 2)
        ))
    }
    
    list(
        cluster_table = cluster_table,
        centroid_info = tibble::as_tibble(centroid_info),
        many_to_one_clustering_qc =
            if (is.null(rep$many_to_one_clustering_qc)) {
                NULL
            } else {
                tibble::as_tibble(rep$many_to_one_clustering_qc)
            },
        umap_fit = fit
    )
}


# Level 1 function definitions -------------------------------------------------


#' Validate availability of suggested package dependencies
#'
#' @noRd
#'
#' @description
#' Internal helper that checks whether all suggested packages required for
#' UMAP-based clustering and spectral graph clustering are installed.
#'
#' This function verifies the availability of the following packages:
#'
#' \itemize{
#'   \item{\pkg{uwot} for UMAP neighborhood graph construction.}
#'   \item{\pkg{Matrix} for sparse matrix representations and operations.}
#'   \item{\pkg{RSpectra} for sparse eigenvalue decomposition used in spectral
#'   clustering.}
#' }
#'
#' These packages are declared under \code{Suggests} rather than
#' \code{Imports} because they are only required when this functionality is
#' invoked. If one or more packages are missing, execution is aborted with an
#' informative error message via \code{rlang::abort()}.
#'
#' @return
#' Invisibly returns \code{TRUE} if all required packages are available.
#' Otherwise, an error is raised.
#' 
.check_cluster_genes_multiomics_dependencies <- function() {
    missing <- character()
    
    if (!requireNamespace("uwot", quietly = TRUE)) {
        missing <- c(missing, "uwot")
    }
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        missing <- c(missing, "Matrix")
    }
    if (!requireNamespace("RSpectra", quietly = TRUE)) {
        missing <- c(missing, "RSpectra")
    }
    
    if (length(missing) > 0) {
        rlang::abort(c(
            "Required suggested packages are not installed.",
            paste0(
                "Missing package(s): ",
                paste(missing, collapse = ", ")
            ),
            "Please install them with:",
            paste0(
                "install.packages(c(",
                paste(sprintf('"%s"', missing), collapse = ", "),
                "))"
            )
        ))
    }
    
    invisible(TRUE)
}


#' Validate inputs for gene-centric multi-omics clustering
#'
#' @noRd
#'
#' @description
#' Internal helper that performs interface-level validation for
#' \code{cluster_genes_multiomics()}. This function checks that:
#'
#' \itemize{
#'   \item{\code{data} follows the expected nested structure
#'   (condition × modality) and contains valid numeric matrices.}
#'   \item{\code{meta} contains the required modality-level metadata columns
#'   and valid entries.}
#'   \item{Core scalar parameters (\code{k}, \code{n_pcs}, \code{knn_k}) are
#'   well-formed integers in valid ranges.}
#'   \item{Neighbor-search and graph parameters (\code{bn_param}, \code{sigma})
#'   are valid.}
#'   \item{\code{data} and \code{meta} agree on the set of modalities.}
#'   \item{\code{verbose} is a valid logical scalar.}
#' }
#'
#' All user-facing errors are raised via \code{rlang::abort()} to provide
#' informative, structured messages. Downstream internal helpers rely on
#' these checks and therefore assume inputs are valid.
#'
#' @param data
#' Named list of conditions. Each condition is a named list of numeric
#' modality matrices. See \code{cluster_genes_multiomics()} for the expected
#' structure.
#'
#' @param meta
#' Data frame with one row per modality and required columns
#' \code{modality}, \code{many_to_one_k}, and \code{modality_w}.
#'
#' @param k
#' Integer scalar giving the number of gene clusters for spectral clustering.
#'
#' @param n_neighbors
#' Integer. Size of the local neighborhood used by UMAP to construct the
#' k-nearest-neighbor graph. Larger values emphasize broader structure, smaller
#' values emphasize local structure.
#'
#' @param verbose
#' Logical scalar indicating whether informative messages should be emitted.
#'
#' @return
#' Invisibly returns \code{TRUE} if all checks pass. Otherwise, raises an
#' error.
#' 
.check_cluster_genes_multiomics_input <- function(
        data,
        meta,
        k,
        n_neighbors,
        verbose
) {
    .check_cluster_genes_multiomics_data(data)
    .check_cluster_genes_multiomics_meta(meta)
    
    if (!is.numeric(k) || length(k) != 1L || is.na(k) ||
        k != floor(k) || k < 2L) {
        rlang::abort(c(
            "Invalid `k`.",
            "i" = "`k` must be a single integer >= 2."
        ))
    }

    .check_cluster_genes_multiomics_cross_args(
        data = data,
        meta = meta
    )
    
    if (!is.numeric(n_neighbors) ||
        length(n_neighbors) != 1L ||
        is.na(n_neighbors) ||
        n_neighbors != floor(n_neighbors) ||
        n_neighbors < 2L) {
        rlang::abort(c(
            "Invalid `n_neighbors`.",
            "i" = "`n_neighbors` must be a single integer >= 2."
        ))
    }
    
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        rlang::abort(c(
            "Invalid `verbose`.",
            "i" = "`verbose` must be a single non-missing logical value."
        ))
    }
    
    invisible(TRUE)
}


#' Build a gene-centric joint feature representation across conditions
#'
#' @noRd
#'
#' @description
#' Internal helper that converts nested multi-condition, multi-modality input
#' data into a single **gene-centric joint feature matrix** suitable for PCA
#' and downstream graph-based clustering.
#'
#' The helper performs the following conceptual steps:
#'
#' \itemize{
#'   \item{Align modality metadata to the modalities present in \code{data}.}
#'   \item{Collapse many-to-one modalities (e.g. sites, probes) into
#'   gene-level pattern signatures using \code{many_to_one_k}.}
#'   \item{Normalize each gene-level modality matrix using
#'   \code{.normalize_modality_mat()} (z-scoring for one-to-one, Hellinger for
#'   many-to-one signatures).}
#'   \item{Harmonize the gene universe across all condition × modality blocks
#'   according to \code{gene_mode} (intersection or union).}
#'   \item{Align, weight, and concatenate all blocks to form the joint matrix
#'   \code{X} with one row per gene.}
#' }
#'
#' The returned \code{X} matrix is the direct input to PCA. Additional return
#' components provide provenance (feature mapping, aligned metadata) and
#' access to the normalized block matrices for centroid/QC computation.
#'
#' @param data
#' Named list of conditions. Each condition is a named list of numeric
#' matrices, one per modality. One-to-one modalities must use gene IDs as
#' row names. Many-to-one modalities must use row names of the form
#' \code{<gene_id>_<feature_id>} so that features can be mapped back to genes.
#'
#' @param meta
#' Data frame with one row per modality describing modality-level parameters,
#' including \code{modality}, \code{many_to_one_k}, and \code{modality_w}.
#' This metadata is condition-independent.
#'
#' @param gene_mode
#' Character string specifying how genes are harmonized across all blocks.
#' \code{"intersection"} retains only genes present in all condition × modality
#' blocks. \code{"union"} retains genes present in at least one block and
#' fills missing blocks with zeros during concatenation.
#'
#' @param verbose
#' Logical scalar indicating whether informative messages should be emitted.
#' Currently unused in this helper, but included for a consistent internal
#' interface.
#'
#' @return
#' A list with the following components:
#'
#' \describe{
#'   \item{\code{X}}{
#'     Numeric matrix (\code{genes x features}) giving the gene-centric joint
#'     feature representation obtained by weighting and concatenating all
#'     condition × modality blocks. Row names are gene IDs. Column names
#'     encode condition and modality provenance.}
#'
#'   \item{\code{gene_ids}}{
#'     Character vector of gene IDs defining the row order used across all
#'     downstream objects derived from \code{X}.}
#'
#'   \item{\code{feature_index}}{
#'     Data frame mapping each column of \code{X} to its originating condition,
#'     modality, and within-block feature name (and the modality weight used).}
#'
#'   \item{\code{mats_norm}}{
#'     Nested list of normalized gene-level matrices, structured as
#'     \code{mats_norm[[condition]][[modality]]}. This is used for downstream
#'     centroid and QC summaries per condition and modality.}
#'
#'   \item{\code{aligned_meta}}{
#'     Modality metadata aligned to the modalities present in \code{data},
#'     including normalized modality weights and many-to-one flags.}
#'
#'   \item{\code{many_to_one_clustering_qc}}{
#'     Data frame of QC diagnostics for many-to-one feature clustering steps,
#'     or \code{NULL} if no many-to-one modalities are present.}
#' }
#' 
.build_gene_centric_representation <- function(
        data,
        meta,
        gene_mode,
        verbose
) {
    cond_ids <- names(data)
    mod_ids <- names(data[[1L]])
    
    aligned_meta <- .align_modality_meta(
        meta = meta,
        modality_ids = mod_ids
    )
    
    res_gene_level <- .collapse_many_to_one_modalities(
        data = data,
        meta = aligned_meta,
        verbose = verbose
    )
    
    mats_gene <- res_gene_level$mats_gene
    many_to_one_qc <- res_gene_level$many_to_one_clustering_qc
    
    mats_norm <- .normalize_gene_level_modalities(
        mats_gene = mats_gene,
        meta = aligned_meta,
        verbose = verbose
    )
    
    gene_ids <- .harmonize_genes_across_blocks(
        mats_norm = mats_norm,
        gene_mode = gene_mode
    )
    
    X_res <- .flatten_blocks_to_joint_matrix(
        mats_norm = mats_norm,
        meta = aligned_meta,
        gene_ids = gene_ids,
        verbose = verbose
    )
    
    list(
        X = X_res$X,
        gene_ids = gene_ids,
        feature_index = X_res$feature_index,
        mats_norm = mats_norm,
        aligned_meta = aligned_meta,
        many_to_one_clustering_qc = many_to_one_qc
    )
}


#' Perform spectral clustering on a gene–gene adjacency matrix
#'
#' @noRd
#'
#' @description
#' Internal helper that performs **spectral clustering** on a weighted,
#' undirected gene–gene graph represented by an adjacency matrix.
#'
#' Starting from the adjacency matrix \code{W}, the function constructs the
#' (unnormalized) graph Laplacian \eqn{L = D - W}, where \eqn{D} is the diagonal
#' degree matrix. The first \code{k} eigenvectors corresponding to the smallest
#' eigenvalues of \eqn{L} define a low-dimensional spectral embedding of the
#' genes.
#'
#' Genes are then clustered by applying k-means to this spectral embedding
#' after row-wise normalization. This procedure identifies groups of genes
#' that are tightly connected in the original kNN graph and therefore exhibit
#' similar behavior in the reduced feature space.
#'
#' This implementation follows the standard unnormalized spectral clustering
#' pipeline and is intended for moderate-sized graphs (e.g. thousands of
#' genes), as typically encountered in gene-centric multi-omics analyses.
#'
#' @param W
#' Numeric, symmetric matrix (\code{genes x genes}) representing the weighted
#' adjacency matrix of the gene–gene graph. Entries must be non-negative, with
#' zeros indicating no edge between gene pairs.
#'
#' @param k
#' Integer scalar giving the number of clusters to compute. This also
#' determines the number of eigenvectors retained for the spectral embedding.
#'
#' @return
#' A list with one component:
#'
#' \describe{
#'   \item{\code{membership}}{
#'     Integer vector of length equal to the number of genes, giving the
#'     cluster assignment for each gene. The order corresponds to the row and
#'     column order of \code{W}.}
#' }
#' 
.spectral_clustering <- function(
        W, 
        k
        ) {
    # Keep sparse if possible
    if (!inherits(W, "Matrix")) W <- Matrix::Matrix(W, sparse = TRUE)
    
    # Symmetrize (important; UMAP graphs can be slightly asymmetric)
    W <- (W + Matrix::t(W)) / 2
    
    # Remove self-loops
    Matrix::diag(W) <- 0
    
    # Degrees
    d <- Matrix::rowSums(W)
    if (any(d < 0)) 
        stop(
            "W has negative weights; spectral clustering assumes",
            "nonnegative weights."
            )
    
    # Handle isolated nodes
    d_inv_sqrt <- 1 / sqrt(pmax(d, .Machine$double.eps))
    D_inv_sqrt <- Matrix::Diagonal(x = d_inv_sqrt)
    
    # Normalized Laplacian: L = I - D^{-1/2} W D^{-1/2}
    n <- nrow(W)
    L <- Matrix::Diagonal(n = n, x = 1) - (D_inv_sqrt %*% W %*% D_inv_sqrt)
    
    # Smallest eigenvectors of L
    # (for k clusters you typically want the k smallest eigenvectors)
    eig <- RSpectra::eigs_sym(
        A     = L,
        k     = k,
        which = "SM"
        )
    U <- eig$vectors
    
    # Row-normalize (Ng, Jordan, Weiss style)
    rs <- sqrt(rowSums(U^2))
    rs[rs == 0] <- 1
    U <- U / rs
    
    cl <- stats::kmeans(
        x       = U,
        centers = k,
        nstart  = 10
        )$cluster
    cl <- as.integer(cl)
    names(cl) <- rownames(W)
    
    cl
}


#' Compute modality- and condition-specific cluster centroids and QC summaries
#'
#' @noRd
#'
#' @description
#' Internal helper that summarizes the resulting gene clusters in the original
#' normalized feature spaces for each condition and modality.
#'
#' For every condition × modality block in \code{mats_norm}, the function:
#'
#' \itemize{
#'   \item{Computes a cluster centroid representation for each gene cluster.}
#'   \item{Computes within-cluster coherence (QC) for each member gene relative
#'   to its cluster centroid using \code{.compute_cluster_centroids_qc()}.}
#'   \item{Tracks coverage, i.e. the fraction of genes in each cluster that
#'   are present in the current condition × modality block. This is relevant
#'   when \code{gene_mode = "union"} was used upstream and some genes are
#'   missing from certain blocks.}
#' }
#'
#' The QC method depends on the modality type:
#'
#' \itemize{
#'   \item{\code{"Pearson R2"}} for one-to-one modalities (trajectory-like
#'   features).
#'   \item{\code{"BC(HD)"}} for many-to-one modalities represented as
#'   Hellinger-transformed signature vectors.
#' }
#'
#' The returned table is designed for downstream reporting and ranking of
#' clusters by coherence.
#'
#' @param mats_norm
#' Nested list of normalized gene-level matrices structured as
#' \code{mats_norm[[condition]][[modality]]}. Row names are gene IDs; columns
#' are modality-specific features.
#'
#' @param cl
#' Integer vector of gene cluster assignments. \code{names(cl)} must be gene
#' identifiers.
#'
#' @param meta
#' Data frame of modality metadata aligned to the modalities used in
#' \code{mats_norm}. Must contain at least \code{modality} and
#' \code{many_to_one_k} to determine modality type.
#'
#' @return
#' A data frame with one row per condition, modality, and gene cluster.
#' Columns include:
#'
#' \describe{
#'   \item{\code{condition}}{Condition identifier.}
#'   \item{\code{modality}}{Modality identifier.}
#'   \item{\code{modality_type}}{Either \code{"one_to_one"} or
#'   \code{"many_to_one"}.}
#'   \item{\code{cluster}}{Cluster label.}
#'   \item{\code{n_genes_cluster}}{Number of genes assigned to the cluster.}
#'   \item{\code{n_genes_used}}{Number of genes from the cluster present in
#'   the current block.}
#'   \item{\code{coverage}}{Fraction \code{n_genes_used / n_genes_cluster}.}
#'   \item{\code{qc_method}}{QC method used (\code{"Pearson R2"} or
#'   \code{"BC(HD)"}).}
#'   \item{\code{mean_qc}}{Mean QC value across member genes.}
#'   \item{\code{sd_qc}}{Standard deviation of QC values across member genes.}
#'   \item{\code{qc_member}}{List-column of per-gene QC values.}
#'   \item{\code{centroid}}{List-column of centroid vectors in the modality
#'   feature space.}
#' }
#' 
.compute_centroids <- function(
        mats_norm,
        cl,
        meta
        
) {
    cond_ids <- names(mats_norm)
    mod_ids <- meta$modality
    clusters <- sort(unique(as.integer(cl)))
    out <- list()
    row_i <- 0L
    
    for (cond in cond_ids) {
        mats_cond <- mats_norm[[cond]]
        
        for (m in mod_ids) {
            mat_m <- mats_cond[[m]]
            
            lk <- meta$many_to_one_k[meta$modality == m][1L]
            modality_type <- if (is.na(lk)) "one_to_one" else "many_to_one"
            
            qc_method <- if (modality_type == "one_to_one") {
                "Pearson R2"
            } else {
                "BC(HD)"
            }
            
            stats_m <- .compute_cluster_centroids_qc(
                X = mat_m,
                cl = cl,
                qc_method = qc_method,
                center_scale_rows = FALSE,
                require = "intersection"
            )
            
            for (c in clusters) {
                genes_cluster <- names(cl)[cl == c]
                n_cluster <- length(genes_cluster)
                
                genes_used <- intersect(genes_cluster, rownames(mat_m))
                n_used <- length(genes_used)
                
                coverage <- if (n_cluster > 0L) n_used / n_cluster else NA_real_
                
                st <- stats_m[stats_m$cluster == c, , drop = FALSE]
                
                if (nrow(st) == 0L) {
                    st <- data.frame(
                        cluster = c,
                        n_used = 0L,
                        qc_method = qc_method,
                        mean_qc = NA_real_,
                        sd_qc = NA_real_,
                        qc_member =
                            I(list(setNames(numeric(0), character(0)))),
                        centroid = I(list(rep(NA_real_, ncol(mat_m)))),
                        stringsAsFactors = FALSE
                    )
                }
                
                row_i <- row_i + 1L
                out[[row_i]] <- data.frame(
                    condition = cond,
                    modality = m,
                    modality_type = modality_type,
                    cluster = c,
                    n_genes_cluster = n_cluster,
                    n_genes_used = n_used,
                    coverage = coverage,
                    qc_method = st$qc_method,
                    mean_qc = st$mean_qc,
                    sd_qc = st$sd_qc,
                    qc_member = st$qc_member,
                    centroid = st$centroid,
                    stringsAsFactors = FALSE
                )
            }
        }
    }
    
    if (length(out) == 0L) {
        return(data.frame(
            condition = character(0),
            modality = character(0),
            modality_type = character(0),
            cluster = integer(0),
            n_genes_cluster = integer(0),
            n_genes_used = integer(0),
            coverage = numeric(0),
            qc_method = character(0),
            mean_qc = numeric(0),
            sd_qc = numeric(0),
            qc_member = I(list()),
            centroid = I(list()),
            stringsAsFactors = FALSE
        ))
    }
    
    do.call(rbind, out)
}


# Level 2 function definitions -------------------------------------------------
## .check_cluster_genes_multiomics_input ---------------------------------------


#' Validate multi-condition, multi-modality input data
#'
#' @noRd
#'
#' @description
#' Internal helper that validates the structure and contents of the
#' \code{data} argument passed to \code{cluster_genes_multiomics()}.
#'
#' The expected structure of \code{data} is a nested list:
#'
#' \itemize{
#'   \item{Top level: named list of conditions.}
#'   \item{Second level: for each condition, a named list of modalities.}
#'   \item{Leaf level: numeric matrices with genes or features in rows.}
#' }
#'
#' This function verifies that:
#'
#' \itemize{
#'   \item{\code{data} is a non-empty named list with unique condition names.}
#'   \item{Each condition contains a non-empty named list of modalities.}
#'   \item{Each modality is a numeric matrix with valid dimensions.}
#'   \item{Row names are present, non-empty, unique, and finite-valued.}
#' }
#'
#' All checks are performed at the interface level so that downstream internal
#' helpers may assume well-formed inputs.
#'
#' @param data
#' Named list of conditions. Each condition is a named list of numeric
#' modality matrices. One-to-one and many-to-one modalities are both allowed.
#'
#' @return
#' Invisibly returns \code{TRUE} if all checks pass. Otherwise, raises an
#' error via \code{rlang::abort()}.
#' 
.check_cluster_genes_multiomics_data <- function(data) {
    if (!is.list(data) || length(data) == 0L) {
        rlang::abort(c(
            "Invalid `data`.",
            "i" = paste0(
                "`data` must be a non-empty named list of conditions, ",
                "each containing a named list of modality matrices."
            )
        ))
    }
    
    cond_names <- names(data)
    if (is.null(cond_names) || anyNA(cond_names) || any(cond_names == "")) {
        rlang::abort(c(
            "Invalid `data`.",
            "i" = "`data` must be a named list (no empty names)."
        ))
    }
    
    if (anyDuplicated(cond_names)) {
        rlang::abort(c(
            "Invalid `data`.",
            "i" = "Condition names must be unique.",
            "x" = paste0(
                "Duplicated names: ",
                paste(
                    unique(cond_names[duplicated(cond_names)]),
                    collapse = ", "
                )
            )
        ))
    }
    
    mod_name_list <- lapply(data, names)
    mod_names <- sort(unique(unlist(mod_name_list, use.names = FALSE)))
    
    if (length(mod_names) == 0L) {
        rlang::abort(c(
            "Invalid `data`.",
            "i" = "Each condition must contain a non-empty modality list."
        ))
    }
    
    for (cond in cond_names) {
        mods <- data[[cond]]
        
        if (!is.list(mods) || length(mods) == 0L) {
            rlang::abort(c(
                "Invalid `data`.",
                "i" =
                    "Each condition must contain a non-empty modality list.",
                "x" = paste0("Condition '", cond, "' is empty.")
            ))
        }
        
        mod_names_cond <- names(mods)
        if (is.null(mod_names_cond) ||
            anyNA(mod_names_cond) ||
            any(mod_names_cond == "")) {
            rlang::abort(c(
                "Invalid `data`.",
                "i" = paste0(
                    "Each condition's modality list must be named ",
                    "(no empty names)."
                ),
                "x" = paste0(
                    "Condition '", cond, "' has invalid modality names."
                )
            ))
        }
        
        if (anyDuplicated(mod_names_cond)) {
            rlang::abort(c(
                "Invalid `data`.",
                "i" =
                    "Modality names must be unique within each condition.",
                "x" = paste0(
                    "Condition '", cond,
                    "' has duplicated modalities: ",
                    paste(
                        unique(
                            mod_names_cond[
                                duplicated(mod_names_cond)
                            ]
                        ),
                        collapse = ", "
                    )
                )
            ))
        }
        
        for (m in mod_names_cond) {
            mat <- mods[[m]]
            
            if (!is.matrix(mat) || !is.numeric(mat)) {
                rlang::abort(c(
                    "Invalid `data` element.",
                    "i" = "Each modality must be a numeric matrix.",
                    "x" = paste0(
                        "Condition '", cond, "', modality '", m,
                        "' is not a numeric matrix."
                    )
                ))
            }
            
            rn <- rownames(mat)
            if (is.null(rn) || anyNA(rn) || any(rn == "")) {
                rlang::abort(c(
                    "Invalid `data` element.",
                    "i" =
                        "Each modality matrix must have non-empty row names.",
                    "x" = paste0(
                        "Condition '", cond, "', modality '", m,
                        "' has missing/empty row names."
                    )
                ))
            }
            
            if (nrow(mat) < 2L) {
                rlang::abort(c(
                    "Invalid `data` element.",
                    "i" =
                        "Each modality matrix must have at least 2 rows.",
                    "x" = paste0(
                        "Condition '", cond, "', modality '", m,
                        "' has nrow = ", nrow(mat), "."
                    )
                ))
            }
            
            if (ncol(mat) < 1L) {
                rlang::abort(c(
                    "Invalid `data` element.",
                    "i" =
                        "Each modality matrix must have at least 1 column.",
                    "x" = paste0(
                        "Condition '", cond, "', modality '", m,
                        "' has ncol = ", ncol(mat), "."
                    )
                ))
            }
            
            if (any(!is.finite(mat))) {
                rlang::abort(c(
                    "Invalid `data` element.",
                    "i" =
                        "Modality matrices must contain finite values.",
                    "x" = paste0(
                        "Condition '", cond, "', modality '", m,
                        "' contains NA/NaN/Inf."
                    )
                ))
            }
            
            if (anyDuplicated(rn)) {
                rlang::abort(c(
                    "Invalid `data` element.",
                    "i" =
                        "Row names must be unique within each matrix.",
                    "x" = paste0(
                        "Condition '", cond, "', modality '", m,
                        "' has duplicated row names."
                    )
                ))
            }
        }
    }
    
    invisible(TRUE)
}


#' Validate modality metadata for multi-omics clustering
#'
#' @noRd
#'
#' @description
#' Internal helper that validates the modality-level metadata supplied to
#' \code{cluster_genes_multiomics()}.
#'
#' The \code{meta} data frame defines how individual data modalities are
#' interpreted and weighted during construction of the gene-centric feature
#' representation. This function performs structural and semantic checks to
#' ensure that downstream steps can rely on a consistent and well-defined
#' modality specification.
#'
#' In particular, this function verifies:
#' \itemize{
#'   \item Presence and correctness of required columns.
#'   \item Uniqueness and validity of modality identifiers.
#'   \item Correct specification of many-to-one collapsing parameters.
#'   \item Valid, strictly positive modality weights.
#' }
#'
#' @param meta
#' A data frame with one row per modality and at least the following columns:
#' \describe{
#'   \item{\code{modality}}{
#'     Character vector of unique modality identifiers. Each entry must match
#'     a modality present in \code{data}.
#'   }
#'   \item{\code{many_to_one_k}}{
#'     Numeric vector specifying the number of feature clusters used to collapse
#'     many-to-one modalities to the gene level. Use \code{NA} for one-to-one
#'     (already gene-level) modalities.
#'   }
#'   \item{\code{modality_w}}{
#'     Numeric vector of strictly positive modality weights. Values are
#'     normalized internally and control each modality's contribution to the
#'     joint feature matrix.
#'   }
#' }
#'
#' @return
#' Invisibly returns \code{TRUE} if \code{meta} is valid. Otherwise, raises an
#' error via \code{rlang::abort()} describing the violation.
#' 
.check_cluster_genes_multiomics_meta <- function(meta) {
    if (!is.data.frame(meta)) {
        rlang::abort(c(
            "Invalid `meta`.",
            "i" = "`meta` must be a data frame."
        ))
    }
    
    if (nrow(meta) < 1L) {
        rlang::abort(c(
            "Invalid `meta`.",
            "i" = "`meta` must have at least one row."
        ))
    }
    
    required_cols <- c("modality", "many_to_one_k", "modality_w")
    missing_cols <- setdiff(required_cols, colnames(meta))
    
    if (length(missing_cols) > 0L) {
        rlang::abort(c(
            "Invalid `meta`.",
            "i" = "Missing required columns.",
            "x" = paste0(
                "Missing: ",
                paste(missing_cols, collapse = ", ")
            )
        ))
    }
    
    modality <- meta$modality
    
    if (!is.character(modality)) {
        rlang::abort(c(
            "Invalid `meta$modality`.",
            "i" = "`modality` must be a character vector."
        ))
    }
    
    if (anyNA(modality) || any(modality == "")) {
        rlang::abort(c(
            "Invalid `meta$modality`.",
            "i" = "Modality names must be non-empty."
        ))
    }
    
    if (anyDuplicated(modality)) {
        rlang::abort(c(
            "Invalid `meta$modality`.",
            "i" = "Each modality must appear only once.",
            "x" = paste0(
                "Duplicated: ",
                paste(
                    unique(modality[duplicated(modality)]),
                    collapse = ", "
                )
            )
        ))
    }
    
    many_to_one_k <- meta$many_to_one_k
    
    if (!is.numeric(many_to_one_k)) {
        rlang::abort(c(
            "Invalid `meta$many_to_one_k`.",
            "i" = "`many_to_one_k` must be numeric (NA for one-to-one)."
        ))
    }
    
    bad_k <- !is.na(many_to_one_k) &
        (many_to_one_k < 2 |
             many_to_one_k != floor(many_to_one_k))
    
    if (any(bad_k)) {
        rlang::abort(c(
            "Invalid `meta$many_to_one_k`.",
            "i" = "Non-NA values must be integers >= 2.",
            "x" = paste0(
                "Invalid for: ",
                paste(modality[bad_k], collapse = ", ")
            )
        ))
    }
    
    modality_w <- meta$modality_w
    
    if (!is.numeric(modality_w) ||
        anyNA(modality_w) ||
        any(!is.finite(modality_w))) {
        rlang::abort(c(
            "Invalid `meta$modality_w`.",
            "i" = "`modality_w` must be finite numeric values."
        ))
    }
    
    if (any(modality_w <= 0)) {
        rlang::abort(c(
            "Invalid `meta$modality_w`.",
            "i" = "`modality_w` must be strictly positive.",
            "x" = paste0(
                "Non-positive for: ",
                paste(modality[modality_w <= 0], collapse = ", ")
            )
        ))
    }
    
    invisible(TRUE)
}


#' Validate consistency between data and modality metadata
#'
#' @noRd
#'
#' @description
#' Internal helper that checks cross-consistency between the modalities
#' present in the input \code{data} object and the rows of the modality
#' metadata table \code{meta}.
#'
#' The function ensures that:
#'
#' \itemize{
#'   \item{Every modality appearing in any condition of \code{data} has a
#'   corresponding entry in \code{meta}.}
#'
#'   \item{Every row in \code{meta} corresponds to a modality that is actually
#'   present in \code{data}.}
#' }
#'
#' This check guarantees that modality-level parameters (weights,
#' many-to-one settings) are well-defined for all data blocks and that no
#' extraneous metadata entries are provided.
#'
#' @param data
#' Named list of conditions. Each condition is a named list of modality
#' matrices.
#'
#' @param meta
#' Data frame of modality metadata. Must include a \code{modality} column
#' naming each modality.
#'
#' @return
#' Invisibly returns \code{TRUE} if \code{data} and \code{meta} are consistent.
#' Otherwise, raises an error via \code{rlang::abort()}.
#' 
.check_cluster_genes_multiomics_cross_args <- function(
        data,
        meta
        ) {
    data_mods <- sort(unique(unlist(lapply(data, names), use.names = FALSE)))
    meta_mods <- meta$modality
    
    missing_in_meta <- setdiff(data_mods, meta_mods)
    if (length(missing_in_meta) > 0L) {
        rlang::abort(c(
            "Invalid inputs.",
            "i" = "All modalities in `data` must have an entry in `meta`.",
            "x" = paste0(
                "Missing in `meta`: ",
                paste(missing_in_meta, collapse = ", ")
            )
        ))
    }
    
    extra_in_meta <- setdiff(meta_mods, data_mods)
    if (length(extra_in_meta) > 0L) {
        rlang::abort(c(
            "Invalid inputs.",
            "i" = "All rows in `meta` must correspond to a modality in `data`.",
            "x" = paste0(
                "Unknown in `meta`: ",
                paste(extra_in_meta, collapse = ", ")
            )
        ))
    }
    
    invisible(TRUE)
}


## .build_gene_centric_representation ------------------------------------------


#' Align modality metadata to the modalities present in the data
#'
#' @noRd
#'
#' @description
#' Internal helper that reorders and augments the modality-level metadata so
#' that it is aligned with the set and order of modalities actually present
#' in the input data.
#'
#' The function subsets \code{meta} to the modalities specified by
#' \code{modality_ids}, preserving their order. It then derives additional
#' modality-level annotations required downstream:
#'
#' \itemize{
#'   \item{A logical flag indicating whether each modality is many-to-one.}
#'   \item{Normalized modality weights that sum to one across modalities.}
#' }
#'
#' The aligned metadata is used consistently for normalization, weighting,
#' feature concatenation, and centroid/QC computation.
#'
#' @param meta
#' Data frame with one row per modality describing modality-level parameters.
#' Must contain at least the columns \code{modality}, \code{many_to_one_k},
#' and \code{modality_w}.
#'
#' @param modality_ids
#' Character vector of modality identifiers defining the desired order and
#' subset of modalities.
#'
#' @return
#' A data frame containing the aligned modality metadata, augmented with:
#'
#' \describe{
#'   \item{\code{is_many_to_one}}{
#'     Logical flag indicating whether the modality is many-to-one.}
#'
#'   \item{\code{modality_w_norm}}{
#'     Numeric vector of modality weights normalized to sum to one.}
#' }
#' 
.align_modality_meta <- function(
        meta,
        modality_ids
) {
    idx <- match(modality_ids, meta$modality)
    meta <- meta[idx, , drop = FALSE]
    meta$is_many_to_one <- !is.na(meta$many_to_one_k)
    w <- meta$modality_w
    meta$modality_w_norm <- w / sum(w)
    rownames(meta) <- NULL
    meta
}


#' Collapse many-to-one modalities to gene-level representations
#'
#' @noRd
#'
#' @description
#' Internal helper that converts all many-to-one modalities in the input data
#' into gene-level matrices by building pattern signature representations.
#'
#' The function processes modalities using the modality metadata in
#' \code{meta}:
#'
#' \itemize{
#'   \item{One-to-one modalities are passed through unchanged for each
#'   condition.}
#'
#'   \item{Many-to-one modalities (e.g. site- or probe-level measurements) are
#'   pooled \emph{across all conditions} for a given modality and collapsed
#'   using \code{.build_site_signatures()}. Feature trajectories are clustered
#'   into \code{many_to_one_k} global temporal archetypes, and each gene is
#'   represented by a fractional signature vector over these shared
#'   archetypes.}
#' }
#'
#' For many-to-one modalities, archetypes are fit once per modality on the
#' pooled feature matrix. Condition-specific gene signatures are obtained by
#' temporarily treating each gene-condition pair as a distinct gene identifier
#' during signature construction and then splitting the resulting signature
#' matrix back into per-condition gene-level matrices.
#'
#' Quality-control summaries from the many-to-one clustering step are
#' collected across many-to-one modalities and returned alongside the gene-
#' level matrices. QC rows corresponding to pooled fits are labeled with
#' \code{condition = "POOLED"}.
#'
#' @param data
#' Named list of conditions. Each condition is a named list of modality
#' matrices. One-to-one modalities must use gene IDs as row names. Many-to-one
#' modalities must use row names that encode gene identity (e.g.
#' \code{<gene>_<feature>}).
#'
#' @param meta
#' Data frame of modality metadata aligned to the modalities present in
#' \code{data}. Must include \code{modality} and \code{many_to_one_k}. A value
#' of \code{NA} in \code{many_to_one_k} indicates a one-to-one modality.
#'
#' @param verbose
#' Logical scalar indicating whether informative messages should be emitted.
#' Currently unused in this helper, but included for a consistent internal
#' interface.
#'
#' @return
#' A list with two components:
#'
#' \describe{
#'   \item{\code{mats_gene}}{
#'     Nested list of gene-level matrices structured as
#'     \code{mats_gene[[condition]][[modality]]}. One-to-one modalities are
#'     unchanged; many-to-one modalities are replaced by gene-level signature
#'     matrices whose columns correspond to shared global archetypes.}
#'
#'   \item{\code{many_to_one_clustering_qc}}{
#'     Data frame of QC summaries for pooled many-to-one feature clustering
#'     steps (one per many-to-one modality), or \code{NULL} if no many-to-one
#'     modalities are present.}
#' }
#' 
.collapse_many_to_one_modalities <- function(
        data,
        meta,
        verbose
) {
    cond_ids <- names(data)
    mod_ids  <- meta$modality

    mats_gene <- setNames(vector(
        "list",
        length(cond_ids)
        ), cond_ids)
    qc_list <- list()
    
    for (cond in cond_ids) {
        mats_gene[[cond]] <- setNames(vector(
            "list",
            length(mod_ids)
            ), mod_ids)
    }
    
    for (i in seq_along(mod_ids)) {
        m  <- mod_ids[[i]]
        lk <- meta$many_to_one_k[[i]]
        
        if (is.na(lk)) { 
            # means m is a one-to-one modality and is skipped
            for (cond in cond_ids) {
                mats_gene[[cond]][[m]] <- data[[cond]][[m]]
            }
            next
        }
        
        mats_m <- lapply(cond_ids, function(cond) data[[cond]][[m]])
        names(mats_m) <- cond_ids
        mats_m <- mats_m[!vapply(mats_m, is.null, logical(1))]
        
        pooled_mat <- do.call(rbind, mats_m)
        cond_of_row <- rep(
            names(mats_m),
            vapply(
                mats_m,
                nrow,
                integer(1)
                )
            )
        
        rn <- rownames(pooled_mat)
        gene_raw <- sub("_.*$", "", rn)
        
        # “Lift” the gene IDs to be condition-specific  
        feature_to_gene <- paste0(
            gene_raw,
            "__",
            cond_of_row
            )
        
        # Fit archetypes globally and build signatures
        sig_obj <- .build_site_signatures(
            modality_mat    = pooled_mat,
            feature_to_gene = feature_to_gene,
            many_to_one_k   = as.integer(lk)
        )
        
        sig_all <- sig_obj$signatures
        
        # Split signatures back into per-condition gene matrices
        for (cond in cond_ids) { 
            suffix <- paste0("__", cond, "$")
            keep <- grepl(suffix, rownames(sig_all))
            
            sig_c <- sig_all[keep, , drop = FALSE]
            rownames(sig_c) <- sub(suffix, "", rownames(sig_c))
            
            mats_gene[[cond]][[m]] <- sig_c[order(rownames(sig_c)), ,
                                            drop = FALSE]
        }
        
        # Collect QC from the pooled clustering run
        qc <- sig_obj$many_to_one_clustering_qc
        if (is.data.frame(qc) && nrow(qc) > 0L) {
            qc$condition <- "POOLED"
            qc$modality <- m
            qc$many_to_one_k <- as.integer(lk)
            qc_list[[length(qc_list) + 1L]] <- qc
        }
    }
    
    many_to_one_clustering_qc <- if (length(qc_list) > 0L) {
        do.call(rbind, qc_list)
    } else {
        NULL
    }
    
    list(
        mats_gene = mats_gene,
        many_to_one_clustering_qc = many_to_one_clustering_qc
    )
}


#' Normalize gene-level modality matrices across conditions
#'
#' @noRd
#'
#' @description
#' Internal helper that applies modality-appropriate normalization to all
#' gene-level matrices prior to gene harmonization and joint feature
#' construction.
#'
#' Normalization is applied independently for each condition × modality block
#' in \code{mats_gene}. The normalization strategy depends on the modality
#' type, as defined in \code{meta}:
#'
#' \itemize{
#'   \item{One-to-one modalities are normalized by row-wise z-scoring, making
#'   clustering sensitive to relative patterns rather than absolute scale.}
#'
#'   \item{Many-to-one modalities (gene-level pattern signatures) are
#'   normalized using a Hellinger transform (square root of non-negative
#'   entries), yielding a geometry appropriate for Euclidean methods on
#'   compositional data.}
#' }
#'
#' The output preserves the original nested structure of the input, enabling
#' downstream alignment, weighting, and centroid/QC computation.
#'
#' @param mats_gene
#' Nested list of gene-level matrices structured as
#' \code{mats_gene[[condition]][[modality]]}. Row names are gene identifiers.
#'
#' @param meta
#' Data frame of modality metadata aligned to the modalities present in
#' \code{mats_gene}. Must include \code{modality} and \code{many_to_one_k}.
#'
#' @param verbose
#' Logical scalar indicating whether informative messages should be emitted.
#' Currently unused in this helper, but included for a consistent internal
#' interface.
#'
#' @return
#' Nested list of normalized gene-level matrices with the same structure as
#' \code{mats_gene}: \code{mats_norm[[condition]][[modality]]}.
#' 
.normalize_gene_level_modalities <- function(
        mats_gene,
        meta,
        verbose
) {
    cond_ids <- names(mats_gene)
    mod_ids <- meta$modality
    
    mats_norm <- vector("list", length(cond_ids))
    names(mats_norm) <- cond_ids
    
    is_m2o <- !is.na(meta$many_to_one_k)
    names(is_m2o) <- mod_ids
    
    for (cond in cond_ids) {
        mats_cond <- vector("list", length(mod_ids))
        names(mats_cond) <- mod_ids
        
        for (m in mod_ids) {
            mats_cond[[m]] <- .normalize_modality_mat(
                mat = mats_gene[[cond]][[m]],
                is_many_to_one = is_m2o[[m]]
            )
        }
        
        mats_norm[[cond]] <- mats_cond
    }
    
    mats_norm
}


#' Harmonize the gene universe across condition × modality blocks
#'
#' @noRd
#'
#' @description
#' Internal helper that determines the set of genes to be included in the
#' gene-centric analysis by harmonizing row identities across all normalized
#' condition × modality blocks in \code{mats_norm}.
#'
#' The gene universe is derived from the row names of every block
#' \code{mats_norm[[condition]][[modality]]}. The harmonization policy is
#' controlled by \code{gene_mode}:
#'
#' \itemize{
#'   \item{\code{"intersection"}}{Retain only genes present in every block.}
#'   \item{\code{"union"}}{Retain genes present in at least one block. Missing
#'   genes in a block are handled downstream during concatenation by filling
#'   zeros.}
#' }
#'
#' The returned vector defines the global gene ordering used for building the
#' joint feature matrix and ensures consistent alignment across downstream
#' objects (PCA scores, kNN graph, adjacency matrix, clustering labels).
#'
#' @param mats_norm
#' Nested list of normalized gene-level matrices structured as
#' \code{mats_norm[[condition]][[modality]]}. Row names are gene identifiers.
#'
#' @param gene_mode
#' Character string specifying the harmonization policy. Must be one of
#' \code{"intersection"} or \code{"union"}.
#'
#' @return
#' Character vector of gene identifiers defining the gene universe and a
#' deterministic ordering for downstream alignment.
#' 
.harmonize_genes_across_blocks <- function(
        mats_norm,
        gene_mode
        ) {
    gene_sets <- unlist(
        lapply(mats_norm, function(mats_cond) {
            lapply(mats_cond, rownames)
        }),
        recursive = FALSE,
        use.names = FALSE
    )
    
    if (gene_mode == "intersection") {
        return(Reduce(intersect, gene_sets))
    }
    
    sort(unique(unlist(gene_sets, use.names = FALSE)))
}


#' Flatten normalized condition × modality blocks into a joint feature matrix
#'
#' @noRd
#'
#' @description
#' Internal helper that constructs the **gene-centric joint feature matrix**
#' used for PCA by aligning and concatenating all normalized gene-level blocks
#' across conditions and modalities.
#'
#' For each condition × modality block in \code{mats_norm}, the function:
#'
#' \itemize{
#'   \item{Aligns the block to the global gene universe \code{gene_ids}. Genes
#'   missing from a block are filled with zeros.}
#'   \item{Applies the modality weight (from \code{meta$modality_w_norm}) to all
#'   columns of that block. Weighting is applied as \code{sqrt(w)} so that
#'   squared Euclidean distances in the concatenated space reflect the intended
#'   modality contribution.}
#'   \item{Prefixes block column names with condition and modality identifiers
#'   (\code{<condition>__<modality>__<feature>}) to preserve provenance.}
#'   \item{Records a feature-level index mapping each output column back to its
#'   originating condition, modality, feature name, and modality weight.}
#' }
#'
#' The resulting matrix has one row per gene and columns spanning all
#' condition-specific representations across all modalities.
#'
#' @param mats_norm
#' Nested list of normalized gene-level matrices structured as
#' \code{mats_norm[[condition]][[modality]]}. Row names are gene identifiers.
#'
#' @param meta
#' Data frame of modality metadata aligned to the modalities present in
#' \code{mats_norm}. Must contain \code{modality} and \code{modality_w_norm}.
#'
#' @param gene_ids
#' Character vector giving the global gene universe and row order for the
#' output matrix. All blocks are aligned to this ordering.
#'
#' @param verbose
#' Logical scalar indicating whether informative messages should be emitted.
#' Currently unused in this helper, but included for a consistent internal
#' interface.
#'
#' @return
#' A list with two components:
#'
#' \describe{
#'   \item{\code{X}}{
#'     Numeric matrix (\code{genes x features}) giving the joint feature matrix
#'     formed by concatenating all weighted condition × modality blocks. Row
#'     names are \code{gene_ids}. Column names encode provenance.}
#'
#'   \item{\code{feature_index}}{
#'     Data frame mapping each column of \code{X} to its originating condition,
#'     modality, within-block feature name, and the modality weight used. The
#'     \code{column} field matches \code{colnames(X)}.}
#' }
#' 
.flatten_blocks_to_joint_matrix <- function(
        mats_norm,
        meta,
        gene_ids,
        verbose
) {
    cond_ids <- names(mats_norm)
    mod_ids <- meta$modality
    
    w <- meta$modality_w_norm
    names(w) <- mod_ids
    
    blocks <- list()
    feat_index <- list()
    
    for (cond in cond_ids) {
        for (m in mod_ids) {
            mat <- mats_norm[[cond]][[m]]
            
            mat_aligned <- matrix(
                0,
                nrow = length(gene_ids),
                ncol = ncol(mat),
                dimnames = list(gene_ids, colnames(mat))
            )
            
            idx <- match(rownames(mat), gene_ids)
            keep <- !is.na(idx)
            
            if (any(keep)) {
                mat_aligned[idx[keep], ] <- mat[keep, , drop = FALSE]
            }
            
            # Block-wise scaling: normalize average squared L2 norm per gene
            # to 1.
            # This stabilizes heterogeneous blocks 
            # (different feature counts / geometries)
            # before applying modality weights and concatenation.
            row_norm2 <- rowSums(mat_aligned^2)
            present <- row_norm2 > 0
            scale_fac <- sqrt(mean(row_norm2[present]))
            if (is.finite(scale_fac) && scale_fac > 0) {
                mat_aligned <- mat_aligned / scale_fac
            }
            
            block_w <- sqrt(w[[m]])
            mat_aligned <- mat_aligned * block_w
            
            cn <- colnames(mat_aligned)
            if (is.null(cn)) {
                cn <- paste0("V", seq_len(ncol(mat_aligned)))
            }
            
            colnames(mat_aligned) <- paste0(cond, "__", m, "__", cn)
            
            blocks[[length(blocks) + 1L]] <- mat_aligned
            
            feat_index[[length(feat_index) + 1L]] <- data.frame(
                condition = cond,
                modality = m,
                feature = cn,
                weight = w[[m]],
                stringsAsFactors = FALSE
            )
        }
    }
    
    X <- do.call(cbind, blocks)
    
    feature_index <- do.call(rbind, feat_index)
    feature_index$column <- colnames(X)
    
    list(
        X = X,
        feature_index = feature_index
    )
}


# Level 3 function definitions -------------------------------------------------


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
        rlang::abort(c(
            "`feature_to_gene` must have length equal to nrow(modality_mat). ",
            "Mapping is positional: element i is the gene ID for row i."
        ))
    }
    
    X_scaled <- scale(modality_mat)

    feature_clusters <- .cluster_feature_matrix(
        feature_mat = X_scaled,
        k           = many_to_one_k
    )

    many_to_one_clustering_qc <- .compute_cluster_centroids_qc(
        X  = X_scaled,
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


# Level 4 function definitions -------------------------------------------------


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