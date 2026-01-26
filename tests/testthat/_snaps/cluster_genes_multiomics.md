# cluster_genes_multiomics works.

    structure(list(gene = c("gene1", "gene2", "gene3", "gene4", "gene5", 
    "gene6"), time_Ctrl = c(1L, 2L, 2L, 2L, 2L, 1L), time_Treat = c(1L, 
    2L, 1L, 1L, 2L, 1L)), class = c("tbl_df", "tbl", "data.frame"
    ), row.names = c(NA, -6L))

---

    structure(list(block = c("time_Ctrl", "time_Ctrl", "time_Ctrl", 
    "time_Ctrl", "time_Treat", "time_Treat", "time_Treat", "time_Treat"
    ), modality = c("rna", "rna", "site", "site", "rna", "rna", "site", 
    "site"), modality_type = c("one_to_one", "one_to_one", "many_to_one", 
    "many_to_one", "one_to_one", "one_to_one", "many_to_one", "many_to_one"
    ), cluster = c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L), n_genes_cluster = c(2L, 
    4L, 2L, 4L, 4L, 2L, 4L, 2L), n_genes_used = c(2L, 4L, 2L, 4L, 
    4L, 2L, 4L, 2L), coverage = c(1, 1, 1, 1, 1, 1, 1, 1), mean_R2 = c(0.10631112086932, 
    0.358883594848579, 0.933012701892219, 0.75, 0.482504199033334, 
    0.353467788837418, 0.95064540160315, NaN), sd_R2 = c(1.27946881663023e-16, 
    0.224450218796353, 0, 0, 0.336116177527721, 1.24126707662364e-16, 
    0.0658061311957994, NA), r2_member = structure(list(c(gene1 = 0.10631112086932, 
    gene6 = 0.10631112086932), c(gene2 = 0.548935732855668, gene3 = 0.0817965723356662, 
    gene4 = 0.533643433808088, gene5 = 0.271158640394894), c(gene1 = 0.933012701892219, 
    gene6 = 0.933012701892219), c(gene2 = 0.75, gene3 = 0.75, gene4 = NA, 
    gene5 = NA), c(gene1 = 0.644497800584004, gene3 = 0.0222647931628325, 
    gene4 = 0.463721269205557, gene6 = 0.799532933180941), c(gene2 = 0.353467788837418, 
    gene5 = 0.353467788837418), c(gene1 = 0.851936204809451, gene3 = 0.98354846720105, 
    gene4 = 0.98354846720105, gene6 = 0.98354846720105), c(gene2 = NA_real_, 
    gene5 = NA_real_)), class = "AsIs"), centroid = structure(list(
        c(0.0183294109362451, 0.38833706059501, -0.437243744247436, 
        0.218332063798552, -0.187754791082371), c(0.494475811037839, 
        -0.228080339759566, 0.716533658758647, -0.538683491305696, 
        -0.444245638731224), c(cluster_1 = -0.288675134594813, cluster_2 = -0.788675134594813, 
        cluster_3 = 1.07735026918963), c(cluster_1 = 0.5, cluster_2 = 0.5, 
        cluster_3 = -1), c(0.701685936927369, -0.237525633850794, 
        -0.763197517251794, 0.60813989529384, -0.309102681118621), 
        c(0.323977795361364, -0.886741841721731, 0.364441032721987, 
        -0.331021695085326, 0.529344708723705), c(cluster_1 = 1.03867513459481, 
        cluster_2 = -0.144337567297406, cluster_3 = -0.894337567297406
        ), c(cluster_1 = NaN, cluster_2 = NaN, cluster_3 = NaN)), class = "AsIs")), class = c("tbl_df", 
    "tbl", "data.frame"), row.names = c(NA, -8L))

