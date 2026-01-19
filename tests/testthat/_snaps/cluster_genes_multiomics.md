# cluster_genes_multiomics snapshot with mixed layers

    structure(list(gene = c("gene1", "gene2", "gene3", "gene4", "gene5", 
    "gene6"), cluster_Ctrl = c(1L, 1L, 2L, 1L, 1L, 2L), cluster_Treat = c(2L, 
    1L, 2L, 2L, 1L, 2L)), class = c("tbl_df", "tbl", "data.frame"
    ), row.names = c(NA, -6L))

---

    structure(list(block = c("time_Ctrl", "time_Ctrl", "time_Ctrl", 
    "time_Ctrl", "time_Treat", "time_Treat", "time_Treat", "time_Treat"
    ), modality = c("rna", "rna", "site", "site", "rna", "rna", "site", 
    "site"), modality_type = c("one_to_one", "one_to_one", "many_to_one", 
    "many_to_one", "one_to_one", "one_to_one", "many_to_one", "many_to_one"
    ), cluster = c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L), n_genes_cluster = c(4L, 
    2L, 4L, 2L, 2L, 4L, 2L, 4L), n_genes_used = c(4L, 2L, 4L, 2L, 
    2L, 4L, 2L, 4L), coverage = c(1, 1, 1, 1, 1, 1, 1, 1), mean_R2 = c(0.454509281530177, 
    0.896641566062548, 0.0669872981077807, 0.25, 0.353467788837418, 
    0.482504199033334, NaN, 0.95064540160315), sd_R2 = c(0.118153903203798, 
    0, 1.67110694432208e-16, 0, 1.24126707662364e-16, 0.336116177527721, 
    NA, 0.0658061311957994), r2_member = structure(list(c(gene1 = 0.378928459568216, 
    gene2 = 0.625655160832676, gene4 = 0.440721765611455, gene5 = 0.372731740108359
    ), c(gene3 = 0.896641566062548, gene6 = 0.896641566062548), c(gene1 = 0.0669872981077806, 
    gene2 = 0.0669872981077808, gene4 = NA, gene5 = NA), c(gene3 = 0.25, 
    gene6 = 0.25), c(gene2 = 0.353467788837418, gene5 = 0.353467788837418
    ), c(gene1 = 0.644497800584004, gene3 = 0.0222647931628325, gene4 = 0.463721269205557, 
    gene6 = 0.799532933180941), c(gene2 = NA_real_, gene5 = NA_real_
    ), c(gene1 = 0.851936204809451, gene3 = 0.98354846720105, gene4 = 0.98354846720105, 
    gene6 = 0.98354846720105)), class = "AsIs"), centroid = structure(list(
        c(-0.104613070502927, 0.0837762408465411, 0.975020186792221, 
        -0.0452338024357303, -0.908949554700105), c(1.21650717401778, 
        -0.235376100617204, -0.954216800314583, -0.768567313941379, 
        0.74165304085539), c(cluster_1 = 0.0773502691896257, cluster_2 = 0.211324865405187, 
        cluster_3 = -0.288675134594813), c(cluster_1 = 0, cluster_2 = -0.5, 
        cluster_3 = 0.5), c(0.323977795361364, -0.886741841721731, 
        0.364441032721987, -0.331021695085326, 0.529344708723705), 
        c(0.701685936927369, -0.237525633850794, -0.763197517251794, 
        0.60813989529384, -0.309102681118621), c(cluster_1 = NaN, 
        cluster_2 = NaN, cluster_3 = NaN), c(cluster_1 = 1.03867513459481, 
        cluster_2 = -0.894337567297406, cluster_3 = -0.144337567297406
        )), class = "AsIs")), class = c("tbl_df", "tbl", "data.frame"
    ), row.names = c(NA, -8L))

