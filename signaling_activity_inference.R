library(Seurat)
library(decoupleR)
library(reshape2)

infer_cell_line_signaling_activity_progeny <- function(
    expr_mat,
    cell_meta,
    cell_line_col = "cell_lines",
    cell_lines_to_keep = NULL,
    organism = "human",
    top = 500,
    assay = "RNA",
    slot = "data",
    minsize = 5
){
    
    # Build Seurat object + normalize
    seu <- Seurat::CreateSeuratObject(counts = expr_mat, meta.data = cell_meta)
    seu <- Seurat::NormalizeData(seu)

    # Subset cell lines
    seu <- subset(seu, subset = get(cell_line_col) %in% cell_lines_to_keep)

    # Get PROGENy network
    net <- decoupleR::get_progeny(organism = organism, top = top)

    # Run MLM
    mat_use <- Seurat::GetAssayData(seu, assay = assay, slot = slot)
    acts <- decoupleR::run_mlm(
        mat = mat_use,
        net = net,
        .source = "source",
        .target = "target",
        .mor = "weight",
        minsize = minsize
    )

    # get mat
    pathway_activity_score_mat <- reshape2::dcast(
        acts,
        formula = condition ~ source,
        value.var = "score"
    )
    rownames(pathway_activity_score_mat) <- pathway_activity_score_mat[, 1]
    pathway_activity_score_mat <- pathway_activity_score_mat[, -1, drop = FALSE]

    # Add to Seurat meta.data (cell-level)
    common_cells <- intersect(colnames(seu), rownames(pathway_activity_score_mat))
    if(length(common_cells) == 0){
        stop("No overlap between Seurat cell names and PROGENy 'condition' names. Check matrix column names.")
    }
    seu@meta.data[common_cells, colnames(pathway_activity_score_mat)] <-
        pathway_activity_score_mat[common_cells, , drop = FALSE]

    # Also return a per-cell-line summary (median activity per pathway)
    if(cell_line_col %in% colnames(seu@meta.data)){
        cell_line_summary <- aggregate(
            seu@meta.data[, colnames(pathway_activity_score_mat), drop = FALSE],
            by = list(cell_line = seu@meta.data[[cell_line_col]]),
            FUN = median,
            na.rm = TRUE
        )
        rownames(cell_line_summary) <- cell_line_summary$cell_line
    } else {
        cell_line_summary <- NULL
    }

    return(list(
        seurat = seu,
        acts_long = acts,
        pathway_activity_score_mat = pathway_activity_score_mat,
        cell_line_summary_median = cell_line_summary,
        net = net
    ))
}