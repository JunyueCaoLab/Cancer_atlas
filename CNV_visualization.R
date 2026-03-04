library(dplyr)

visualize_inferCNV_heatmap <- function(inferCNV_out_df, chr_len_df, cell_line_to_visualize_in_order){

    # deps: dplyr only if you keep dplyr verbs; below is base-safe.

    chrom_to_keep <- paste0("chr", 1:22)

    chr_len_df <- chr_len_df[chr_len_df$V1 %in% chrom_to_keep, , drop = FALSE]
    chr_len_df <- chr_len_df[match(chrom_to_keep, chr_len_df$V1), , drop = FALSE]  # enforce order

    # Parse cell line name
    inferCNV_out_df$cell_line_names <- vapply(
        as.character(inferCNV_out_df$cell_group_name),
        function(x) strsplit(x, ".", fixed = TRUE)[[1]][1],
        character(1)
    )

    inferCNV_out_df <- inferCNV_out_df[inferCNV_out_df$cell_line_names %in% cell_line_to_visualize_in_order, , drop = FALSE]

    # Map inferCNV state -> "actual_status"
    state_map <- c(`1`=0, `2`=0.5, `3`=1, `4`=1.5, `5`=2, `6`=3)
    inferCNV_out_df$actual_status <- unname(state_map[as.character(inferCNV_out_df$state)])

    # Build chromosome skeleton with cumulative coords
    # start_offset: 0-based genome coordinate where each chromosome starts
    chr_len <- as.numeric(chr_len_df$V2)
    start_offset <- c(0, cumsum(chr_len)[-length(chr_len)])  # chr1 starts 0, chr2 starts len(chr1), ...
    end_coord <- start_offset + chr_len                       # end (exclusive) if you want half-open intervals

    chr_skeleton_df <- data.frame(
        V1 = chrom_to_keep,
        V2 = end_coord,               # keep your original column naming, but now it's cumulative end
        start_offset = start_offset,
        chr_len = chr_len,
        stringsAsFactors = FALSE
    )

    # Vectorized coordinate shift using chromosome start_offset
    idx <- match(inferCNV_out_df$chr, chr_skeleton_df$V1)
    inferCNV_out_df$start <- inferCNV_out_df$start + chr_skeleton_df$start_offset[idx]
    inferCNV_out_df$end   <- inferCNV_out_df$end   + chr_skeleton_df$start_offset[idx]

    # Factor ordering for plotting (reverse so first is on top)
    inferCNV_out_df$cell_lines_fct <- factor(
        inferCNV_out_df$cell_line_names,
        levels = rev(cell_line_to_visualize_in_order)
    )

    return(list(
        CNV_block_df = inferCNV_out_df,
        chr_skeleton_df = chr_skeleton_df
    ))
}