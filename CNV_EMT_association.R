library("parallel")
library("dplyr")
library("rtracklayer")
library("GenomeInfoDb")

manhattan_plot_prep.v2 <- function(
    gene_input_df,
    gene_symbol_col,
    chromosome_length_df,
    gtf_path){
    
    ### import the gtf
    annotation_gtf <- rtracklayer::import(gtf_path, format="gtf")
    annotation_gtf <- annotation_gtf[annotation_gtf$type == "gene"]
    annotation_gtf <- GenomeInfoDb::keepStandardChromosomes(annotation_gtf, pruning.mode="coarse")
    annotation_gtf <- annotation_gtf[!duplicated(annotation_gtf$gene_name)]

    ### assign chromosome and coordinates
    intersected_genes <- intersect(annotation_gtf$gene_name, as.character(gene_input_df[,gene_symbol_col]))
    trimmed_gene_df <- gene_input_df[
        as.character(gene_input_df[,gene_symbol_col]) %in% intersected_genes,
    ]

    trimmed_gene_df$chromosome <- as.character(seqnames(annotation_gtf))[
        match(as.character(trimmed_gene_df[,gene_symbol_col]), annotation_gtf$gene_name)
    ]

    trimmed_gene_df$start <- start(annotation_gtf)[
        match(as.character(trimmed_gene_df[,gene_symbol_col]), annotation_gtf$gene_name)
    ]

    trimmed_gene_df$end <- end(annotation_gtf)[
        match(as.character(trimmed_gene_df[,gene_symbol_col]), annotation_gtf$gene_name)
    ]

    trimmed_gene_df$mean_coordinate <- rowMeans(
        cbind(trimmed_gene_df$start, trimmed_gene_df$end)
    )

    ### order each chromosome
    chromosomes_to_keep <- paste("chr", 1:22, sep="")

    chromosome_ordered_df <- lapply(chromosomes_to_keep, function(chr){
        subdf <- dplyr::filter(trimmed_gene_df, chromosome == chr)
        subdf[order(subdf$mean_coordinate), ]
    })

    reordered_chromosome_FDR_df <- do.call(rbind, chromosome_ordered_df)
    reordered_chromosome_FDR_df$order <- seq_len(nrow(reordered_chromosome_FDR_df))

    ### build genomic coordinates across chromosomes
    total_len <- 0
    places <- c()
    new_df_list <- list()

    for(chr in chromosomes_to_keep){

        subdf <- dplyr::filter(trimmed_gene_df, chromosome == chr)

        coordinates <- total_len + subdf$mean_coordinate
        subdf$coordinates <- coordinates

        new_df_list[[chr]] <- subdf

        total_len <- total_len + chromosome_length_df[[chr]]
        places <- c(places, total_len)
    }

    new_df <- do.call(rbind, new_df_list)

    return(list(
        final_reordered_chromosome_FDR_df = new_df,
        places = places
    ))
}

infer_CNV_ATAC_EMT_manhattan <- function(
    median_EMT_df,
    gene_copy_number_mat,
    gene_activity_mat = NULL,
    chromosome_length_df,
    mc.cores = 16,
    gtf_path = "/index/GTF/human-latest-gencode-release-38/gencode.v38.primary_assembly.annotation.gtf.gz",
    gene_symbol_col = "amplified_gene",
    genes_to_examine = NULL,
    manhattan_pval_col = "cnv_term_pval",
    add_manhattan_y = TRUE
){

    ### CNV ~ EMT
    gene_wise_association <- parallel::mclapply(
        rownames(gene_copy_number_mat),
        mc.cores = mc.cores,
        function(x){
            EMT_pseudotime_CNV_test_df <- median_EMT_df
            EMT_pseudotime_CNV_test_df$cnv_val <- as.numeric(
                gene_copy_number_mat[x, EMT_pseudotime_CNV_test_df$cell_lines]
            )

            lm_model <- lm(data = EMT_pseudotime_CNV_test_df, formula = median_pseudotime ~ cnv_val)
            coef_df <- summary(lm_model)$coefficients

            if(nrow(coef_df) == 2){
                data.frame(
                    amplified_gene     = x,
                    intercept_estimate = as.numeric(coef_df[1, 1]),
                    cnv_term_estimate  = as.numeric(coef_df[2, 1]),
                    intercept_pval     = as.numeric(coef_df[1, 4]),
                    cnv_term_pval      = as.numeric(coef_df[2, 4]),
                    stringsAsFactors = FALSE
                )
            } else {
                NULL
            }
        }
    )
    cnv_df <- do.call(rbind, gene_wise_association)
    if(is.null(cnv_df) || nrow(cnv_df) == 0) stop("CNV association returned empty results.")

    ### decide genes for ATAC validation
    if(is.null(genes_to_examine)) {
        genes_for_atac <- cnv_df[[gene_symbol_col]]
    } else {
        genes_for_atac <- genes_to_examine
    }

    ### ATAC gene activity ~ EMT
    atac_df <- NULL
    if(!is.null(gene_activity_mat)) {
        genes_for_atac <- intersect(genes_for_atac, rownames(gene_activity_mat))

        atac_list <- lapply(genes_for_atac, function(x){
            EMT_df <- median_EMT_df
            EMT_df$gene_activity_val <- as.numeric(
                gene_activity_mat[x, EMT_df$cell_lines]
            )

            # cor() can be NA if constant / missing; guard lightly
            cor_coef <- suppressWarnings(cor(EMT_df$median_pseudotime, EMT_df$gene_activity_val, use = "pairwise.complete.obs"))
            cor_pval <- suppressWarnings(cor.test(EMT_df$median_pseudotime, EMT_df$gene_activity_val)$p.value)

            data.frame(
                amplified_gene = x,
                cor_coef = cor_coef,
                cor_pval = cor_pval,
                stringsAsFactors = FALSE
            )
        })
        atac_df <- do.call(rbind, atac_list)
    }

    ### merge CNV + ATAC
    merged_df <- cnv_df
    if(!is.null(atac_df) && nrow(atac_df) > 0) {
        merged_df <- merge(cnv_df, atac_df, by = "amplified_gene", all.x = TRUE)
    }

    ### make the df for visualization
    manhattan_in <- merged_df

    # add -log10(p) column if requested
    if(add_manhattan_y) {
        if(!manhattan_pval_col %in% colnames(manhattan_in)) {
            stop(sprintf("manhattan_pval_col '%s' not found in merged_df.", manhattan_pval_col))
        }
        manhattan_in$manhattan_y <- -log10(pmax(as.numeric(manhattan_in[[manhattan_pval_col]]), .Machine$double.xmin))
    }

    manhattan_out <- manhattan_plot_prep.v2(
        gene_input_df = manhattan_in,
        gene_symbol_col = gene_symbol_col,
        chromosome_length_df = chromosome_length_df,
        gtf_path = gtf_path
    )

    return(list(
        cnv_df = cnv_df,
        atac_df = atac_df,
        merged_df = merged_df,
        manhattan_df = manhattan_out$final_reordered_chromosome_FDR_df,
        manhattan_places = manhattan_out$places
    ))
}