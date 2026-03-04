library(parallel)
library(dplyr)

check_gene_wise_drug_eff_cnv_association <- function(
    genes_of_interest,
    cell_line_gene_cnv_df,
    cell_line_drug_trt_df,
    mc.cores = 16,
    drug_cell_line_start_col = 5,
    require_three_categories = TRUE,
    drop_nonpositive_ic50 = TRUE
){

    genes_of_interest_trimmed <- intersect(genes_of_interest, rownames(cell_line_gene_cnv_df))

    ### Precompute drug response matrix portion and its colnames
    drug_cell_lines <- colnames(cell_line_drug_trt_df)[drug_cell_line_start_col:ncol(cell_line_drug_trt_df)]

    gene_wise <- lapply(genes_of_interest_trimmed, function(x){

        message(x)

        cnv_vec <- as.numeric(cell_line_gene_cnv_df[x, drug_cell_lines, drop = TRUE]) # align to drug cols if possible
        ### If CNV df has different cols than drug df, fall back to all cnv cols
        if(all(is.na(cnv_vec))) {
            cnv_vec <- as.numeric(cell_line_gene_cnv_df[x, , drop = TRUE])
            cnv_cell_lines <- colnames(cell_line_gene_cnv_df)
        } else {
            cnv_cell_lines <- drug_cell_lines
        }

        amp_cell_lines     <- cnv_cell_lines[cnv_vec > 0]
        diploid_cell_lines <- cnv_cell_lines[cnv_vec == 0]
        del_cell_lines     <- cnv_cell_lines[cnv_vec < 0]

        cell_gene_cnv_df <- data.frame(
            cell_lines = c(amp_cell_lines, diploid_cell_lines, del_cell_lines),
            category   = c(rep("amp", length(amp_cell_lines)),
                           rep("diploid", length(diploid_cell_lines)),
                           rep("del", length(del_cell_lines))),
            stringsAsFactors = FALSE
        )

        ### per-drug analysis
        drug_wise_examination <- parallel::mclapply(seq_len(nrow(cell_line_drug_trt_df)), mc.cores = mc.cores, function(y){

            ### IC50 values for this drug across cell lines
            ic50_raw <- suppressWarnings(as.numeric(cell_line_drug_trt_df[y, drug_cell_lines, drop = TRUE]))
            names(ic50_raw) <- drug_cell_lines

            ok <- !is.na(ic50_raw)
            if(drop_nonpositive_ic50) ok <- ok & (ic50_raw > 0)

            intersected_cell_lines <- intersect(names(ic50_raw)[ok], cell_gene_cnv_df$cell_lines)

            base_row <- data.frame(
                gene = x,
                drug = as.character(cell_line_drug_trt_df[y, 1]),
                feature = as.character(cell_line_drug_trt_df[y, 4]),
                anova_pval = NA_real_,
                intercept_estimate = NA_real_,
                cnv_term_estimate = NA_real_,
                intercept_pval = NA_real_,
                cnv_term_pval = NA_real_,
                del_ic50_median = NA_real_,
                diploid_ic50_median = NA_real_,
                amp_ic50_median = NA_real_,
                stringsAsFactors = FALSE
            )

            if(length(intersected_cell_lines) == 0) return(base_row)

            drug_df <- cell_gene_cnv_df[cell_gene_cnv_df$cell_lines %in% intersected_cell_lines, , drop = FALSE]

            ### category requirement
            if(require_three_categories && length(unique(drug_df$category)) != 3) return(base_row)

            #### attach response
            drug_df$log10_ic50 <- log10(ic50_raw[drug_df$cell_lines])

            #### ANOVA
            anova_pval <- tryCatch(
                summary(aov(log10_ic50 ~ category, data = drug_df))[[1]][["Pr(>F)"]][1],
                error = function(e) NA_real_
            )

            #### ordinal CNV coding
            drug_df$copy_number <- 0
            drug_df$copy_number[drug_df$category == "amp"] <-  2
            drug_df$copy_number[drug_df$category == "del"] <- -2

            lm_model <- tryCatch(lm(log10_ic50 ~ copy_number, data = drug_df), error = function(e) NULL)
            if(is.null(lm_model)) {
                base_row$anova_pval <- anova_pval
                return(base_row)
            }

            coef_df <- summary(lm_model)$coefficients

            ### medians
            amp_med <- median(drug_df$log10_ic50[drug_df$category == "amp"], na.rm = TRUE)
            dip_med <- median(drug_df$log10_ic50[drug_df$category == "diploid"], na.rm = TRUE)
            del_med <- median(drug_df$log10_ic50[drug_df$category == "del"], na.rm = TRUE)

            base_row$anova_pval <- anova_pval
            if(nrow(coef_df) == 2){
                base_row$intercept_estimate <- as.numeric(coef_df[1, 1])
                base_row$cnv_term_estimate  <- as.numeric(coef_df[2, 1])
                base_row$intercept_pval     <- as.numeric(coef_df[1, 4])
                base_row$cnv_term_pval      <- as.numeric(coef_df[2, 4])
                base_row$del_ic50_median    <- del_med
                base_row$diploid_ic50_median<- dip_med
                base_row$amp_ic50_median    <- amp_med
            }
            base_row
        })

        do.call(rbind, drug_wise_examination)
    })

    do.call(rbind, gene_wise)
}