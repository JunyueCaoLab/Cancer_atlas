library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(clusterProfiler)
library(ArchR)
library(SummarizedExperiment)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
library(chromVAR)
library(Matrix)
library(destiny)
library(BiocParallel)

run_chromVAR_fullset_motifs <- function(
    atac_meta,
    atac_count_mat,
    motif_collection,
    out_rdata_path = NULL,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    n_cores = 16
){

    from_char_to_GRanges <- function(input_peak_char){
        peak_chrom <- vapply(input_peak_char, function(x) strsplit(x, ":", fixed = TRUE)[[1]][1], character(1))
        peak_start <- vapply(input_peak_char, function(x) as.numeric(strsplit(strsplit(x, ":", fixed = TRUE)[[1]][2], "-", fixed = TRUE)[[1]][1]), numeric(1))
        peak_end   <- vapply(input_peak_char, function(x) as.numeric(strsplit(strsplit(x, ":", fixed = TRUE)[[1]][2], "-", fixed = TRUE)[[1]][2]), numeric(1))

        GenomicRanges::GRanges(
            seqnames = peak_chrom,
            ranges = IRanges::IRanges(start = peak_start, end = peak_end),
            names = input_peak_char
        )
    }

    peak_to_motif_cell_mat_fullset <- function(ATAC_meta, ATAC_count_mat){
        rowranges <- from_char_to_GRanges(rownames(ATAC_count_mat))

        atac_se <- SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = ATAC_count_mat),
            rowRanges = rowranges,
            colData = ATAC_meta
        )

        motif_ix <- motifmatchr::matchMotifs(
            pwms = motif_collection,
            subject = atac_se,
            genome = genome,
            out = "scores"
        )

        list(
            ATAC_meta = ATAC_meta,
            motif_ix = motif_ix,
            ATAC_SE_obj = atac_se
        )
    }

    ### peak to motif scores
    peak_to_motif_conversion_fullset <- peak_to_motif_cell_mat_fullset(
        ATAC_meta = atac_meta,
        ATAC_count_mat = atac_count_mat
    )

    #### add bg
    fullset_bg_added_counts <- chromVAR::addGCBias(
        peak_to_motif_conversion_fullset$ATAC_SE_obj,
        genome = genome
    )

    #### deviation calculation
    BiocParallel::register(BiocParallel::MulticoreParam(workers = n_cores))
    fullset_motif_dev <- chromVAR::computeDeviations(
        object = fullset_bg_added_counts,
        annotations = peak_to_motif_conversion_fullset$motif_ix
    )

    out <- list(
        peak_to_motif_conversion_fullset = peak_to_motif_conversion_fullset,
        fullset_bg_added_counts = fullset_bg_added_counts,
        fullset_motif_dev = fullset_motif_dev
    )

    ### final save
    if(!is.null(out_rdata_path)){
        save(
            list = names(out),
            file = out_rdata_path,
            envir = list2env(out, parent = emptyenv())
        )
    }

    return(out)
}