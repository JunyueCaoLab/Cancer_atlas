import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


def compute_program_enrichment_fast(
        adata,
        gene_program_dict,
        cell_line_col="cell_lines"
):

    ### score gene programs
    for program_name, gene_list in gene_program_dict.items():
        genes = [g for g in gene_list if g in adata.var_names]
        sc.tl.score_genes(adata, genes, score_name=program_name)

    df = adata.obs

    cell_lines = df[cell_line_col].unique()
    programs = list(gene_program_dict.keys())

    results = []

    for program in programs:

        program_scores = df[program]

        for cl in cell_lines:

            mask = df[cell_line_col] == cl

            group1 = program_scores[mask]
            group2 = program_scores[~mask]

            stat, p = mannwhitneyu(
                group1,
                group2,
                alternative="two-sided"
            )

            median1 = group1.median()
            median2 = group2.median()

            results.append({
                "cell_line": cl,
                "program": program,
                "median_score": median1,
                "background_median": median2,
                "score_diff": median1 - median2,
                "p_value": p
            })

    results_df = pd.DataFrame(results)

    ### FDR correction
    results_df["fdr"] = multipletests(
        results_df["p_value"],
        method="fdr_bh"
    )[1]

    results_df = results_df.sort_values("fdr")

    return results_df