import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

def map_cell_lines_to_reference_atlas(
    ref_h5ad_path: str,
    query_h5ad_path: str,
    cell_lines_to_map,
    ref_cell_type_col: str = "cell_type",
    query_cell_line_col: str = "standardized_cell_line_names",
    query_counts_layer: str = "counts",
    ref_feature_col: str = "feature_name",
    query_gene_symbol_col: str = "gene_symbol",
    min_cells_ref: int = 350,
    n_top_genes: int = 4000,
    n_pcs: int = 30,
    do_umap: bool = True
):

    ### load + prep reference
    ref = sc.read_h5ad(ref_h5ad_path)

    ref_raw = ref.raw.to_adata()
    ref_raw.layers["raw_counts"] = ref_raw.X.copy()

    ### filter + HVG on raw counts
    sc.pp.filter_genes(ref_raw, min_cells=min_cells_ref)
    sc.pp.highly_variable_genes(ref_raw, layer="raw_counts", flavor="seurat_v3", n_top_genes=n_top_genes)

    ref_hvg = ref_raw[:, ref_raw.var["highly_variable"]].copy()

    ### clean gene symbols on reference
    ref_hvg.var["clean_gene_symbol"] = ref_hvg.var[ref_feature_col].astype(str).map(lambda x: x.split("_", 1)[0])
    ref_hvg = ref_hvg[:, ~ref_hvg.var["clean_gene_symbol"].duplicated()].copy()
    ref_hvg.var_names = ref_hvg.var["clean_gene_symbol"].astype(str).values
    ref_hvg.var_names_make_unique()

    ### normalize/log for reference embedding+ingest
    sc.pp.normalize_total(ref_hvg)
    sc.pp.log1p(ref_hvg)
    ref_hvg.layers["log1p_norm_counts"] = ref_hvg.X.copy()

    ### load + prep query
    q = sc.read_h5ad(query_h5ad_path)

    if query_gene_symbol_col not in q.var.columns:
        raise KeyError(f"Query `.var` missing '{query_gene_symbol_col}'. Available: {list(q.var.columns)[:20]} ...")

    q.var[query_gene_symbol_col] = q.var[query_gene_symbol_col].astype(str)
    q = q[:, ~q.var[query_gene_symbol_col].duplicated()].copy()
    q.var_names = q.var[query_gene_symbol_col].astype(str).values
    q.var_names_make_unique()

    if query_cell_line_col not in q.obs.columns:
        raise KeyError(f"Query `.obs` missing '{query_cell_line_col}'. Available: {list(q.obs.columns)[:20]} ...")

    ### subset cell lines
    q_sub = q[q.obs[query_cell_line_col].isin(list(cell_lines_to_map)), :].copy()
    q_sub.layers["raw_counts"] = q_sub.layers[query_counts_layer].copy()

    ### normalize/log for query
    sc.pp.normalize_total(q_sub)
    sc.pp.log1p(q_sub)
    q_sub.layers["log1p_norm_counts"] = q_sub.X.copy()

    ### intersect genes
    intersected = sorted(set(ref_hvg.var_names) & set(q_sub.var_names))
    if len(intersected) == 0:
        raise ValueError("No intersected genes between reference HVGs and query genes after symbol cleaning.")

    ref_hvg = ref_hvg[:, intersected].copy()
    q_sub = q_sub[:, intersected].copy()

    ### build ref embedding
    sc.tl.pca(ref_hvg, n_comps=n_pcs)
    sc.pp.neighbors(ref_hvg)
    if do_umap:
        sc.tl.umap(ref_hvg)

    ### build query PCA
    sc.tl.pca(q_sub, n_comps=n_pcs)
    sc.pp.neighbors(q_sub)
    if do_umap:
        sc.tl.umap(q_sub)

    ### ingest labels
    sc.tl.ingest(adata=q_sub, adata_ref=ref_hvg, obs=ref_cell_type_col)

    ### concat for joint visualization
    concat = ad.concat([ref_hvg, q_sub], label="batch", keys=["ref", "new"])

    ### add a unified cell-line column for plotting
    concat.obs["cell_line_names"] = "Ref_atlas"
    concat.obs.loc[q_sub.obs_names, "cell_line_names"] = q_sub.obs[query_cell_line_col].astype(str).values

    ### percent mapped cell type per cell line
    pct_df = pd.crosstab(
        q_sub.obs[query_cell_line_col].astype(str),
        q_sub.obs[ref_cell_type_col].astype(str),
        normalize="index"
    ) * 100.0

    return concat, pct_df
