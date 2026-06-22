import os
import scanpy as sc
import pandas as pd
import numpy as np
import scglue 
import copy
import math
import gzip
import re
import pybedtools
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.io import mmread
from scipy.sparse import csr_matrix
from itertools import chain
import anndata as ad
import itertools
import networkx as nx
import seaborn as sns
from matplotlib import rcParams
import sys

def multicore_cell_line_wise_RNA_ATAC_integration_glue(main_folder, full_RNA_sc_obj_path, cell_line_ATAC_w_DM_sc_obj_path, cores=1):
    
    ###create a parental folder for every cell line
    os.mkdir(main_folder)
    
    ###read in the full RNA obj
    full_RNA_sc_obj = sc.read_h5ad(full_RNA_sc_obj_path)
    
    ###get the cell line name list
    cell_line_name_list = list(full_RNA_sc_obj.obs["standardized_cell_line_names"].unique())
    print(str(len(cell_line_name_list)) + " cell lines are going to be processed.")
    
    for each_cell_line in cell_line_name_list:
        single_cell_line_RNA_ATAC_integration_glue(cell_line_name=each_cell_line,
                                                   full_RNA_sc_obj = full_RNA_sc_obj, 
                                                   cell_line_ATAC_w_DM_sc_obj_path = cell_line_ATAC_w_DM_sc_obj_path, 
                                                   output_path=main_folder)
    
    return("All done!")

def single_cell_line_RNA_ATAC_integration_glue(cell_line_name, full_RNA_sc_obj, cell_line_ATAC_w_DM_sc_obj_path, output_path, gtf_path="gencode.v38.primary_assembly.annotation.gtf.gz"):
    
    print("Start processing " + cell_line_name)
    
    ###step1: create a folder under the parental folder
    cell_line_spec_folder=output_path+"/scglue."+cell_line_name
    os.mkdir(cell_line_spec_folder)
    
    ###step2: get the single cell line RNA and ATAC data
    cell_line_RNA_sc_obj = full_RNA_sc_obj[full_RNA_sc_obj.obs["standardized_cell_line_names"] == cell_line_name,:].copy()
    cell_line_ATAC_sc_obj = sc.read_h5ad(cell_line_ATAC_w_DM_sc_obj_path+"/snapatac2." + cell_line_name + ".spectral.h5ad")
    
    ###step3: preprocess RNA obj with a raw gene count matrix in X
    sc.pp.highly_variable_genes(cell_line_RNA_sc_obj, n_top_genes=3000, flavor="seurat_v3")
    cell_line_RNA_sc_obj.layers["counts"] = cell_line_RNA_sc_obj.X.copy()
    
    sc.pp.normalize_total(cell_line_RNA_sc_obj)
    sc.pp.log1p(cell_line_RNA_sc_obj)
    cell_line_RNA_sc_obj = cell_line_RNA_sc_obj[:, cell_line_RNA_sc_obj.var.highly_variable].copy()
    sc.tl.pca(cell_line_RNA_sc_obj, n_comps=30, svd_solver="auto")
    
    scglue.data.get_gene_annotation(cell_line_RNA_sc_obj,
                                    gtf=gtf_path,
                                    gtf_by="gene_name", 
                                    var_by="gene_symbol")
    cell_line_RNA_sc_obj = cell_line_RNA_sc_obj[:,np.isnan(cell_line_RNA_sc_obj.var["chromStart"])==False].copy()
    
    ###step4: preprocess ATAC obj with a raw peak count matrix in X and X_spectral representing dimension reduction
    split = cell_line_ATAC_sc_obj.var_names.str.split(r"[:-]")
    cell_line_ATAC_sc_obj.var["chrom"] = split.map(lambda x: x[0])
    cell_line_ATAC_sc_obj.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
    cell_line_ATAC_sc_obj.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
    
    ###step5: scglue training
    guidance = scglue.genomics.rna_anchored_guidance_graph(cell_line_RNA_sc_obj, cell_line_ATAC_sc_obj)
    scglue.graph.check_graph(guidance, [cell_line_RNA_sc_obj, cell_line_ATAC_sc_obj])
    
    scglue.models.configure_dataset(cell_line_RNA_sc_obj, "NB", use_highly_variable=True,use_layer="counts", use_rep="X_pca")
    scglue.models.configure_dataset(cell_line_ATAC_sc_obj, "NB", use_highly_variable=True,use_rep="X_spectral")
    guidance_hvf = guidance.subgraph(chain(cell_line_RNA_sc_obj.var.query("highly_variable").index,cell_line_ATAC_sc_obj.var.query("highly_variable").index)).copy()
    
    glue = scglue.models.fit_SCGLUE({"rna": cell_line_RNA_sc_obj, "atac": cell_line_ATAC_sc_obj}, 
                                    guidance_hvf,
                                    fit_kws={"directory": "glue"})
    
    ###step6: combine RNA and ATAC modalities for co-embedding
    cell_line_RNA_sc_obj.obsm["X_glue"] = glue.encode_data("rna", cell_line_RNA_sc_obj)
    cell_line_ATAC_sc_obj.obsm["X_glue"] = glue.encode_data("atac", cell_line_ATAC_sc_obj)
    combined = ad.concat([cell_line_RNA_sc_obj, cell_line_ATAC_sc_obj])
    
    combined.obs["modality"] = "RNA"
    combined.obs.loc[np.array([("ATAC_Cells" in x) for x in combined.obs.index]),"modality"] = "ATAC"
    
    combined.obs["sample"]=combined.obs["sample"].astype("str")
    
    sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine", n_neighbors=40)
    sc.tl.umap(combined, min_dist=0.3)
    
    ###step7: export related files
    combined.write_h5ad(cell_line_spec_folder+"/scglue."+cell_line_name+".combined.h5ad")
    nx.write_graphml(guidance, cell_line_spec_folder+"/scglue."+cell_line_name+".guidance.graphml.gz")
    nx.write_graphml(guidance_hvf, cell_line_spec_folder+"/scglue."+cell_line_name+".guidance_hvf.graphml.gz")
    glue.save(cell_line_spec_folder+"/scglue."+cell_line_name+".glue.drill")
    
    print("Done processing " + cell_line_name)
    
    return(0)

#####execute the function
if __name__ == "__main__":
    multicore_cell_line_wise_RNA_ATAC_integration_glue(main_folder="/out/path", 
                                                       full_RNA_sc_obj_path="RNA_obj.h5ad", 
                                                       cell_line_ATAC_w_DM_sc_obj_path="/ATAC_DM_objs/path")