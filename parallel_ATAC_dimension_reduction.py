from multiprocessing import Pool
from multiprocessing import *
import sys
from functools import partial
import snapatac2 as snap
import pandas as pd
import numpy as np
import os
import re
from collections import Counter
import copy
import scanpy as sc
import matplotlib as plt

def multicore_cell_line_wise_single_cell_ATAC_preprocessing(full_ATAC_sc_obj_path, parent_folder, cell_line_name_col="standardized_cell_line_names", cores=4):
    
    ###make the folder
    os.mkdir(parent_folder)
    
    ###load the full data
    full_ATAC_sc_obj = snap.read(full_ATAC_sc_obj_path, backed=None)
    
    ###get the cell line name list
    cell_line_name_list = list(full_ATAC_sc_obj.obs[cell_line_name_col].unique())
    
    ###parallelization
    p = Pool(processes = int(cores))
    print("Processing core number: ", cores)
    
    func = partial(single_cell_ATAC_cell_line_preprocessing, full_ATAC_sc_obj = full_ATAC_sc_obj, out_folder_path=parent_folder)
    result = p.map(func, cell_line_name_list)
    p.close()
    p.join()
    
    return("All done!")
    

def single_cell_ATAC_cell_line_preprocessing(cell_line_name, full_ATAC_sc_obj, out_folder_path, n_top_peaks=300000):
    
    ###define the out folder
    out_folder_path_full = out_folder_path + "/snapatac2." + cell_line_name + ".spectral.h5ad"
    
    ###split the cell line ATAC obj
    cell_line_ATAC_sc_obj = full_ATAC_sc_obj[full_ATAC_sc_obj.obs["standardized_cell_line_names"] == cell_line_name,:].copy()
    
    ###perform dimension reduction
    snap.pp.select_features(cell_line_ATAC_sc_obj, n_features=n_top_peaks)
    snap.tl.spectral(cell_line_ATAC_sc_obj)
    snap.tl.umap(cell_line_ATAC_sc_obj)
    
    ###export 
    cell_line_ATAC_sc_obj.write_h5ad(out_folder_path_full)
    
    return("SnapATAC2 preprocessing on cell line: " + cell_line_name + ", is done.")

#####execute the function
if __name__ == "__main__":
    multicore_cell_line_wise_single_cell_ATAC_preprocessing(full_ATAC_sc_obj_path="/the/path/pan_cancer_snapatac2_peak_counts.h5ad", 
                                                            parent_folder="/out/path")