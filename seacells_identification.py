import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import os
import sys
from scipy.io import mmread
from scipy.sparse import csc_matrix

def sea_cell_assignment(glue_combined_obj_path, seacell_assigned_out_path, n_cell_per_metacell=75):
    
    ###load the integrated obj
    integrated_obj = sc.read_h5ad(glue_combined_obj_path)
    cell_line_name = integrated_obj.obs["standardized_cell_line_names"].unique()[0]
    
    ###determine the output path
    cell_line_spec_out_path = seacell_assigned_out_path+"/cell_line." + cell_line_name
    
    ###estimate the number of meta cells
    n_meta_cells = round(integrated_obj.shape[0]/75)
    
    ###start building meta cells
    model = SEACells.core.SEACells(integrated_obj, 
                              build_kernel_on='X_glue', 
                              n_SEACells=n_meta_cells, 
                              n_waypoint_eigs=10,
                              convergence_epsilon = 1e-5)
    
    model.construct_kernel_matrix()
    M = model.kernel_matrix
    model.initialize_archetypes()
    model.fit(min_iter=10, max_iter=70)
    
    ###visualize statistics and seacell centers
    plt.figure(figsize=(3,2))
    sns.distplot((model.A_.T > 0.1).sum(axis=1), kde=False)
    plt.title(f'Non-trivial (> 0.1) assignments per cell')
    plt.xlabel('# Non-trivial SEACell Assignments')
    plt.ylabel('# Cells')
    plt.show()

    plt.figure(figsize=(3,2))
    b = np.partition(model.A_.T, -5)    
    sns.heatmap(np.sort(b[:,-5:])[:, ::-1], cmap='viridis', vmin=0)
    plt.title('Strength of top 5 strongest assignments')
    plt.xlabel('$n^{th}$ strongest assignment')
    plt.show()
    
    labels,weights = model.get_soft_assignments()
    
    SEACells.plot.plot_2D(integrated_obj, key='X_umap', colour_metacells=True)
    
    ###check ncells
    ncell_in_meta_cells = list(integrated_obj.obs['SEACell'].value_counts())
    
    print("The size of meta cells are: " + str(ncell_in_meta_cells), file=sys.stderr)
    
    ###check the balance of RNA and ATAC modalities
    all_seacell_names = list(integrated_obj.obs['SEACell'].unique())
    #n_RNA_cells_col = []
    #n_ATAC_cells_col = []
    
    filtered_seacell_names = []
    
    for each_meta_cell in all_seacell_names:
        modality_cells = integrated_obj.obs.loc[integrated_obj.obs['SEACell'] == each_meta_cell,"modality"]
        RNA_ncells = sum(modality_cells == "RNA")
        ATAC_ncells = sum(modality_cells == "ATAC")
        
        if RNA_ncells/len(modality_cells) <= 0.8 and RNA_ncells/len(modality_cells) >= 0.2 and RNA_ncells >= 5 and ATAC_ncells >= 5:
            filtered_seacell_names.append(each_meta_cell)
    
    n_initial_meta_cells = len(all_seacell_names)
    n_meta_cells_keep = len(filtered_seacell_names)
    
    print("Meta cells identified: " + str(n_initial_meta_cells), file=sys.stderr)
    print("Meta cells kept: " + str(n_meta_cells_keep), file=sys.stderr)
            
    ###reformat meta data and obj
    ready_for_pseudobulk_cell_meta = integrated_obj.obs.copy()
    ready_for_pseudobulk_cell_meta = ready_for_pseudobulk_cell_meta.loc[ready_for_pseudobulk_cell_meta["SEACell"].isin(filtered_seacell_names), :]
    
    os.mkdir(cell_line_spec_out_path)
    
    ready_for_pseudobulk_cell_meta.to_csv(cell_line_spec_out_path + "/cell_line_spec.filtered_seacell.csv")
    integrated_obj.write_h5ad(cell_line_spec_out_path + "/cell_line_spec.filtered_seacell.h5ad")
    