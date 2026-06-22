## Overview This repository contains analysis scripts used to construct and analyze a pan-cancer single-cell multiomic atlas of 60 cancer cell lines spanning 16 tissue origins and 20 cancer types. The atlas integrates snRNA-seq and snATAC-seq profiles to characterize cancer-intrinsic regulatory programs, cell-state heterogeneity, EMT-associated regulatory features, copy-number-driven state reprogramming, and subtype-specific melanoma regulatory landscapes.

## Script descriptions

### `parallel_ATAC_dimension_reduction.py`

This script performs cell line–wise preprocessing and dimension reduction for single-cell ATAC-seq data using SnapATAC2. For each cell line, the script selects highly variable peaks, performs spectral embedding, computes UMAP coordinates, and saves the processed ATAC object.

**Input:** A full single-cell ATAC `.h5ad` object containing peak count information and cell line annotations.

**Output:** One processed SnapATAC2 `.h5ad` object per cell line, containing selected features, spectral embeddings, and UMAP coordinates.

---

### `parallel_RNA_ATAC_integration.py`

This script performs cell line–wise integration of matched single-cell RNA-seq and single-cell ATAC-seq data using SCGLUE. For each cell line, RNA and ATAC profiles are preprocessed, linked through a gene–peak guidance graph, embedded into a shared latent space, and combined for joint visualization and downstream analysis.

**Input:** A full RNA `.h5ad` object, cell line–specific ATAC `.h5ad` objects, and a gene annotation GTF file.

**Output:** For each cell line, the script outputs a combined RNA/ATAC `.h5ad` object with GLUE embeddings and UMAP coordinates, the full and highly variable feature guidance graphs, and the trained SCGLUE model object.

---

### `seacells_identification.py`

This script identifies SEACell metacells from the integrated RNA/ATAC embedding for each cell line. It uses the GLUE latent space to group cells into metacells and filters metacells to retain those with sufficient representation from both RNA and ATAC modalities.

**Input:** A cell line–specific integrated RNA/ATAC `.h5ad` object containing GLUE embeddings.

**Output:** A filtered cell metadata table listing SEACell assignments and a corresponding `.h5ad` object containing cells assigned to retained metacells.

---

### `TF_regulon_scoring.py`

This script scores predefined transcription factor regulon or gene-program signatures in single-cell RNA-seq data. It calculates per-cell gene-program scores and tests whether each program is enriched in individual cell lines compared with all other cells.

**Input:** An AnnData object and a dictionary of gene programs or TF regulons.

**Output:** A summary table containing cell line–specific program scores, background scores, score differences, nominal p-values, and FDR-adjusted p-values.

---

### `chromVAR_TF_deviation.R`

This script computes transcription factor motif deviation scores from single-cell ATAC-seq peak count matrices using chromVAR. It converts genomic peak coordinates into genomic ranges, matches peaks to motif collections, corrects for GC bias, and calculates motif deviation activity across cells.

**Input:** An ATAC cell metadata table, an ATAC peak-by-cell count matrix, and a motif collection.

**Output:** A chromVAR result object containing peak-to-motif annotations, GC-bias-corrected counts, and motif deviation scores. The results can also be saved as an `.RData` file.

---

### `signaling_activity_inference.R`

This script infers pathway activity from single-cell RNA-seq expression data using PROGENy through the `decoupleR` framework. It estimates signaling pathway activities at the single-cell level and summarizes pathway activity across cell lines.

**Input:** A gene expression matrix, cell metadata, and cell line annotations.

**Output:** A Seurat object with inferred pathway activity scores added to metadata, a long-format pathway activity table, a pathway activity matrix, a cell line–level median activity summary, and the PROGENy regulatory network used for inference.

---

### `CNV_visualization.R`

This script prepares inferCNV output for chromosome-level visualization. It maps inferCNV states to copy-number status values, adjusts genomic coordinates across chromosomes, and formats cell line ordering for heatmap-style visualization.

**Input:** An inferCNV output table, chromosome length information, and an ordered list of cell lines to visualize.

**Output:** A formatted CNV block table with cumulative genomic coordinates and a chromosome skeleton table for plotting chromosome boundaries.

---

### `CNV_EMT_association.R`

This script tests associations between gene-level copy-number variation and EMT-related phenotypes across cell lines. It performs gene-wise linear modeling between copy-number values and median EMT pseudotime, optionally evaluates the relationship between ATAC gene activity and EMT, and prepares results for Manhattan-style visualization.

**Input:** A cell line–level EMT score table, a gene-by-cell line copy-number matrix, optional gene activity matrix, chromosome length information, and a GTF annotation file.

**Output:** Tables containing CNV–EMT association statistics, optional gene activity–EMT correlation results, merged association results, and a formatted Manhattan plot table with genomic coordinates.

---

### `CNV_drug_response_association.R`

This script evaluates associations between gene-level copy-number alterations and drug response across cancer cell lines. For each gene and drug, cell lines are grouped by amplification, deletion, or diploid status, and drug sensitivity is tested using ANOVA and linear modeling.

**Input:** A list of genes of interest, a gene-by-cell line CNV matrix, and a drug response matrix containing IC50 values across cell lines.

**Output:** A gene-by-drug association table containing ANOVA p-values, linear model estimates, CNV-term p-values, and median drug response values for amplified, diploid, and deleted groups.

---

### `cell_phenotype_mapping.py`

This script maps query cancer cell line profiles onto a reference single-cell atlas. It preprocesses the reference and query datasets, identifies shared highly variable genes, builds reference and query embeddings, transfers reference cell-type labels to query cells using Scanpy ingest, and summarizes the mapped cell-type composition for each cell line.

**Input:** A reference `.h5ad` object, a query `.h5ad` object, selected cell lines to map, and metadata columns specifying reference cell types and query cell line labels.

**Output:** A combined AnnData object for visualization of reference and query cells, and a percentage table summarizing the inferred reference cell-type composition of each mapped cell line.
