# GenePathRx

A gene-gene interaction model-based drug sensitivity prediction algorithm using public database information and patient/cell line multi-omics data for predictions. GenePathRx was validated using the GDSC2 cancer cell line dataset as well as breast cancer patient-derived organoid data and has been demonstrated to provide high accuracy predictions.##

Due to its simple and explicit model, GenePathRx is readily applicable to any data (e.g. cell lines, PDO, PDX) and can provide reliable predictions of drug sensitivity for targeted- and chemotherapy drugs alike. Moreover, GenePathRx can be further utilized to understand the working mechanism of drugs through pathway analysis and has the potential to be used as a computational tool for biomarker discovery. 

## Running:
python3 genepathrx.py drug_file vulnerability_file

## Output files:
TVIS_results_drugx.txt - Description of interaction paths of drug x target genes with every other genes in the gene interaction graph
FINAL_results_drugx.txt - TVIR scores and overall DTCS scores of drug x
DTCR_results_drugx.txt - Overall DTCR score and DTCS scores of background calculations

## Necessary files to run GenePathRx in the same path as genepathrx.py:
NCBI_gene_symbols.txt - List of genes in the gene interaction network
parsed_DGIDB.txt - Drug-target interactions of DGIdb (v4.0)
gene_interaction_network.pickle - Gene-gene interaction network from KEGG pathway database

## Contact:
Szabolcs Tarapcsak, PhD
szabolcs.tarapcsak@genetics.utah.edu
Eccles Institute of Human Genetics
University of Utah, Salt Lake City, UT
