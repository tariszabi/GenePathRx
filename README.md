# GenePathRx

**A gene-gene interaction model-based drug sensitivity prediction algorithm using public database information and tumor/patient-derived organoid/cell-line multi-omics data for predictions. GenePathRx was validated using the GDSC2 cancer cell line dataset as well as breast cancer patient-derived organoid data and has been demonstrated to provide high accuracy predictions.**<br>
> See: Tarapcsak et al. - GenePathRx: Model-based cancer therapy selection by linking genomic and transcriptomic tumor vulnerabilities to drug mechanisms.<br>

Due to its simple and explicit model, GenePathRx is readily applicable to any data (e.g. cell-lines, PDO, PDX, patient tumor) and can provide reliable predictions of drug sensitivity for targeted- and chemotherapy drugs alike (using DGIdb). Moreover, GenePathRx can be further utilized to understand the working mechanism of drugs through pathway analysis (using KEGG pathway database) and has the potential to be used as a computational tool for biomarker discovery.

## Running:
Before running GenePathRx make sure NetworkX is installed:
```
pip install networkx
```
To run GenePathRx in the terminal, python3 is needed:
```
python3 genepathrx.py /path/to/drug_file path/to/vulnerability_file /output/folder/path
```
To test if GenePathRx is running on your machine, download the genepathrx folder on your machine and run the following lines of code:
```
cd genepathrx
python3 genepathrx.py ./test/test_drug_file.txt .test/test_vulnerability_file.txt ./test
```

## Output files:
- **TVIS_results_drugx.txt** Description of interaction paths of drug x target genes with every other genes in the gene interaction graph
- **FINAL_results_drugx.txt** TVIR scores and overall DTCS scores of drug x
- **DTCR_results_drugx.txt** Overall DTCR score and DTCS scores of background calculations

## Necessary files to run GenePathRx in the same folder as genepathrx.py:
| FILE                   | DESCRIPTION                                                      |
|------------------------|------------------------------------------------------------------|
| **NCBI_gene_symbols.txt** | List of genes in the gene interaction network|
| **parsed_DGIDB.txt**| Drug-target interactions of DGIdb (v4.0)|
| **gene_interaction_network.pickle** | Gene-gene interaction network from KEGG pathway database|

## Contact:
**Szabolcs Tarapcsak, PhD**<br>
_**szabolcs.tarapcsak@genetics.utah.edu**_<br>
_Eccles Institute of Human Genetics_<br>
_University of Utah, Salt Lake City, UT, USA_<br>
