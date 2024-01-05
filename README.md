# GenePathRx

**A gene-gene interaction model-based drug sensitivity prediction algorithm using public database information and tumor/patient-derived organoid/cell-line multi-omics data for predictions. GenePathRx was validated using the GDSC2 cancer cell line dataset as well as breast cancer patient-derived organoid data and has been demonstrated to provide high accuracy predictions.**<br>
> See: Tarapcsak et al. - GenePathRx: Model-based cancer therapy selection by linking genomic and transcriptomic tumor vulnerabilities to drug mechanisms.<br>

Due to its simple and explicit model, GenePathRx is readily applicable to any data (e.g. cell-lines, PDO, PDX, patient tumor) and can provide reliable predictions of drug sensitivity for targeted- and chemotherapy drugs alike (using DGIdb). Moreover, GenePathRx can be further utilized to understand the working mechanism of drugs through pathway analysis (using KEGG pathway database) and has the potential to be used as a computational tool for biomarker discovery.

## Running:
Before running GenePathRx make sure NetworkX is installed (version 3.2.1):
```
pip install networkx
```
GenePathRx can be run using a single line of python code (python3 is needed):
```
python3 genepathrx.py /path/to/drug_file path/to/vulnerability_file output_folder_name num_of_background_dist_values
```
<br>

## Testing:
To test if GenePathRx is running on your machine, download the genepathrx folder and run the following lines of code. Sample data of two patient samples are provided in non_responder_patient_test and responder_patient_test folders.<br>
In this example, you calculate DTCR scores for a targeted drug, M3814. M3814 is a DNA-dependent protein kinase (DNA-PK) inhibitor that has one target gene in DGIdb, PRKDC, the catalytic subunit of DNA-PK.
To run GenePathRx on the M3814 sensitive patient sample, use the following code:
```
cd genepathrx
python3 genepathrx.py ./responder_patient_test/M3814.txt ./responder_patient_test/responder_vulnerabilities.txt responder_patient_test 500
```
And to run GenePathRx on the M3814 non-responder patient sample:
```
python3 genepathrx.py ./non_responder_patient_test/M3814.txt ./non_responder_patient_test/non_responder_vulnerabilities.txt non_responder_patient_test 500
```
Each calculation take ~ 10 min to finish.
DTCR scores will be calculated using a background distribution of 500 random samples. Due to random sampling, the calculated DTCR scores might differ from DTCR scores in the sample output files.
<br>

## Input files:
You can see examples in the responder_patient_test and non_responder_patient_test folders:
- **drug_file.txt** &emsp;_- Name of the drug/drugs to be analyzed from DGIdb. Each line in the input file is one drug name_
- **vulnerability_file.txt** &emsp;_- Name of the genes identified as genomic/transcriptomic vulnerabilites (vulnerabilities are identified by the user). Each line in the input file is one gene name (NCBI)_
<br>

## Output files:
You can see examples in the sample_output_files folder:
- **TVIS_results_drugx.txt** &emsp;_- Description of interaction paths of drug x target genes with every other genes in the gene interaction graph_
- **DTCS_results_drugx.txt** &emsp;_- TVIR statistics and DTCS score of drug x_
- **DTCR_results_drugx.txt** &emsp;_- DTCR score and DTCS scores of background calculations_
<br>

## Necessary files to run GenePathRx (in genepathrx folder):
- **NCBI_gene_symbols.txt** &emsp;_- List of genes in the gene interaction network_
- **parsed_DGIDB.txt** &emsp;_- Drug-target interactions of DGIdb (v4.0)_
- **gene_interaction_network.pickle** &emsp;_- Gene-gene interaction network from KEGG pathway database_
<br>

## Contact:
**Szabolcs Tarapcsak, PhD**<br>
_**szabolcs.tarapcsak@servier.com**_<br>
Developed at:<br>
_Eccles Institute of Human Genetics_<br>
_University of Utah, Salt Lake City, UT, USA_<br>
