#!/usr/bin/env

# Copyright 2024 by Szabolcs Tarapcsak PhD
# Developed at Eccles Institute of Human Genetics, University of Utah,
# Salt Lake City, Utah, USA.
# All rights reserved.
# GenePathRx, A gene-gene interaction model-based drug sensitivity prediction 
# algorithm using public database information and tumor/patient-derived 
# organoid/cell-line multi-omics data for predictions.

import networkx as nx
import sys
import random
from functools import reduce
from statistics import median
import time
import pickle
import os.path

start = time.time()

print("\n*** GenePathRx ***\n")
print("Starting analysis\n")
cwd = os.path.dirname(os.path.abspath(__file__))
argument1 = sys.argv[1]
argument2 = sys.argv[2]
argument3 = sys.argv[3]
argument4 = sys.argv[4]
with open('gene_interaction_network.pickle', 'rb') as f:
    G = pickle.load(f)

# Opening files with the lists of drugs and identified vulnerability genes and the KEGG gene network file:
print("Step 1: Loading information from input files")
drug_file = open(str(argument1), "r")
drugs = drug_file.read().splitlines()
print("\tAnalysis is performing on drugs: " + str(drugs))
print("\tSample ID: " + str(argument3))
drug_file.close()
genes_file = open(str(argument2), "r")
pre_genes = genes_file.read().splitlines()
genes_file.close()
genes = []
for g in pre_genes:
    if g in G.nodes():
        genes.append(g)
print("\tThe sample has the following vulnerability genes within the gene network: " + str(genes) + "\n")

# Parsing DGIdb to identify target genes of each drug and filtering out targets not found in the network:
print("Step 2: Identifying drug target genes")
for listed_drug in drugs:
    DGIDB_file = open(cwd + "/" + "parsed_DGIDB.txt", "r")
    drug_target_pairs = DGIDB_file.read().splitlines()
    DGIDB_file.close()
    drug_name = str(listed_drug + "\t")
    pre_target_genes = []
    for line in drug_target_pairs:
        if drug_name in line:
            pre_target_genes.append(line.split("\t")[1] + "," + line.split("\t")[2])
    target_genes = []
    for k in pre_target_genes:
        if k.split(",")[0] in G.nodes():
            target_genes.append(k)
    print("\tTarget genes and associated DGIdb Interaction scores of " + str(listed_drug) + " are: "
          + str(target_genes) + "\n")

# Generation of an input file for Target-Vulnerability Interaction Score (TVIS) calculation using NCBI genes:
    print("Step 3: Generation of drug target gene - NCBI-gene pairs for each gene in gene network graph\n")
    print("Step 4: Calculating Target-Vulnerability Interaction Scores (TVIS)\n")
    NCBI_genes = open(cwd + "/" + "NCBI_gene_symbol.txt", "r")
    NCBI_gene_pre_list = NCBI_genes.read().splitlines()
    NCBI_genes.close()
    NCBI_gene_list = []
    for ge in NCBI_gene_pre_list:
        if ge in G.nodes():
            NCBI_gene_list.append(ge)
    target_NCBI_gene_pairs = []
    for x in target_genes:
        target = x.split(",")[0]
        DGIDB_score = x.split(",")[1]
        for NCBI_gene in NCBI_gene_list:
            target_NCBI_gene_pairs.append(target + "\t" + NCBI_gene.replace("\n", "") + "\t" + DGIDB_score)

# Target-Vulnerability Interaction Score (TVIS) calculation based on provided pickle
# file (stores the gene network created from KEGG pathway database):
    product = lambda x: reduce(lambda a, b: a*b, x)

    def path_score(G, path):
        return 1.0 / product([G.degree[g] - 1 for g in path[1:]]) * pow(2, -(len(path) - 2))

    def weighted_paths_score(G, paths):
        return sum([path_score(G, p) for p in paths])

    calculated_scores = []
    for tv_pair in target_NCBI_gene_pairs:
        pair = tv_pair.rstrip().split("\t")
        if pair[0] == pair[1]:
            calculated_scores.append(pair[0] + "\t" + pair[1] + "\t" + str(0) + "\t" + str(1) + "\t" + str(1)
                                     + "\t" + pair[2])
        elif nx.has_path(G, pair[0], pair[1]):
            shortest_paths = [p for p in nx.all_shortest_paths(G, pair[0], pair[1])]
            shortest_path_len = len(shortest_paths[0]) - 1
            shortest_paths_count = len([p for p in shortest_paths])
            normalized_score = weighted_paths_score(G, shortest_paths)
            calculated_scores.append(pair[0] + "\t" + pair[1] + "\t" + str(shortest_path_len) + "\t" +
                                     str(shortest_paths_count) + "\t" + str(normalized_score) + "\t" + pair[2])
        if not nx.has_path(G, pair[0], pair[1]):
            calculated_scores.append(pair[0] + "\t" + pair[1] + "\t" + "0" + "\t" + "0" + "\t" + "0" + "\t" + pair[2])

# Writing the calculated results into a TVIS_result_"name of the drug".txt file:
    print("Step 5: Saving TVIS scores in TVIS results file for drug " + str(listed_drug) + "\n")
    output_file = open(cwd + "/" + str(argument3) + "/" + 'TVIS_results_{0}.txt'.format(listed_drug), "w")
    output_file.write("Target_gene" + "\t" + "Vulnerability_gene" + "\t" + "Length_of_shortest_path" + "\t"
                      + "Number_of_shortest_paths" + "\t" + "TVIS_score" + "\t" + "DGIdb_interaction_score" + "\n")
    for line in calculated_scores:
        output_file.write(line + "\n")
    output_file.close()

# Parsing data needed for TVIR score calculation of the drug-gene pair
# using data in TVIS_result_"name of the drug".txt file:
    target_vulnerability_pairs = []
    for target_gene in target_genes:
        for vulnerability in genes:
            target_vulnerability_pairs.append(target_gene.split(",")[0] + "\t" + vulnerability + "\t"
                                              + str(target_gene.split(",")[1]))
    print("Step 6: Calculating Target-Vulnerability Interaction Rank (TVIR) scores\n")
    TVIR_scores = []
    for target_vulnerability_pair in target_vulnerability_pairs:
        target = target_vulnerability_pair.split("\t")[0]
        gene = target_vulnerability_pair.split("\t")[1]
        int_score = float(target_vulnerability_pair.split("\t")[2])
        input_file = open(cwd + "/" + str(argument3) + "/" + 'TVIS_results_{0}.txt'.format(listed_drug), "r")
        input_data = input_file.read().splitlines()
        input_file.close()
        target_background = []
        for line in input_data[1:]:
            if target in line:
                target_background.append(line)
        data = {}
        for line in target_background:
            line = line.split("\t")
            data[line[1]] = float(line[4])

# Percentile (Target-Vulnerability Interaction Rank (TVIR)) calculation:
        ls_genes = []
        for e in data.values():
            if e <= data[gene]:
                ls_genes.append(e)
        percentile = (len(ls_genes) / len(data.values()))
        TVIR_scores.append(target + "\t" + gene + "\t" + str(percentile) + "\t" + str(percentile*int_score))

# Drug-vulnerability correlation score (DVCS)) and TVIR statistics calculation:
    print("Step 7: Calculating Drug-Tumor Correlation Scores (DTCS)\n")
    DTCS_score_Data = []
    new_DTCS_score_Data = []
    for item in TVIR_scores:
        DTCS_score_Data.append(float(item.split("\t")[2]))
    for element in TVIR_scores:
        new_DTCS_score_Data.append(float(element.split("\t")[3]))
    summ_DTCS = sum(DTCS_score_Data)
    summ_new_DTCS = sum(new_DTCS_score_Data)
    average_TVIR = summ_DTCS / len(target_vulnerability_pairs)
    median_TVIR = median(DTCS_score_Data)

# Writing results (TVIR statistics and DTCS Score) into "DTCS_result_"name of the drug".txt file:
    print("Step 8: Saving DTCS scores in DTCS results file for drug " + str(listed_drug) + "\n")
    output_file_2 = open(cwd + "/" + str(argument3) + "/" + 'DTCS_result_{0}.txt'.format(listed_drug), "w")
    output_file_2.write("Calculated Drug-Tumor Correlation Score (DTCS) for " + str(listed_drug) + ":" + "\n" + "\n")
    output_file_2.write("DTCS Score: " + str(summ_DTCS) + "\n" + "DGIdb Interaction Score normalized DTCS Score: "
                        + str(summ_new_DTCS) + "\n" + "Average TVIR Score: " + str(average_TVIR) + "\n"
                        + "Median TVIR Score: " + str(median_TVIR) + "\n" + "\n" +
                        "Calculated TVIR scores for " + str(listed_drug) + ":" + "\n" + "\n")
    output_file_2.write("Target_gene" + "\t" + "Vulnerability_gene" + "\t" + "TVIR_Score" + "\t"
                        + "normalized_TVIR_Score" + "\n")
    for data_point in TVIR_scores:
        output_file_2.write(data_point + "\n")
    output_file_2.close()


# 2nd Major Part of GenePathRx algorithm: Calculating background distributions for DTCR score calculations:
# A user-specified number of hypothetical samples are created, with same number of randomly selected vulnerability genes
# as the sample, and background DTCS scores are calculated for these random samples.
# Identifying drug targets and random vulnerabilities:
    print("Step 9: Calculating DTCS scores using " + str(argument4) + " random set of vulnerabilities for "
          + listed_drug)
    final_data = []
    counter = 1
    while counter < int(argument4) + 1:
        print("Number of rounds: " + str(counter) + "/" + str(argument4) + " completed")
        random_genes = random.sample(NCBI_gene_list, len(genes))

        target_vulnerability_pairs_2 = []
        for target_gene in target_genes:
            for vulnerability in random_genes:
                target_vulnerability_pairs_2.append(target_gene.split(",")[0] + "\t"
                                                    + vulnerability + "\t" + str(target_gene.split(",")[1]))

        TVIR_scores = []
        for target_vulnerability_pair in target_vulnerability_pairs_2:
            target = target_vulnerability_pair.split("\t")[0]
            gene = target_vulnerability_pair.split("\t")[1]
            int_score = float(target_vulnerability_pair.split("\t")[2])
            input_file = open(cwd + "/" + str(argument3) + "/" + "TVIS_results_{0}.txt".format(listed_drug), "r")
            input_data = input_file.read().splitlines()
            input_file.close()
            target_background = []
            for line in input_data[1:]:
                if target in line:
                    target_background.append(line)
            data = {}
            for line in target_background:
                line = line.split("\t")
                data[line[1]] = float(line[4])

# Target-Vulnerability Interaction Rank (TVIR) calculation:
            ls_genes = []
            for e in data.values():
                if e <= data[gene]:
                    ls_genes.append(e)
            percentile = (len(ls_genes) / len(data.values()))
            TVIR_scores.append(target + "\t" + gene + "\t" + str(percentile) + "\t" + str(percentile * int_score))
            
        # Drug-Tumor Correlation Score (DTCS)) calculation:
        DTCS_score_Data = []
        new_DTCS_score_Data = []
        for item in TVIR_scores:
            DTCS_score_Data.append(float(item.split("\t")[2]))
        for element in TVIR_scores:
            new_DTCS_score_Data.append(float(element.split("\t")[3]))
        summ_DTCS = sum(DTCS_score_Data)
        summ_new_DTCS = sum(new_DTCS_score_Data)
        final_data.append(str(listed_drug) + "\t" + str(counter) + "\t" + str(summ_DTCS) + "\t" + str(summ_new_DTCS))
        counter += 1

# Comparing previously calculated DTCS values to the background distribution:
    print("\nStep 10: Calculating Drug-Tumor Correlation Rank (DTCR) value for " + listed_drug + "\n")
    summ_DTCS_list = []
    summ_new_DTCS_list = []
    for line in final_data:
        summ_DTCS_list.append(float(line.split("\t")[2]))
        summ_new_DTCS_list.append(float(line.split("\t")[3]))
    real_data_file = open(cwd + "/" + str(argument3) + "/" + "DTCS_result_{0}.txt".format(listed_drug), "r")
    real_data_lines = real_data_file.read().splitlines()
    real_data_file.close()
    real_sum_DTCS = float(real_data_lines[2].split(": ")[1])
    real_sum_new_DTCS = float(real_data_lines[3].split(": ")[1])
    summ_DTCS_genes = []
    for e in summ_DTCS_list:
        if e <= real_sum_DTCS:
            summ_DTCS_genes.append(e)
    summ_DTCS_percentile = (len(summ_DTCS_genes) / len(summ_DTCS_list))
    summ_new_DTCS_genes = []
    for e in summ_new_DTCS_list:
        if e <= real_sum_new_DTCS:
            summ_new_DTCS_genes.append(e)
    summ_new_DTCS_percentile = (len(summ_new_DTCS_genes) / len(summ_new_DTCS_list))
    normalized_DTCS_genes = []

# Writing final results into "DTCR_result_"name of the drug".txt file:
    print("Step 11: Saving DTCR scores in DTCR results file for drug " + str(listed_drug) + "\n")
    output_file_2 = open(cwd + "/" + str(argument3) + "/" + 'DTCR_results_{0}.txt'.format(listed_drug), "w")
    output_file_2.write("DTCR Score: " + str(summ_DTCS_percentile) + "\n" + "\n")
    output_file_2.write("drug" + "\t" + "#_of_rounds" + "\t" + "DTCS_score" + "\t" + "normalized_DTCS_score" + "\n")
    for line in final_data:
        output_file_2.write(line + "\n")
    output_file_2.close()

end = time.time()
res = end - start
final_res = res / 60
print("Analysis COMPLETE")
print("Analysis finished in", final_res, 'minutes\n')
print("\n*** GenePathRx ***\n")
print("Developed at Eccles Institute of Human Genetics, University of Utah,")
print("Salt Lake City, Utah, USA by Szabolcs Tarapcsak PhD")
print("All rights reserved. ")
print("2024\n")
