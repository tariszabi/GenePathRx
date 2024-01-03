#!/usr/bin/env

# Copyright 2024 by Szabolcs Tarapcsak, Department of Human Genetics, University of Utah.
# All rights reserved.
# GenePathRx, A gene-gene interaction model-based drug sensitivity prediction 
# algorithm using public database information and tumor/patient-derived 
# organoid/cell-line multi-omics data for predictions.

import networkx as nx
from itertools import combinations
import sys
import random
from functools import reduce
from statistics import median
import time

start = time.time()

argument1 = sys.argv[1]
argument2 = sys.argv[2]
argument3 = sys.argv[3]
G = nx.read_gpickle("gene_interaction_network.pickle")

# Opening files with the lists of drugs and identified vulnerability genes and network file:
print("Step 1: Loading information from input files")
drug_file = open(str(argument1), "r")
drugs = drug_file.read().splitlines()
drug_file.close()
genes_file = open(str(argument2), "r")
pre_genes = genes_file.read().splitlines()
genes_file.close()
genes = []
for g in pre_genes:
    if g in G.nodes():
        genes.append(g)

# Parsing DGIdb to identify target genes of each drug and filtering out targets not found in the network:
print("Step 2: Identifying drug target genes")
for listed_drug in drugs:
    DGIDB_file = open("parsed_DGIDB.txt", "r")
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

# Generation of an input file for target-vulnerability interaction score (TVIS) calculation using NCBI genes:
    print("Step 3: Generation of target gene - NCBI gene pairs")
    print("Step 4: Target-vulnerability interaction score (TVIS) calculation")
    NCBI_genes = open("NCBI_gene_symbol.txt", "r")
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

# Target-vulnerability interaction score (TVIS) calculation based on "graph.pickle" file previously generated (e.g. KEGG):
    product = lambda x: reduce(lambda a, b: a*b, x)

    def path_score(G, path):
        return 1.0 / product([G.degree[g] - 1 for g in path[1:]]) * pow(2, -(len(path) - 2))

    def weighted_paths_score(G, paths):
        return sum([path_score(G, p) for p in paths])

    calculated_scores = []
    for tv_pair in target_NCBI_gene_pairs:
        pair = tv_pair.rstrip().split("\t")
        if pair[0] == pair[1]:
            calculated_scores.append(pair[0] + "\t" + pair[1] + "\t" + str(0) + "\t" + str(1) + "\t" + str(1) + "\t" + pair[2])
        elif nx.has_path(G, pair[0], pair[1]):
            shortest_paths = [p for p in nx.all_shortest_paths(G, pair[0], pair[1])]
            shortest_path_len = len(shortest_paths[0]) - 1
            shortest_paths_count = len([p for p in shortest_paths])
            normalized_score = weighted_paths_score(G, shortest_paths)
            calculated_scores.append(pair[0] + "\t" + pair[1] + "\t" + str(shortest_path_len) + "\t" + str(shortest_paths_count) + "\t" + str(normalized_score) + "\t" + pair[2])
        if not nx.has_path(G, pair[0], pair[1]):
             calculated_scores.append(pair[0] + "\t" + pair[1] + "\t" + "0" + "\t" + "0" + "\t" + "0" + "\t" + pair[2])

# Writing the calculated results into a TVIS_result_"name of the drug".txt file for each drug:
    print("Step 5: Target-vulnerability interaction score (TVIS) output file generation")
    output_file = open(str(argument3) + 'TVIS_results_{0}.txt'.format(listed_drug), "w")
    output_file.write("Target_gene" + "\t" + "Vulnerability_gene" + "\t" + "Length_of_shortest_path" + "\t" + "Number_of_shortest_paths" + "\t" + "TVIS_score" + "\t" + "DGIdb_interaction_score" + "\n")
    for line in calculated_scores:
        output_file.write(line + "\n")
    output_file.close()

# Parsing data needed for percentile value calculation of drug-gene pair using data in TVIS_result_"name of the drug".txt file:
    target_vulnerability_pairs = []
    for target_gene in target_genes:
        for vulnerability in genes:
            target_vulnerability_pairs.append(target_gene.split(",")[0] + "\t" + vulnerability + "\t" + str(target_gene.split(",")[1]))

    TVIR_scores = []
    for target_vulnerability_pair in target_vulnerability_pairs:
        target = target_vulnerability_pair.split("\t")[0]
        gene = target_vulnerability_pair.split("\t")[1]
        int_score = float(target_vulnerability_pair.split("\t")[2])
        input_file = open(str(argument3) + 'TVIS_results_{0}.txt'.format(listed_drug), "r")
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

# Percentile (Target-vulnerability interaction rank (TVIR)) calculation:
        print("Step 6: Target-vulnerability interaction rank (TVIR) calculation")
        ls_genes = []
        for e in data.values():
            if e <= data[gene]:
                ls_genes.append(e)
        percentile = (len(ls_genes) / len(data.values()))
        TVIR_scores.append(target + "\t" + gene + "\t" + str(percentile) + "\t" + str(percentile*int_score))

# Drug-vulnerability correlation score (DVCS)) calculation:
    print("Step 7: Drug-vulnerability correlation score (DVCS) calculation")
    DVCS_score_Data = []
    new_DVCS_score_Data = []
    for l in TVIR_scores:
        DVCS_score_Data.append(float(l.split("\t")[2]))
    for f in TVIR_scores:
        new_DVCS_score_Data.append(float(f.split("\t")[3]))
    summ_DVCS = sum(DVCS_score_Data)
    summ_new_DVCS = sum(new_DVCS_score_Data)
    normalized_DVCS = summ_DVCS / len(target_genes)
    average_DVCS = summ_DVCS / len(target_vulnerability_pairs)
    median_DVCS = median(DVCS_score_Data)

# Writing final results (TVIR and DVCS Scores) into "FINAL_result_"name of the drug".txt file:
    print("Step 8: Saving final results")
    output_file_2 = open(str(argument3) + 'FINAL_result_{0}.txt'.format(listed_drug), "w")
    output_file_2.write("Calculated drug-vulnerability interaction score (DVCS) for " + str(listed_drug) + ":" + "\n" + "\n")
    output_file_2.write("Sum of DVCS Scores: " + str(summ_DVCS) + "\n" + "Sum of normalized DVCS Scores: " + str(summ_new_DVCS) + "\n" + "Target normalized DVCS Score: " + str(normalized_DVCS) + "\n" + "Average DVCS Score: " + str(average_DVCS) + "\n" + "Median DVCS Score: " + str(median_DVCS) + "\n" + "\n" + "Calculated TVIR scores for " + str(listed_drug) + ":" + "\n" + "\n")
    output_file_2.write("Target_gene" + "\t" + "Vulnerability_gene" + "\t" + "TVIR_Score" + "\t" + "DGIdb_normalized_Score" + "\n")
    for data_point in TVIR_scores:
        output_file_2.write(data_point + "\n")
    output_file_2.close()

# Identifying drug targets and random vulnerabilities:
    print("Step 9: Identifying drug targets and random gene vulnerabilities and calculation of DTCS scores using random set of vulnerabilities for " + listed_drug)
    final_data = []
    counter = 1
    while counter < 501:
        print("Number of rounds: " + str(counter) + "/500")
        random_genes = random.sample(NCBI_gene_list, len(genes))

        target_vulnerability_pairs_2 = []
        for target_gene in target_genes:
            for vulnerability in random_genes:
                target_vulnerability_pairs_2.append(target_gene.split(",")[0] + "\t" + vulnerability + "\t" + str(target_gene.split(",")[1]))

        TVIR_scores = []
        for target_vulnerability_pair in target_vulnerability_pairs_2:
            target = target_vulnerability_pair.split("\t")[0]
            gene = target_vulnerability_pair.split("\t")[1]
            int_score = float(target_vulnerability_pair.split("\t")[2])
            input_file = open("TVIS_results_{0}.txt".format(listed_drug), "r")
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

# Percentile (Target-vulnerability interaction rank (TVIR)) calculation:
            ls_genes = []
            for e in data.values():
                if e <= data[gene]:
                    ls_genes.append(e)
            percentile = (len(ls_genes) / len(data.values()))
            TVIR_scores.append(target + "\t" + gene + "\t" + str(percentile) + "\t" + str(percentile * int_score))
            
        # Drug-vulnerability correlation score (DVCS)) calculation:
        DVCS_score_Data = []
        new_DVCS_score_Data = []
        for l in TVIR_scores:
            DVCS_score_Data.append(float(l.split("\t")[2]))
        for f in TVIR_scores:
            new_DVCS_score_Data.append(float(f.split("\t")[3]))
        summ_DVCS = sum(DVCS_score_Data)
        summ_new_DVCS = sum(new_DVCS_score_Data)
        normalized_DVCS = summ_DVCS / len(target_genes)
        average_DVCS = summ_DVCS / len(target_vulnerability_pairs)
        median_DVCS = median(DVCS_score_Data)
        final_data.append(str(listed_drug) + "\t" + str(counter) + "\t" + str(summ_DVCS) + "\t" + str(summ_new_DVCS)+ "\t" + str(normalized_DVCS) + "\t" + str(median_DVCS))
        counter += 1

# Comparing the results of the random values with the actual values:
    print("Step 3: Calculating Drug-Tumor Correlation Rank (DTCR) values for " + listed_drug)
    summ_DTCS_list = []
    summ_new_DTCS_list = []
    normalized_DTCS_list = []
    median_DTCS_list = []
    for line in final_data:
        summ_DTCS_list.append(float(line.split("\t")[2]))
        summ_new_DTCS_list.append(float(line.split("\t")[3]))
        normalized_DTCS_list.append(float(line.split("\t")[4]))
        median_DTCS_list.append(float(line.split("\t")[5]))
    real_data_file = open(str(argument3) + "FINAL_result_{0}.txt".format(listed_drug), "r")
    real_data_lines = real_data_file.read().splitlines()
    real_data_file.close()
    real_sum_DTCS = float(real_data_lines[2].split(": ")[1])
    real_sum_new_DTCS = float(real_data_lines[3].split(": ")[1])
    real_normalized_DTCS = float(real_data_lines[4].split(": ")[1])
    real_median_DTCS = float(real_data_lines[6].split(": ")[1])
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
    for e in normalized_DTCS_list:
        if e <= real_normalized_DTCS:
            normalized_DTCS_genes.append(e)
    normalized_DTCS_percentile = (len(normalized_DTCS_genes) / len(normalized_DTCS_list))
    median_DTCS_genes = []
    for e in median_DTCS_list:
        if e <= real_median_DTCS:
            median_DTCS_genes.append(e)
    median_DTCS_percentile = (len(median_DTCS_genes) / len(median_DTCS_list))

# Writing final results (TVIR and DVCS Scores) into "FINAL_result_"name of the drug".txt file:
    print("Step 4: Saving results")
    output_file_2 = open(str(argument3) + 'DTCR_results_{0}.txt'.format(listed_drug), "w")
    output_file_2.write("DTCR Score: " + str(summ_DTCS_percentile) + "\n" + "\n")
    output_file_2.write("drug" + "\t" + "#_of_rounds" + "\t" + "DTCS_score" + "\t" + "normalized_DTCS_score" +"\t" + "median_DTCS" + "\t" + "normalized_median_DTCS" + "\t")
    for l in final_data:
        output_file_2.write(l + "\n")
    output_file_2.close()

end = time.time()
print(end - start)
