import sys
from collections import defaultdict

gtf_file, gene_file = sys.argv[1:]
gene_dict = {}
gene_dict_two = {}
transcript_output_file = "ensemble_to_gene_name.txt"
weird_genes_output = "genes_output.txt"
gene_list = []
new_gene_dict = defaultdict(list)

with open(gene_file, "r") as input_file:
	for line in input_file:
		gene_id, es_1, es_2 = line.strip().split("\t")
		gene_list.append(gene_id)

with open(gtf_file, "r") as transcript_file:
	for line in transcript_file:
		for gene in gene_list:
			if gene in line and "gene" in line:
				new_gene_dict[gene].append(line)

for item in new_gene_dict:
	print(item)
	for gene_names in new_gene_dict[item]:
		print(gene_names)
