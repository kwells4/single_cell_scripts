import sys
import os
import re
from collections import defaultdict

gtf_file = sys.argv[1]
gene_dict = {}
gene_dict_two = {}
transcript_output_file = "ensemble_to_gene_name.txt"
weird_genes_output = "genes_output.txt"
weird_gene_dict = {}


with open(gtf_file, "r") as transcript_file:
	
	for line in transcript_file:
		if "gene_status" in line:
			gene_known = False
			line_dict = {}
			line = line.rstrip().split('\t')
			gene_info = line[8].split(";")
			
			for descriptor in gene_info:
				descriptor = descriptor.strip().split(" ")
				if len(descriptor) > 1:
					line_dict[descriptor[0]]=descriptor[1]
			
			gene_id = line_dict['gene_id'][1:-1]

			if line_dict["gene_status"] == '"KNOWN"':
				gene_known = True
			
			if gene_known:
				gene_name = line_dict['gene_name'][1:-1]
				if gene_name in gene_dict_two and gene_id != gene_dict_two[gene_name]:
					weird_gene_dict[gene_name] = [gene_id, gene_dict_two[gene_name]]
					
				gene_dict[gene_id] = gene_name
				gene_dict_two[gene_name] = gene_id
				

print(len(gene_dict))
print(len(gene_dict_two))

with open(weird_genes_output, "w") as output_file:
	for gene in weird_gene_dict:
		output_file.write(gene + "\t" + "\t".join(weird_gene_dict[gene]) + "\n")

# with open(transcript_output_file, "w") as output_file:
# 	for gene in gene_dict:
# 		output_file.write(gene + "\t" + gene_dict[gene] + "\n")