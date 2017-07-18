import sys
import os
import re
from collections import defaultdict

gtf_file = sys.argv[1]
gene_dict = {}
transcript_output_file = "ensemble_to_gene_name.txt"



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
				gene_dict[gene_id] = gene_name
				

with open(transcript_output_file, "w") as output_file:
	for gene in gene_dict:
		output_file.write(gene + "\t" + gene_dict[gene] + "\n")