import sys
import os
import re
from collections import defaultdict

data_file_directory, date = sys.argv[1:]
mito_list = []
gene_dict = defaultdict(list)
sample_list = ["sample"]
file_count = 0
output_file_name = "HTSeq_gene_file_2.txt"

# read in files from the command line
INPUT_FILE_NAME = r'KW-\d+_HTSeq\.txt'  # regex for file name
SAMPLE_NAME = r'KW-\d+'  # regex for sample namec

for file_name in os.listdir(data_file_directory):

	# only read files with the format of the sorted bam file
	if re.match(INPUT_FILE_NAME, file_name):
		print(file_name)
		file_count += 1
		
		# rename the complete filepath to include the name of the file
		file_path = os.path.join(data_file_directory, file_name)

		# # name the input file
		txt_input = file_path

		# name the sample based on the file name
		sample_name = ''.join(re.findall(SAMPLE_NAME, file_name))
		sample_title = sample_name + "_" + date

		with open(file_path, "r") as HTSeq_file:

			line_count = 0
			for line in HTSeq_file:
				line_count += 1
				line = line.strip().split('\t')
				gene_name, count = line[0:]

				if line_count == 1:
					sample_list.append(sample_title)
				gene_dict[gene_name].append(count)
		
				

with open(output_file_name, "w") as output_file:
	output_file.write("\t".join(sample_list) + "\n")
	for gene in gene_dict:
		print(len(gene_dict[gene]))
		output_file.write(gene + "\t" + "\t".join(gene_dict[gene]) + "\n")



		