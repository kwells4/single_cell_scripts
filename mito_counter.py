import sys
import os
import re
from collections import defaultdict

gtf_file, data_file_directory = sys.argv[1:]
mito_list = []
mito_percent_dict = defaultdict(float)

# read in files from the command line
INPUT_FILE_NAME = r'KW-\d+_HTSeq\.txt'  # regex for file name
SAMPLE_NAME = r'KW-\d+'  # regex for sample namec

with open(gtf_file, "r") as transcript_file:
	for line in transcript_file:
		if "chrM" in line:
			line = line.rstrip().split('\t')
			gene_info = line[8].split(";")
			gene_name = gene_info[0]
			gene_name_len = len(gene_name) - 1
			gene_name = gene_name[9:gene_name_len]
			if gene_name not in mito_list:
				mito_list.append(gene_name)

print(mito_list)


for file_name in os.listdir(data_file_directory):

	# only read files with the format of the sorted bam file
	if re.match(INPUT_FILE_NAME, file_name):
		print(file_name)
		
		# rename the complete filepath to include the name of the file
		file_path = os.path.join(data_file_directory, file_name)

		# # name the input file
		txt_input = file_path

		# name the sample based on the file name
		sample_name = ''.join(re.findall(SAMPLE_NAME, file_name))

		with open(file_path, "r") as HTSeq_file:
			mito_count = 0
			total_count = 0

			for line in HTSeq_file:
				line = line.strip().split('\t')
				count = int(line[1])

				if line[0] in mito_list:
					mito_count += count

				elif "ENS" in line[0]:
					total_count += count

			total_aligned = total_count + mito_count
			mito_percent = mito_count / total_aligned * 100
			mito_percent_dict[sample_name] = mito_percent

print(mito_percent_dict)