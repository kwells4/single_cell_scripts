import sys
import os
import re

data_file_directory = sys.argv[1]

output_file_name = "HTSeq_gene_file.txt"
file_count = 0

# read in files from the command line
INPUT_FILE_NAME = r'KW-\d+_HTSeq\.txt'  # regex for file name

with open(output_file_name, "w") as output_file:

	for file_name in os.listdir(data_file_directory):

		# only read files with the format of the sorted bam file
		if re.match(INPUT_FILE_NAME, file_name):
			print(file_name)
			file_count += 1
			
			# rename the complete filepath to include the name of the file
			file_path = os.path.join(data_file_directory, file_name)

			with open(file_path, "r") as HTSeq_file:
				if file_count > 1:
					next(HTSeq_file)
				for line in HTSeq_file:
					output_file.write(line)
