import sys
import os
import re
from collections import defaultdict
		
sample_dict = defaultdict(list)

# read in files from the command line
INPUT_FILE_NAME = r'KW-\d+_HTSeq\.txt'  # regex for file name
SAMPLE_NAME = r'KW-\d+'  # regex for sample namec

output_file_name = "HTSeq_count_file.txt"
header = ["sample_ID", "> 10", ">20", "ERCC >10", "ERCC > 20"]
data_file_directory = sys.argv[1]

# read each file from the directory
for file_name in os.listdir(data_file_directory):

	# only read files with the format of the sorted bam file
	if re.match(INPUT_FILE_NAME, file_name):
		print(file_name)
		ten_count = 0
		twenty_count = 0
		ERCC_ten_count = 0
		ERCC_twenty_count = 0

		# rename the complete filepath to include the name of the file
		file_path = os.path.join(data_file_directory, file_name)

		# # name the input file
		txt_input = file_path

		# name the sample based on the file name
		sample_name = ''.join(re.findall(SAMPLE_NAME, file_name))

		with open(file_path, "r") as HTSeq_file:
			for line in HTSeq_file:
				line = line.strip().split('\t')

				if int(line[1]) >= 10 and "ERCC" not in line[0]:
					ten_count += 1

				if int(line[1]) >= 20 and "ERCC" not in line[0]:
					twenty_count += 1

				if int(line[1]) >= 10 and "ERCC" in line[0]:
					ERCC_ten_count += 1

				if int(line[1]) >= 20 and "ERCC" in line[0]:
					ERCC_twenty_count += 1

		count_list = [str(ten_count), str(twenty_count), str(ERCC_ten_count), str(ERCC_twenty_count)]

		sample_dict[sample_name].extend(count_list)		

with open(output_file_name, 'w') as output_file:
	output_file.write('\t'.join(header) + '\n')
	for item in sample_dict:
		output_file.write(item + '\t' + '\t'.join(sample_dict[item]) + '\n')


