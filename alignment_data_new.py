import sys
import os
import re
from collections import defaultdict

data_file_directory, output_file, date = sys.argv[1:]

# read in files from the command line
INPUT_FILE_NAME = r'KW-\d+Log\.final\.out'  # regex for file name
SAMPLE_NAME = r'KW-\d+'  # regex for sample namec
file_number = 0
header_list = ['sample_ID']
values_dict = defaultdict(list)

# read each file from the directory
for file_name in os.listdir(data_file_directory):

	# only read files with the format of the sorted bam file
	if re.match(INPUT_FILE_NAME, file_name):
		file_number += 1

		# rename the complete filepath to include the name of the file
		file_path = os.path.join(data_file_directory, file_name)

		# name the sample based on the file name
		sample_name = ''.join(re.findall(SAMPLE_NAME, file_name)) + "_" + date

		#print(sample_name)

		# name the input file
		log_final = file_path

		with open(log_final, "r") as data_file:
			for line in data_file:
				line = line.strip().split("\t")

				if len(line) > 1:
					if file_number is 1:
						header_list.append(line[0][:-2])
				
					values_dict[sample_name].append(line[1])


with open(output_file, "w") as stats_file:
	stats_file.write('\t'.join(header_list) + '\n')
	for item in values_dict:
		stats_file.write(item + '\t' + '\t'.join(values_dict[item]) + '\n')