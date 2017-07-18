import subprocess
import sys
import os
import re

data_file_directory = sys.argv[1]

# read in files from the command line
data_file_directory = sys.argv[1]



INPUT_FILE_NAME = r'^[KW].+[gz]$'  # regex for file name
SAMPLE_NAME = r'KW-\d+'  # regex for sample name

# read each file from the directory
for file_name in os.listdir(data_file_directory):

	# only read files with the format of the sorted bam file
	if re.match(INPUT_FILE_NAME, file_name):
		print(file_name)

		# rename the complete filepath to include the name of the file
		file_path = os.path.join(data_file_directory, file_name)

		# name the sample based on the file name
		sample_name = ''.join(re.findall(SAMPLE_NAME, file_name))
		print(sample_name)