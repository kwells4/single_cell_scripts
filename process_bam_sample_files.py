'''
Kristen Wells

Takes a sorted bam file and processes it to be used in calling ATAC-seq
transcription factor peaks.

The bam file is first indexed using samtools. Duplicates are then removed from
the sorted and indexed bam file using picard tools MarkDuplicates. Next, 
samtools is again used to filter the reads based on samflag and quality score.
Finally these files are indexed.

The output is an indexed file for the original sorted bam file, a bam file with
no duplicates (that can then be filtered using a different flag if desired), a
filterd bam file and an index for the filtered bam file.
'''

import subprocess
import sys
import os
import re
import time

# read in files from the command line
data_file_directory, my_output_file = sys.argv[1:]


# define all constants
REMOVE_DUPLICATES = 'true'  # make remove duplicates true or flase
ASSUME_SORTED = 'true'  # make assume sorted true or false

# full location of picard tools
MARK_DUPLICATES_LOCATION = '/home/wuwei5/software/picard/picard/MarkDuplicates.jar'

QUALITY = 20  # quality score for processing bam file
FLAG_SORT = 1804  # flag used to process bam file
INPUT_FILE_NAME = r'novo\.sorted\.s\d+_[a-z]+\.bam$'  # regex for file name
SAMPLE_NAME = r's\d+_[a-z]+'  # regex for sample name

# empty dictionary to keep track of completed tasks for each sample
my_output_dict = {}


def remove_duplicate_reads(dup_input_file, dup_output_file, metrics_file):
	'''
	Calls the mark duplicates command from the command line using subprocess
	call. Uses commands for sorted, metrics file, input, output, and location
	specified by variables defined at the top of the program. Will output
	a metrics file and a bam file named based on the output file name. Requires
	a input file, an output file, and a metrics file.
	'''

	# setup the command to remove duplicates using picard tools
	remove_duplicates_call = (["java", "-Xmx4G", "-jar",
		MARK_DUPLICATES_LOCATION, "INPUT=" + bam_input,
		"OUTPUT=" + remove_dup_output, "REMOVE_DUPLICATES=" + REMOVE_DUPLICATES,
		"ASSUME_SORTED=" + ASSUME_SORTED, "METRICS_FILE=" + metrics_file])

	# call the remove duplicates command
	subprocess.call(remove_duplicates_call)

	# return that duplicates were removed
	return 'removed dup'


def filter_for_mapping_quality(map_input_file, map_output_file):
	'''
	Uses samtools view to fileter reads based on quality and samflag. The
	quality score and samflag are specified by variables defined at the top
	of the program. Will output a filtered bam file named based on the output
	file name. Requires an input and an output file.
	'''

	# define the map quality command required to filter reads using samtools
	map_quality = (['samtools', 'view', '-q', str(QUALITY), '-F',
		str(FLAG_SORT), '-b', map_input_file])

	# call the map quality command using the output file as the stdout
	subprocess.call(map_quality, stdout=map_output_file)

	# return that the file was filtered
	return 'map quality'


def index_bam(index_input_file):

	'''
	Uses samtools to index a bam file. Will output an index file based on
	the name of the bam input file. Requires a bam input file.
	'''

	# define the index command required to index a bam file
	index_call = (['samtools', 'index', index_input_file])

	# call the index command
	subprocess.call(index_call)

	# return that the bam file was indexed
	return 'indexed'


# read each file from the directory
for file_name in os.listdir(data_file_directory):

	# only read files with the format of the sorted bam file
	if re.match(INPUT_FILE_NAME, file_name):

		# rename the complete filepath to include the name of the file
		file_path = os.path.join(data_file_directory, file_name)

		# name the sample based on the file name
		sample_name = ''.join(re.findall(SAMPLE_NAME, file_name))

		# name the input file
		bam_input = file_name

		# name the metrics based on the sample name
		metrics_file = ('metrics.' + sample_name + '.txt')

		# name the output of the remove duplicates file based on the sample
		# name
		remove_dup_output = ('novo.sorted.no.dup.' + sample_name + '.bam')

		# name the output of the bam file created from filtering the original
		# bam file with duplicates removed
		sam_flag_output = ('novo.sorted.f1804.' + sample_name + '.bam')

		# call the function to index the original bam file
		index_bam_bai = index_bam(bam_input)

		# call the function to remove duplicates from the original bam file
		remove_duplicates_BAM = remove_duplicate_reads(bam_input,
			remove_dup_output, metrics_file)

		# open the sam file to write to
		with open(sam_flag_output, 'wb') as sam_output_file:

			# call the function to filter the file with duplciates removed
			BAM_map_quality = filter_for_mapping_quality(remove_dup_output, 
				sam_output_file)
		
		# index the filetered bam file
		index_bam_quality_bai = index_bam(sam_flag_output)

		# take the return value from each function called and add it to a list
		# that will keep track of what was completed for each original
		# bam file
		my_output_dict[sample_name] = [index_bam_bai, remove_duplicates_BAM,
			BAM_map_quality, index_bam_quality_bai]

# open an output file to keep track of what was done to each file
with open(my_output_file, 'w') as output_file:

	# write what was done to each original bam file to the output file
	output_file.write(str(my_output_dict) + '\n')