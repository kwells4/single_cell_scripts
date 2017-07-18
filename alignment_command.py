'''
Kristen Wells

Takes a fastq file and aligns to the mouse genome (including ERCCs) using STAR

Then inputs the STAR output bam into HTSeq to count the number of times each gene appears.
'''
import subprocess
import sys
import os
import re

# read in files from the command line
data_file_directory = sys.argv[1]

INPUT_FILE_NAME = r'^[KW].+[gz]$'  # regex for file name
SAMPLE_NAME = r'KW-\d+'  # regex for sample name
gtf_path ="/home/kwells4/star/genome/mm10_ercc92/gencode.vM11.primary_assembly.annotation.ercc92.gtf"


def align_reads(fastq_input, outfile_prefix):
	'''
	Calls STAR aligner
	'''

	star_call = (["STAR", "--runThreadN", "10", "--genomeDir",
		"/home/kwells4/star/genome/mm10_ercc92/index/", " --sjdbGTFfile",
		"/home/kwells4/star/genome/mm10_ercc92/gencode.vM11.primary_assembly.annotation.gtf",
		"--readFilesIn", fastq_input, "--readFilesCommand",
		"zcat", "--outSAMtype", "BAM", "SortedByCoordinate", "--outFileNamePrefix", outfile_prefix])

	subprocess.call(star_call)

def HTseq(bam_input, output_file):
	'''
	Calls HTseq
	'''

	htseq_call = (["htseq-count", 
		"-f", "bam", "-s", "no", bam_input, gtf_path])

	subprocess.call(htseq_call, stdout=output_file)

# read each file from the directory
for file_name in os.listdir(data_file_directory):

	# only read files with the format of the sorted bam file
	if re.match(INPUT_FILE_NAME, file_name):

		# rename the complete filepath to include the name of the file
		file_path = os.path.join(data_file_directory, file_name)

		# name the sample based on the file name
		sample_name = ''.join(re.findall(SAMPLE_NAME, file_name))

		print(sample_name)

		# name the input file
		fastq_file = file_path

		#name star output file
		bam_file = (sample_name + "Aligned.sortedByCoord.out.bam")

		#call the star aligner
		align_reads(fastq_file, sample_name)

		#determine working directory
		working_directory = os.getcwd() + "/"

		#name file to be used for HTseq
		star_output_file = os.path.join(working_directory, bam_file)

		#name output file
		output_name = sample_name + "_HTSeq.txt"

		#open output file for HTseq
		with open(output_name, 'w') as output_file:
				
			#call HTseq
			HTseq(star_output_file, output_file)


		