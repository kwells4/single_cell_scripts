'''
Kristen Wells

Determines the location of transcription factor motifs and random motifs in 
the promoters of genes of interest.

Promoter sequences are downloaded from ensemble using a flexible distance up
and downstream of genes of interest. Transcription factor motifs from JASPAR
are then searched for within the promoter sequences. These motifs and positions
are returned in an output file. A flexible number of random motifs of a flexible
length are then generated and their positions are also found within the promoter
sequences. These motifs and their positions are returned in a second output
file.

This takes a file including the list of genes, a file with the list of gene 
ranges cooresponding to the gene list (there can be extra genes in this file),
a file containing all motifs of interest in the psm format from JASPAR, and two
output files, one for transcription factor positions and one for random
motif positions.
'''

import requests
import sys
import re
import time
import gzip
from Bio import motifs
from Bio.Seq import Seq
import copy
import random

# Set background nucleotide percentages
BACKGROUND = {'A': 0.25, 'T': 0.25, 'C': 0.25, 'G': 0.25}

# read in files from the command line
(cluster_file, gene_range_file, motif_file,
	transcription_factor_positions, random_motif_positions) = sys.argv[1:]

# connect to ensembl to download promoter sequences
server = 'http://rest.ensembl.org'

# define paramaters for ensemble. If working with non-mouse, change 
# mus_musculus to your species name
ext = '/sequence/region/mus_musculus'

# define headers for calling promoters from ensembl
headers = {"Content-Type": "application/json", "Accept": "application/json"}

# define headers for the transcription factor output file
TF_HEADER = ['Chromosome', 'Gene_ID', 'TF', 'TF_Seq', 'TF_Start', 'TF_End']
MOTIF_HEADER = ['Chromosome', 'Gene_ID', 'MOTIF_SEQ', 'TF_Start', 'TF_End']

# define empty lists and dictionaries
gene_list = []
data_dict = {}
data_list = []
sequence_list = []
position_dict = {}

# define position of chromosome, gene, motif start, and motif end
MOTIF_CHROMOSOME = 0
GENE = 0
MOTIF_START_POS = 1
MOTIF_END_POS = 1

PSEUDOCOUNTS = 0.5  # determines the value to normalize the JASPAR matrix
PRECISION = 10**4  # determines the precision for the distribution
BALANCE = 1000  # determines the value use to balance FDR and FNR

RANGE_START = 0  # the start of the range for number of random motifs generated
RANGE_END = 1000  # the end of the range for number of random motifs generated
RANDOM_MOTIF_LENGTH = 7  # length of random motifs created

SLEEP_TIME = 2  # determines the amount of time to wait between fetching sequences

input_5_prime = 1000  # number of basepairs upstream of gene start
input_3_prime = 100  # number of basepairs downstream of gene start


def find_promoters(promoter_chromosome, start_position, end_position, gene_id):
	'''
	connects to the ensembl server and downloads the sequence for a region of
	interest. Takes chromosome, start position, end position, and gene id as
	arguments. Passes these arguments along with the pre-set upstream and 
	downstream parameters to ensembl and returns the sequence of the given
	region.
	'''

	# defines the list of information to be passed to ensembl and puts
	# it into a format readible by ensembl
	my_line = ('{ "regions" : ["' + promoter_chromosome + ':' +
				str(start_position) + '..' + str(end_position) + ':1"]}')

	# request sequence for the given position from ensembl
	r = requests.post(server+ext, headers=headers,
		data=my_line)
	if not r.ok:
		r.raise_for_status()
		sys.exit()

	# sends the request to ensembl, the return is a list of dictionaries
	sequence_info = r.json()

	# for each dictionary in the list of dictionaries returned, add
	# the gene name, start position, end position, and chromsome to the
	# dictionary
	for item in sequence_info:
		item['gene'] = gene_id
		item['start_position'] = start_position
		item['end_position'] = end_position
		item['chromosome'] = promoter_chromosome
	
		# append the entre dictionary to the list of all sequence
		# dictionaries
		sequence_list.append(item)

	# wait for two seconds before calling more items from ensembl
	time.sleep(SLEEP_TIME)

	# return the list of all sequence dictionaries
	return sequence_list


def motif_search(motif_matrix, promoter_sequence_list, position_dictionary):
	'''
	takes a JASPAR file formatted matrix and a sequence and determines the
	location of the given motif in each sequence provided. The positions 
	in the motif matrix are normalized and a threshold is determined based
	on false negative and false positive rates. Only positions where the
	threshold is met are returned in a dictionary of the gene, the position
	of the transcription factor within the promoter, and the chromosome of
	interest. This function takes a motif matrix, a list of sequences, and
	an empty position dictionary as arguments
	'''
	
	# determine the length of the motif. This will be used to determine
	# the end position of the motif in future calculations
	motif_length = len(motif.degenerate_consensus)
	
	# normalize the counts using a position weight matrix a pseudocount
	# avoids overfitting of the position-weight matrix to the limited number
	# of instances in the alignment
	pwm = motif.counts.normalize(pseudocounts=PSEUDOCOUNTS)

	# determine the log odds ration for each position in the motif
	# determines the log odds of a nucleotide coming from the motif vs
	# the background
	pssm = pwm.log_odds()

	# determine the distribution of the log odds score
	distribution = pssm.distribution(background=BACKGROUND,
							precision=PRECISION)

	# set a threshold based on both false negative and false positive rate
	threshold = distribution.threshold_balanced(BALANCE)
	
	# read through each item in the promoter sequences
	for item in promoter_sequence_list:

		# define the promoter sequence as a sequence
		promoter_sequence = Seq(item['seq'], motif.alphabet)

		# define the location of gene name, start position, end position,
		# and chromosome
		gene_name = item['gene']
		promoter_start_position = item['start_position']
		promoter_end_position = item['end_position']
		gene_chromosome = item['chromosome']

		# for each position in the sequence, search for the tf motif of interest
		# and only return the position if the pssm score is above the threshold
		for position, score in pssm.search(promoter_sequence,
			threshold=threshold):

			# define the motif position start site using the 
			# sequence start, and the position of the tf motif
			# within the sequence. Define the motif position end 
			# site by using the motif start site and the length
			# of the motif
			motif_position_start = promoter_start_position + position
			motif_position_end = motif_position_start + motif_length

			# add inormation to the position dictionary with the gene
			# ensembl id as the key and the chromosome, start, and stop
			# positions as the values
			position_dictionary[gene_name] = [gene_chromosome,
				motif_position_start, motif_position_end]

	# return the dictionary of positions
	return position_dictionary

						
def random_motif_search(random_motif, promoter_sequence_list,
	random_position_dictionary):
	'''
	Takes a random motif generated and a sequence and determines the
	location of the given motif in each sequence provided. Locations where
	the motif is found are returned in a dictionary of the gene, the position
	of the motif within the promoter, and the chromosome of
	interest. This function takes a motif, a list of sequences, and
	an empty position dictionary as arguments
	'''

	# determine the motif length
	motif_length = len(random_motif)

	# compile the motif to make a motif patter
	motif_pattern = re.compile(random_motif)

	# iterate through the promoter sequences created previously
	for item in promoter_sequence_list:

		# define the promoter sequence from the list
		promoter_sequence = item['seq']

		# define the location of gene name, start position, end position,
		# and chromosome
		gene_name = item['gene']
		promoter_start_position = item['start_position']
		promoter_end_position = item['end_position']
		gene_chromosome = item['chromosome']

		# find all occurances of an exact match of the motif within the 
		# promoter. Return the start and end positions of the promoter
		for sequence in motif_pattern.finditer(promoter_sequence):
			position = sequence.start()
			motif_position_start = promoter_start_position + position
			motif_position_end = motif_position_start + motif_length

			# add inormation to the position dictionary with the gene
			# ensembl id as the key and the chromosome, start, and stop
			# positions as the values
			random_position_dictionary[gene_name,
				motif_position_start] = [gene_chromosome,
				motif_position_end]

	# return the created dictionary
	return random_position_dictionary

# open the cluster file of interest and read line by line
with open(cluster_file, 'r') as transcript_file:
	for line in transcript_file:

		# only work with the line if it is associated with an Ensembl id
		if re.search(r'ENS', line):

			# define the names of each column in the file
			transcript_number, ensemble_id_cluster = line.rstrip().split(',')
			
			# make the ensembl id into a string that doesn't include ""
			new_ensemble_id = re.sub(r'\"', '', ensemble_id_cluster)
			new_ensemble_id_str = ''.join(new_ensemble_id)
			
			# make a list of all ensembl ids
			gene_list.append(new_ensemble_id_str)

# open the gene range file and read line by line
with open(gene_range_file, 'r') as range_file:
	for line in range_file:

		# only read the line if it is associated with an Ensembl id
		if re.search(r'ENS', line):

			# name each column in the line
			(range_transcript_number, start_pos, end_pos, width,
				ensemble_id, strand, chromosome) = line.rstrip().split(',')
			
			# read each gene from the list created from the cluster file
			for cluster_gene in gene_list:

				# if the gene from the cluster file is the same as the gene
				# in the ranges file, continue working with the line.
				# if no genes in the cluster file match the line in the ranges
				# file, move onto the next line
				if ensemble_id == cluster_gene:

					# print the id to keep track of the progress through
					# the list of genes
					print(ensemble_id)

					# define the gene start and end as integers
					gene_start = int(start_pos)
					gene_end = int(end_pos)
					gene_chromosome = str(chromosome)

					# determine the start and end sites of the promoter
					# sequence using user define 5 prime and 3 prime positions
					# and the start and end position of the chromosome
					five_prime_start = gene_start - input_5_prime
					three_prime_start = gene_start + input_3_prime
					three_prime_end = gene_end + input_3_prime

					# call the find_promoter function that sends position
					# information to ensembl and returns a gene sequence
					promoter_sequence_list = find_promoters(gene_chromosome,
						five_prime_start, three_prime_start, cluster_gene)
						
# open the output file
with open(transcription_factor_positions, 'w') as output_file:

	# write the header to the output file
	output_file.write('\t'.join(TF_HEADER) + '\n')

	# open the motif file in jaspar format
	with open(motif_file) as handle:

		# read each motif in the jaspar file individually
		for motif in motifs.parse(handle, 'jaspar'):

			# define the tf_name and the consensus sequence
			tf_name = motif.name

			# print the tf name to keep track of the progress trhough the list
			# of tfs.
			print(tf_name)

			# assign a variavle to the tf sequence
			tf_sequence = str(motif.degenerate_consensus)

			# make a copy of the tf_pos_dict so the same dict can be reused
			tf_position_dict = copy.deepcopy(position_dict)

			# call the motif_search function which takes sequences and searches
			# for each motif within those squences
			tf_location_dictionary = motif_search(motif, promoter_sequence_list,
				tf_position_dict) 

			# read the tf_location_dictionary key by key
			for gene in tf_location_dictionary:

				# provide a variable for the gene id, chromosome, and position
				gene_identification = gene
				gene_chr = str(tf_location_dictionary[gene][MOTIF_CHROMOSOME])
				tf_start_pos = str(tf_location_dictionary[gene]
					[MOTIF_START_POS])
				tf_end_pos = str(tf_location_dictionary[gene][MOTIF_END_POS])

				# define the list of variables to write to a line in the output
				# file
				line_to_write = [gene_chr, gene_identification, tf_name,
					tf_sequence, tf_start_pos, tf_end_pos]

				# write the line to the output file in tab delimited format
				output_file.write('\t'.join(line_to_write) + '\n')

# open the file to write random motif positions to
with open(random_motif_positions, 'w') as random_output_file:

	# write the header to the output file
	random_output_file.write('\t'.join(MOTIF_HEADER) + '\n')

	# create a pre-defined number of random motifs
	for number in range(RANGE_START, RANGE_END):

		# create and empty DNA string
		DNA = ""

		# make DNA strings the pre-determined length and randomly pick between
		# the four nucleotides
		for count in range(RANDOM_MOTIF_LENGTH):
			DNA += random.choice("GCTA")

		# make a copy of the motif_pos_dict so the same dict can be reused
		motif_position_dict = copy.deepcopy(position_dict)

		# call the random motif search function using the motif created, the
		# list of prmoter sequences, and the blank motif position dictionary
		random_motif_location_dictionary = random_motif_search(DNA,
			promoter_sequence_list, motif_position_dict)

		if random_motif_location_dictionary is not None:

			# read the motif_location_dictionary key by key
			for gene in random_motif_location_dictionary:

				# provide a variable for the gene id, chromosome, and position
				gene_identification = gene[GENE]
				gene_chr = str(random_motif_location_dictionary[gene][
					MOTIF_CHROMOSOME])
				motif_start_pos = str(gene[MOTIF_START_POS])
				motif_end_pos = str(random_motif_location_dictionary[gene][
					MOTIF_END_POS])

				# define the list of variables to write to a line in the output
				# file
				line_to_write = [gene_chr, gene_identification, DNA,
					motif_start_pos, motif_end_pos]

				# write the line to the output file in tab delimited format
				random_output_file.write('\t'.join(line_to_write) + '\n')