'''
Takes a bam file containing sequence information for a region and determines
the average coverage at each position in regions surrounding a transcription
factor. Requires an indexed bam file, a list of transcription factors and their
positions, a list of random motifs and their positions, a transcription
factor output and a random motif output file as arguments. The list of
transcription factors and positions and the list of random motifs and positions
should be pre-compiled from searching the promoter regions of genes of 
interest.

For each transcription factor site, an area around the transcription factor is
defined (ex. 100). A pileup is generated from the bam file
for each position from 100 (or some other pre-determined distance) basepairs
upstream to 100 basepairs downstream of the transcription factor site. A
dictionary is made for each transcription factor (or random motif) that
contains a position as the key and a list of all pileups at that position
as the value. For each transcription factor (or random motif) the average
at each position is then taken. For the random motifs, the average at each
position for every motif is determined. 

The output is the average pileup at each position for each transcription factor.
A second output is the average pileup at each postion for every random motif
combined
'''

import sys
import pysam
import copy
import re
import math
import numpy

# read in files from the command line
(bam_file, transcription_factors, rand_motifs, pileup_counts,
	rand_pileup_counts) = sys.argv[1:]

# name the bamfile using pysam
bamfile = pysam.AlignmentFile(bam_file, "rb")

#make a header for the output files
HEADER = ['Transcription Factor', 'Genes Containg TF']
RAND_HEADER = []

# determine how large of a region around the tf to look at
REGION_SUR_TF = 100
AVERAGE = 2 # define two to determine the average
ZERO = 0 # define the value necessary to exceed to work with the motif pileups
ONE = 1 # used to adjust the region around tf because it is 0 based

# empty dictionary for pileup positions and counts
coverage_dict = {}
gene_dictionary = {}
coverage_count_dict = {}
rand_coverage_count_dict = {}
average_coverage_dict = {}
total_motif_dict = {}
random_gene_dictionary = {}
percent_list = []

def determine_coverage(pileup_dictionary, tf_dict, gene_dict, tf_mid_point,
	tf_id, gene_id):
	'''
	Determines the average coverage for each around a transcription factor
	motif. The middle of the transcription factor motif is defined as 0. The
	number of nucleotides on either side of the transcription factor is
	determined by a pre_determined value set at the top of this code. A 
	dictionary is compiled with the position as the key and the coverage at
	that position as the value. This dictionary is then put as the value
	for another dictionary that has the transcription factor as the key
	This allows each the coverage surrounding each transcription factor to
	be evaluate. This function takes a pileup dictionary, a transcription factor
	dictionary, a gene dictionary, a transcription factor midpoint, a
	transcription factor, and a gene name.
	'''

	# if the transcription factor is already in the gene dictionary, add the
	# new counts to the existing key.
	if tf_id in tf_dict:

		# for each position in the pileup dictionary, define the position
		# relative to the transcription factor.
		for position in pileup_dictionary:
			tf_relative_pos = position - tf_mid_point

			# if the position is already in the dictionary, append the pileup
			# value to the list of values
			if tf_relative_pos in tf_dict[tf_id]:

				tf_dict[tf_id][tf_relative_pos].append(pileup_dictionary[
					position])

			# if the position is within the pre-defined areas around the 
			# transcription factor and the position is not already in the 
			# dictionary, add the position as the key and the pileup value
			# as a list as the value.
			if (-REGION_SUR_TF <= tf_relative_pos <= REGION_SUR_TF and
				tf_relative_pos not in tf_dict[tf_id]):
				tf_dict[tf_id][tf_relative_pos] = [pileup_dictionary[position]]

	# if the transcription factor is not already in the gene dictionary create
	# add the transcription factor as the key and an empty dictionary as the
	# value in the transcription factor dictionary.
	if tf_id not in tf_dict:
		tf_dict[tf_id] = {}

		# for each position in the pileup dictionary, redefine the position
		# relative to the midpoint of the transcription factor.
		for position in pileup_dictionary:
			tf_relative_pos = position - tf_mid_point

			# if the position is within th pre-determined region relative
			# to the midpoint of the transcription factor, add the position as
			# the key and the pileup value as a list as the value.
			if -REGION_SUR_TF <= tf_relative_pos <= REGION_SUR_TF:
				tf_dict[tf_id][tf_relative_pos] = [pileup_dictionary[position]]
			
	# if the transcription factor is not already in the gene dictionary,
	# add the transcription factor as key key and the gene as the start
	# of a list as the value
	if tf_id not in gene_dict:
		gene_dict[tf_id] = [gene_id]

	# if the transcription factor is already in the gene dictionary, append
	# the gene to the list of genes that are the value for the transcription
	# factor key
	if tf_id in gene_dict:
		gene_dict[tf_id].append(gene_id)

	# return both the transcription factor dictionary and the gene dictionary
	return(tf_dict, gene_dict)

def determine_total_motif_coverage(pileup_pos_dict, complete_pileup_dict):
	'''
	This function takes the coverage for each motif and combines it into
	one list containing the coverage for each position. This function is only
	used with the random motifs. It takes a pileup pos dictionary (that is 
	organized with the position as the key and a list of the pileups as the
	values) 
	'''

	# look at every position in the pileup dictionary
	for position in pileup_pos_dict:

		# if the position is already in the complete dictionary, extend
		# the list of pileups.
		if position in complete_pileup_dict:
			complete_pileup_dict[position].extend(pileup_pos_dict[position])

		# if the position is not in the complete dictionary and the position
		# is within the pre-defined region, make the position a key in the 
		# dictionary and the list of pileups as the value.
		if (-REGION_SUR_TF <= position <= REGION_SUR_TF and position not in
			complete_pileup_dict:
			complete_pileup_dict[position] = pileup_pos_dict[position]

	# return the complete dictionary
	return complete_pileup_dict

def average_coverage(total_coverage_dict):
	'''
	Determines the average coverage for each position relative to a motif.
	Takes a coverage dictionary with positions as keys and lists of pileups
	as values.
	'''

	# create an empty averge coverage dictionary
	average_coverage_dict = {}

	# look at each position in the total coverage dictionar
	for position in total_coverage_dict:
		
		# for each position, take the average of every pileup value in the list
		average_coverage = numpy.average(total_coverage_dict[position])

		# add the position as the key and the average coverage as the value to
		# the average coverage dictionary
		average_coverage_dict[position] = average_coverage 

	# return the average coverage dictionary
	return average_coverage_dict



# open the transcription factor position list file
with open(transcription_factors, 'r') as region_file:

	# read the file line by line
	for line in region_file:

		# for each line, make a new empty copy of the coverage dict
		tf_coverage = copy.deepcopy(coverage_dict)
		
		# only read the line if there a position values in the line
		if re.search(r'\d', line):

			# provide variable names for each column in the file
			(gene_chr, gene_identification, tf_name, tf_sequence, tf_start,
				tf_end) = line.strip().split('\t')

			# make the tf start and stop into integers and find the middle
			# of the tf
			tf_start_int = int(tf_start)
			tf_end_int = int(tf_end)
			mid_point = math.floor((tf_start_int + tf_end_int) / AVERAGE)

			# determine the region to look for counts using a pre-defined
			# number of nucluetides surrounding the motif
			start = mid_point - REGION_SUR_TF
			stop = mid_point + REGION_SUR_TF

			# run pysam pileup on the bamfile using the chromosome from the
			# file and the start and stop positions defined based on the
			# tf motif position
			for pileup_column in bamfile.pileup(gene_chr, start, stop):
				pileup_position = pileup_column.pos
				pileup_count = pileup_column.n

				# add each pileup position and count to the tf_coverage dic
				tf_coverage[pileup_position] = pileup_count

			# call the determine coverage function using the tf coverage
			# dict, coverage count dict, gene dict, mid point, tf and gene
			total_coverage = determine_coverage(tf_coverage,
				coverage_count_dict, gene_dictionary, mid_point, tf_name,
				gene_identification)

			# if the return of the function is a value provide names for
			# each part of the return value
			if total_coverage is not None:
				coverage_count_dict, gene_dictionary = total_coverage

# open the transcription factor position list file
with open(rand_motifs, 'r') as rand_region_file:
	
	# read the file line by line
	for line in rand_region_file:

		# for each line, make a new empty copy of the coverage dict
		motif_coverage = copy.deepcopy(coverage_dict)
		
		# only read the line if there a position values in the line
		if re.search(r'\d', line):

			# provide variable names for each column in the file
			(gene_chr, gene_identification, motif_sequence, motif_start,
				motif_end) = line.strip().split('\t')

			# make the tf start and stop into integers and define the midpoint
			# of the motif
			motif_start_int = int(motif_start)
			motif_end_int = int(motif_end)
			motif_mid_point = math.floor((motif_start_int + motif_end_int) / 
				AVERAGE)

		
			# determine the region to look for counts using a pre-defined
			# number of nucluetides surrounding the motif
			pileup_start = motif_mid_point - REGION_SUR_TF
			pileup_stop = motif_mid_point + REGION_SUR_TF
			


			# run pysam pileup on the bamfile using the chromosome from the
			# file and the start and stop positions defined based on the
			# tf motif position
			for pileup_column in bamfile.pileup(gene_chr, pileup_start,
				pileup_stop):
				pileup_position = pileup_column.pos
				pileup_count = pileup_column.n


				# add each pileup position and count to the motif_coverage dic
				motif_coverage[pileup_position] = pileup_count

			# if the motif_coverage dictionary contains pileup information
			# continue working with the motif
			if len(motif_coverage.keys()) > ZERO:

				# call the determine coverage function using the motif coverage
				# dict, the random coverage count dict, the random gene dict,
				# the motif mid point, the motif sequence, and the gene id
				rand_total_coverage = determine_coverage(motif_coverage,
					rand_coverage_count_dict, random_gene_dictionary,
					motif_mid_point, motif_sequence, gene_identification)

				# if the determine coverage function retuns a value, provide
				# variable names for each part of the function
				if rand_total_coverage is not None:

					(rand_coverage_count_dict,
						rand_gene_dictionary) = rand_total_coverage
					
# for each motif in the previously made random coverage count dictionary
for motif in rand_coverage_count_dict:

	# call the determin total motif coverage function which returns a list
	# of every pileup at each position (rather than each pileup being split
	# into each motif)
	total_motif_coverage = determine_total_motif_coverage(
		rand_coverage_count_dict[motif], total_motif_dict)

# open up the pileup counts file to write to
with open(pileup_counts, 'w') as pileup_file:

	# make a header for the pileup counts that includes pre-determined info
	# plus a number for each region of interest (ie -100 to 100).
	for number in range(-REGION_SUR_TF, REGION_SUR_TF + ONE):
		HEADER.append(str(number))


	# write the tab delimited header to the pileup file
	pileup_file.write('\t'.join(HEADER) + '\n')

	# cycle through all the transcription factors in the coverage count
	# dictionary
	for tf in coverage_count_dict:

		# call the average coverage function to determine the average
		# pileup for each position
		tf_average_coverage = average_coverage(coverage_count_dict[tf])

		# determine the number of genes containing that tf motif in the 
		# promoter region
		number_of_genes_containing_tf = str(len(gene_dictionary[tf]))

		# define the list to print
		print_list = [tf, number_of_genes_containing_tf]

		# if a value is returned from the determine coverage function
		if tf_average_coverage is not None:

			# define an empty average coverage list
			average_coverage_list = []

			# iterate through every number in the region surrounding the 
			# transcription factor
			for number in range (-REGION_SUR_TF, REGION_SUR_TF + ONE):

				# if the number is already in the average coverage list, append
				# the average pileup value to the list
				if number in tf_average_coverage:
					average_coverage_list.append(str(
						tf_average_coverage[number]))

				# if the number is not already in the average coverage list,
				# append 0 in the place of the number in the list
				if number not in tf_average_coverage:
					average_coverage_list.append('0')

		# add all values from the average coverage to the print list
		print_list.extend(average_coverage_list)


		# write the line to the output file
		pileup_file.write('\t'.join(print_list) +'\n')


# open up the pileup counts file to write to
with open(rand_pileup_counts, 'w') as rand_pileup_file:

	# make a header for the pileup counts that includes a number for each
	# region of interest (ie -100 to 100).
	for number in range(-REGION_SUR_TF, REGION_SUR_TF + ONE):
		RAND_HEADER.append(str(number))


		# write the tab delimited header to the pileup file
		rand_pileup_file.write('\t'.join(RAND_HEADER) + '\n')

		# call the average coverage function that will take an average of
		# every pileup for every position for all randomly generated motifs
		rand_average_coverage = average_coverage(total_motif_coverage)

		# if a value is returned from the determine coverage function
		if rand_average_coverage is not None:
			
			# create an empty coverage list
			rand_average_coverage_list = []

			# iterate through every number in the region surrounding the 
			# motif
			for number in range (-REGION_SUR_TF, REGION_SUR_TF + ONE):

				# if the number is in the random average coverage add that
				# number to the average coverage list
				if number in rand_average_coverage:
					rand_average_coverage_list.append(str(
						rand_average_coverage[number]))

				# if the number is not in random average coverage, add 0 
				# in place of that position
				if number not in rand_average_coverage:
					rand_average_coverage_list.append('0')

			# write the average coverage list to the output file
			rand_pileup_file.write('\t'.join(rand_average_coverage_list) +'\n')
