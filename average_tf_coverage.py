'''
Kristen Wells

Takes a bam file containing sequence information for a region and determines
the average coverage in regions surrounding a transcription factor. Requires
an indexed bam file, a list of transcription factors and their positions, and
an output file as arguments. The list of transcription factors and positions
should be pre-compiled from searching the promoter regions of genes of 
interest.

For each transcription factor site, an area around the transcription factor is
defined (ex. 100). A pileup is generated from the bam file
for each position from 100 (or some other pre-determined distance) basepairs
upstream to 100 basepairs downstream of the transcription factor site. The
pileup is averaged both inside the transcription factor motif and around
the transcription factor motif. This treats every occurance of the trancription
factor motif as independent. Next, an average is taken of every occurance of
each transcription factor. This means that for every transcription factor there
will be an average of around the motif and within the motif that represents
every instance of the transcription factor in the promoters of the genes of 
intrest. The average for outside the motif is then divided by the average within
the motif. If the difference in coverage between the motif and outside the motif
for each transcription factor is large enough, the transcription factor and it's
average coverage is printed in the output file.
'''

import sys
import pysam
import copy
import re

# read in files from the command line
bam_file, transcription_factors, pileup_counts = sys.argv[1:]

# name the bamfile using pysam
bamfile = pysam.AlignmentFile(bam_file, "rb")

# make a header for the output file
HEADER = ['Transcription Factor', 'Average Coverage',
	'Average Motif coverage', 'percent diff in coverage']

# determine how large of a region around the tf to look at
REGION_SUR_TF = 100

ONE = 1  # define 1 to be used as a counter
ZERO = 0  # define 0 to prevent from dividing by zero
DIFF_COVERAGE_BOUNDARY = 0.4  # the minimum percent difference in coverage

# the lowest coverage accepted for a region surrounding the tf
COVERAGE_BOUNDARY = 2
# the lowest coverage accepted for a region surrounding
					   # the tf

# empty dictionary for pileup positions and counts
coverage_dict = {}
gene_dictionary = {}
coverage_count_dict = {}


def determine_coverage(pileup_dictionary, tf_dict, gene_dict, tf_start_pos,
	tf_end_pos, tf_id, gene_id):
	'''
	Determines the average coverage within a motif and around a motif.
	Takes a pileup dictionary, a start position, and an end position for
	the transcription factor.
	'''

	# if a transcription factor is already in the dictionary add new information
	# to the exisiting key in the dictionary
	if tf_id in tf_dict:

		# for every key in the pileup dictionary, determine if the position is
		# within the tf motif
		for position in pileup_dictionary:

			# if the position is in the tf motif, add the coverage value for
			# that position to the total coverage for within the motif and add
			# one to the count for positions in the motif
			if tf_start_pos <= position <= tf_end_pos:
				tf_dict[tf_id]['tf_motif_coverage'] += pileup_dictionary[
					position]
				tf_dict[tf_id]['tf_motif_count'] += ONE

			# if the poisition is not in the tf motif, add the coverage value
			# for that position to the total coverage outside of the motif and
			# add one to the count for positions
			else:
				tf_dict[tf_id]['outside_tf_coverage'] += pileup_dictionary[
					position]
				tf_dict[tf_id]['outside_tf_count'] += ONE

	# if a transcription factor is not in the dictionary, create a key in the
	# dictionary for that transcription factor
	else:

		# set all values for coverage and count to zero
		coverage_counts = {'outside_tf_coverage': 0, 'outside_tf_count': 0,
			'tf_motif_coverage': 0, 'tf_motif_count': 0}

		# for every key in the pileup dictionary, determine if the position is
		# within the tf motif
		for position in pileup_dictionary:

			# if the position is in the tf motif, add the coverage value for
			# that position to the total coverage for within the motif and add
			# one to the count for positions in the motif
			if tf_start_pos <= position <= tf_end_pos:
				coverage_counts['tf_motif_coverage'] += pileup_dictionary[
					position]
				coverage_counts['tf_motif_count'] += ONE

			# if the poisition is not in the tf motif, add the coverage value
			# for that position to the total coverage outside of the motif and
			# add one to the count for positions
			else:
				coverage_counts['outside_tf_coverage'] += pileup_dictionary[
					position]
				coverage_counts['outside_tf_count'] += ONE

		# with the transcription factor as the key, add the coverage counts
		# dictionary as the value
		tf_dict[tf_id] = coverage_counts

	# add transcription factor as the key and the gene as the value to the 
	# gene dictionary
	gene_dict[tf_id] = gene_id

	# return both the transcirption factor dictionary and the gene dictionary
	return(tf_dict, gene_dict)


def average_coverage(total_coverage_dict):
	'''
	Takes the total coverages for within and around a motif and returns the
	average coverage for each transcription factor. Takes a dictionary
	containing the total count and total coverage.
	'''
	
	# only determine average coverage for the tf motif if there is sequence data
	# for the region
	if total_coverage_dict['tf_motif_count'] > ZERO:

		# determine average coverage by dividing total coverage by number of
		# positions
		average_coverage_in_motif = (total_coverage_dict['tf_motif_coverage']
			/ total_coverage_dict['tf_motif_count'])

	# if there is no sequence data for the region, coverage is zero
	if total_coverage_dict['tf_motif_count'] == ZERO:

		average_coverage_in_motif = ZERO

	# only determine average coverage outside of the tf if there is sequence
	# data for the region
	if total_coverage_dict['outside_tf_count'] > ZERO:

		# determine average coverage by dividing total coverage by number of
		# positions
		average_coverage_around_motif = (total_coverage_dict[
			'outside_tf_coverage'] / total_coverage_dict['outside_tf_count'])

	# if there is no sequence data for the region, coverage is zero
	if total_coverage_dict['outside_tf_count'] == ZERO:
		average_coverage_around_motif = ZERO

	# if average coverage surrounding the motif is greater than zero, the 
	# percent difference in coverage is found by dividing the average coverage
	# in the motif by the average coverage around the motif
	if average_coverage_around_motif > ZERO:

		percent_diff_in_coverage = (average_coverage_in_motif / 
			average_coverage_around_motif)

	# if average coverage surrounding the motif is less than zeron, the
	# percent difference in coverage is 0.
	if average_coverage_around_motif == ZERO:
		percent_diff_in_coverage = ZERO

	# if there is a great enough difference between the coverage of the motif
	# and surrounding area and the coverage in the area surrounding the motif
	# is high enough, return the coverages and the difference in coverages
	if (percent_diff_in_coverage < DIFF_COVERAGE_BOUNDARY and
		average_coverage_around_motif > COVERAGE_BOUNDARY):
		return [average_coverage_around_motif, average_coverage_in_motif,
			percent_diff_in_coverage]


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

			# make the tf start and stop into integers
			tf_start_int = int(tf_start)
			tf_end_int = int(tf_end)

			# determine the region to look for counts using a pre-defined
			# number of nucluetides surrounding the motif
			start = tf_start_int - REGION_SUR_TF
			stop = tf_end_int + REGION_SUR_TF
				
			# run pysam pileup on the bamfile using the chromosome from the
			# file and the start and stop positions defined based on the
			# tf motif position
			for pileup_column in bamfile.pileup(gene_chr, start, stop):
				pileup_position = pileup_column.pos
				pileup_count = pileup_column.n

				# add each pileup position and count to the tf_coverage dic
				tf_coverage[pileup_position] = pileup_count

			# call the determine coverage function using the tf coverage
			# dict, tf start, and tf end
			total_coverage = determine_coverage(tf_coverage,
				coverage_count_dict, gene_dictionary, tf_start_int, tf_end_int,
				tf_name, gene_identification)

			# rename the two variables returned from the determine coverage
			# function
			coverage_count_dict, gene_dictionary = total_coverage

# open up the pileup counts file to write to
with open(pileup_counts, 'w') as pileup_file:

	# write the tab delimited header to the pileup file
	pileup_file.write('\t'.join(HEADER) + '\n')

	# read each transcription facotor individually
	for tf in coverage_count_dict:

		# determine the average coverage for each transcription factor
		coverage_difference = average_coverage(coverage_count_dict[tf])

		# if a value is returned from the determine coverage function
		if coverage_difference is not None:

			# provide variable names for each part of the return value
			# from the function
			(coverage_around_motif, coverage_in_motif,
				percent_diff) = coverage_difference

			# Make a gene list of all genes contining the transcription factor
			tf_gene_list = gene_dictionary[tf]

			# define all the variables to write to a line of the output
			# file
			line_to_write = [tf, str(coverage_around_motif),
				str(coverage_in_motif), str(percent_diff)]

			# write the line to the output file
			pileup_file.write('\t'.join(line_to_write) + '\n')
