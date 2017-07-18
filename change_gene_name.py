import sys

transcript_file, count_file, new_count_file = sys.argv[1:]
gene_dict = {}

with open(transcript_file, "r") as gene_file:
	for line in gene_file:
		gene_id, gene_name = line.strip().split("\t")
		gene_dict[gene_id] = gene_name

with open(new_count_file, "w") as output_file:
	with open(count_file, "r") as input_file:
		for line in input_file:
			line = line.strip().split("\t")
			gene_id = line[0]
			if gene_id in gene_dict:
				line[0] = gene_dict[gene_id]

				output_file.write("\t".join(line) + "\n")
			elif "sample" in line[0]:
				output_file.write("\t".join(line) + "\n")
			elif "ERCC" in line[0]:
				output_file.write("\t".join(line) + "\n")
	
