import sys

gene_file, count_file = sys.argv[1:]
ERCC_list = []

with open(gene_file, "r") as input_file:
	for line in input_file:
		if "ERCC" in line:
			line = line.strip()
			ERCC_list.append(line)

with open(count_file, "r") as input_file:
	for line in input_file:
		line = line.strip().split("\t")
		print(line)
		if line[0] in ERCC_list:
			print(line)