import sys

gene_file = sys.argv[1]
count = 0

with open(gene_file, "r") as input_file:
	for line in input_file:
		line = line.strip().split("\t")
		count += 1
		if count == 1:
			cell_num = len(line)
		line_len = len(line)
		if cell_num != line_len:
			print("false")

		if "ERCC" in line[0]:
			print("true")

		if "_" in line[0]:
			print("sample?")