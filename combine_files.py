import sys

file_one, file_two, output_file = sys.argv[1:]
HEADER = True

with open(output_file, "w") as outfile:
	with open(file_one, "r") as first_file:
		for line in first_file:
			outfile.write(line)

	with open(file_two, "r") as second_file:
		if HEADER:
			second_file.readline()
		for line in second_file:
			outfile.write(line)