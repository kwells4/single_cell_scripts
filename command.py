import subprocess
import os
import sys

rootdir = sys.argv[1]


for subdir, dirs, files in os.walk(rootdir):
	for file in files:
		file_dir = os.path.join(subdir, file)
		print(file_dir)
		subprocess.call(["fastqc", file_dir])