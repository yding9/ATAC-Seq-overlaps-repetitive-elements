from __future__ import division
import sys

if len(sys.argv) != 4:
	print("Invalid input. Please specify input and output files. (command example: python inputfile outputfile replaceName[1/0])")
	sys.exit(0)

#open the input file
inputfile = open(sys.argv[1])
#open the output file
outputfile = open(sys.argv[2], "w")

if sys.argv[3] == '1':
	replace_name = True
elif sys.argv[3] == '0':
	replace_name = False
else:
	print("Wrong Third Arugment")
	sys.exit(2)
#print header arrays
outputfile.write("chr\tstart\tend\tstrand\tScore\n")

chr_name_replace_dict = {"NC_000001.11":"chr1","NC_000002.12":"chr2","NC_000003.12":"chr3", "NC_000004.12":"chr4", \
"NC_000005.10":"chr5","NC_000006.12":"chr6","NC_000007.14":"chr7","NC_000008.11":"chr8","NC_000009.12":"chr9",\
"NC_000010.11":"chr10","NC_000011.10":"chr11","NC_000012.12":"chr12","NC_000013.11":"chr13","NC_000014.9":"chr14",\
"NC_000015.10":"chr15","NC_000016.10":"chr16","NC_000017.11":"chr17","NC_000018.10":"chr18","NC_000019.10":"chr19",\
"NC_000020.11":"chr20","NC_000021.9":"chr21","NC_000022.11":"chr22","NC_000023.11":"chrX","NC_000024.10":"chrY",\
"NC_012920.1":"chrM"}

for line in inputfile:
	# skip the line start with '#'
	if line[0] == '#':
		continue
	# split the line into an array
	array = line.rstrip("\n").split('\t')
	chr_name = array[1]
	if replace_name == True :
		try:
			chr_name = chr_name_replace_dict[chr_name]
		except:
			print("Key does not exist: ", chr_name)

	outputfile.write(chr_name + '\t' + array[2] + '\t' + array[3] + '\t' + array[4] + '\t' + array[7] + '\n')


#close all the files.	
outputfile.close()
inputfile.close()
