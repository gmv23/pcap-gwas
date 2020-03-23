#!/usr/bin/py

#This script takes a vcf file and returns a matrix with the ratio of the major allele depth / the total read depth
#at each genotype

#First row is sample names
#First column is chromosome
#Next column is SNP

import sys, os

vcf_path, out_path, min_read_depth  = sys.argv[1:]
min_read_depth = int(min_read_depth)

#Get allele frequencies with vcftools; use to figure out which is major allele

sys.call = "vcftools --vcf %s --freq --out allele_frequencies" %vcf_path
os.system(sys.call)
which_is_major = []

fh = open("allele_frequencies.frq", "r")
fh.readline()

for line in fh:
	line = line.strip().split("\t")
	a1, a2 = [float(x[2:]) for x in line[4:6]]
	if a1 > a2:
		which_is_major.append(0)
	else:
		which_is_major.append(1)

fh.close()

#Begin parsing VCF file

fh = open(vcf_path, "r")
out = open(out_path, "w")

line_count = 0

for line in fh:

	#Create header with chromosome, position, and individual names
	if line[0:6] == "#CHROM":
		out.write("CHROM\tPOS\t" + "\t".join(line.strip().split("\t")[9:]))

	elif line[0] != "#": #Skip other header lines

		out.write("\n")
		line = line.strip().split("\t")

		#write chromosome and position for each marker
		chrom, pos = line[0:2]
		out.write(chrom + "\t" + pos + "\t")

		#loop through genotypes across individuals for marker (ie go across the row of vcf file)
		for geno_field in line[9:]:

			#Genotype is first entry in colon separated values
			genotype = geno_field.strip().split(":")[0]

			#Only calculate ratios for heterozygotes
			if genotype == "1/0" or genotype == "0/1":

				reads = geno_field.strip().split(":")[1] # Number of reads for each allele is second entry
				depths = [int(x) for x in reads.strip().split(",")]
				#Only calculate for read depths greater than min_read_depth
				if sum(depths) > min_read_depth:
					ref_or_alt = which_is_major[line_count]
					major = depths[ref_or_alt] #get major allele read count
					ratio = major/float(sum(depths)) #Divide major allele by total number of reads
				else:
					ratio = "NA" #Write NA if read depth less than 12
				out.write(str(ratio) + "\t")	

			else:
				out.write("NA\t") #Write NA if not heterozygous genotype
		line_count += 1
fh.close()
out.close()
