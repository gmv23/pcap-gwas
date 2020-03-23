#!/usr/bin/py

#This script takes a vcf file and returns a matrix with the total read depth at each genotype
#at each genotype

#First row is sample names
#First column is chromosome
#Next column is SNP

import sys, os

vcf_path, out_path  = sys.argv[1:]

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

			depth = geno_field.strip().split(":")[2] # Total reads is third entry
			out.write(depth + "\t")	
		line_count += 1
fh.close()
out.close()
