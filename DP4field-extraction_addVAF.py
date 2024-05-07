########################################################
## jens.theine@uni-bielefeld.de / V:07.05.2024
########################################################

################### config #############################

vcf = "calls_q30filteredvariants.vcf"
out = "calls_q30filteredvariants_biSNPs_AN__Chr_Pos_DP_VAF.tsv"
out2 = "calls_q30filteredvariants_biSNPs_AN__VAF.vcf"

########################################################


open(out, 'w').close()
open(out2, 'w').close()

with open(vcf, 'r') as infile:
	for i in infile:											# linewise
		onerow = i.rstrip("\n")
		onerow = onerow.split("\t")
		if len(onerow) > 3:
			if len(onerow[3]) == 1 and len(onerow[4]) == 1: 	# SNPs only
				chr = onerow[0]
				pos = onerow[1]
				SNPnotes = onerow[7]
				if "DP4=" in SNPnotes:
					parts = SNPnotes.split(";")
					DP = parts[0][3:]
					for n in range(len(parts)):
						if "DP4=" in parts[n]:
							VAF = parts[n]
							VAF_parts = VAF.split(",")
							DP4_1 = float(VAF_parts[0][4:])
							DP4_2 = float(VAF_parts[1]) 
							DP4_3 = float(VAF_parts[2]) 
							DP4_4 = float(VAF_parts[3]) 
							#print(DP4_1)
							#print(DP4_2)
							#print(DP4_3)
							#print(DP4_4)
							VAF = (DP4_4 + DP4_3) / (DP4_1+DP4_2+DP4_3+DP4_4)
							#print(VAF)
							with open(out, 'a') as outfile:
								outfile.write(str(chr) + "\t" + str(pos) + "\t" + str(DP) + "\t" + str(VAF)+ "\n")
							SNPnotes= SNPnotes.replace("AF=1","AF="+str(VAF))
							SNPnotes= SNPnotes.replace("AF=0.5","AF="+str(VAF))
							with open(out2, 'a') as outfile2:
								outfile2.write(onerow[0] + "\t" + onerow[1] + "\t" + onerow[2] + "\t" + onerow[3] + "\t" + onerow[4] + "\t" + onerow[5] + "\t" + onerow[6] + "\t" + SNPnotes + "\t" + onerow[8] + "\t" + onerow[9] + "\n")


infile.close()
outfile.close()

print("\ndone.")