#!/usr/bin/python

import sys, getopt


def main(argv):

        try:
               opts, args = getopt.getopt(argv,"hi:o:",["help","file1=","ofile="])
        except getopt.GetoptError:
                print 'vcf2vcf.py -i <vcf_file> -o <outfile>'
                sys.exit(2)

        for opt, arg in opts:
                if opt == '-h':
                        print 'vcf2vcf.py -i <vcf_file> -o <outfile>'
                        sys.exit()
                elif opt in ("-i", "--file1"):
                        file1 = open( arg, 'r')
                elif opt in ("-o", "--ofile"):
                        outputfile = open( arg, 'w')



        for line in file1:
                if (line[0] != "#"):
                        entry = line.rstrip().split("\t")
                        ref = entry[3]
                        alt = entry[4]
			info = entry[7]
			genotype = entry[9]
			
			if ( len(alt.split(",")) > 1 ):		# if there are more than one alternative genotypes: take only the occuring alternative; change also INFO and GT column to have only the values for the existing alternative genotype
				ac_entries = info.split("AC=")

				alternatives = alt.split(",")
				if(int(genotype[0]) == 0):
					alt = alternatives[0]
					genotype = "0/0:0:0"
					info = ac_entries[0] + "AC=" + ac_entries[1].split(",")[0]
				else:
					alt = alternatives[int(genotype[0])-1]
					info = ac_entries[0] + "AC=" + ac_entries[1].split(",")[int(genotype[0])-1]
					genotype = "1/1:" + entry[9].split(":")[1].split(",")[int(genotype[0])-1] + ":" + entry[9].split(":")[2].split(",")[int(genotype[0])-1]


			

			if (ref[0] == alt[0] and len(ref) != len(alt) and len(ref) == 1):		# insertion (first base equal in ref and alt)
				outputfile.write(entry[0] + "\t" + entry[1] + "\t" + entry[2] + "\t" + ref + "\t" + alt + "\t" + entry[5] + "\t" + entry[6] + "\t" + info + "\t" + entry[8] + "\t" + genotype + "\n")


			elif (ref[0] == alt[0] and len(ref) != len(alt) and len(alt) == 1):		# deletion (first base equal in ref and alt)
				outputfile.write(entry[0] + "\t" + entry[1] + "\t" + entry[2] + "\t" + ref + "\t" + alt + "\t" + entry[5] + "\t" + entry[6] + "\t" + info + "\t" + entry[8] + "\t" + genotype + "\n")


			elif (ref[0] == alt[0] and len(ref) == len(alt)):				# CC -> CA (MBS)
				#print line
				for i in range(1,len(ref)):
	                        	outputfile.write(entry[0] + "\t" + str(int(entry[1]) + i) + "\t" + entry[2] + "\t" + ref[i] + "\t" + alt[i] + "\t" + entry[5] + "\t" + entry[6] + "\t" + info + "\t" + entry[8] + "\t" + genotype + "\n")


			elif (ref[0] == alt[0] and len(ref) > len(alt)):				# CCGTC -> CA
				#print line
				for i in range(1,len(alt)):
					outputfile.write(entry[0] + "\t" + str(int(entry[1]) + i) + "\t" + entry[2] + "\t" + ref[i] + "\t" + alt[i] + "\t" + entry[5] + "\t" + entry[6] + "\t" + info + "\t" + entry[8] + "\t" + genotype + "\n")

				outputfile.write(entry[0] + "\t" + str(int(entry[1]) + len(alt) - 1) + "\t" + entry[2] + "\tX" + ref[len(alt):len(ref)] + "\tX\t" + entry[5] + "\t" + entry[6] + "\t" + info + "\t" + entry[8] + "\t" + genotype + "\n")


			elif (ref[0] == alt[0] and len(ref) < len(alt)):				# CCGTC -> CAATGCATA
				#print line
				for j in range(1,len(ref)):
					outputfile.write(entry[0] + "\t" + str(int(entry[1]) + j) + "\t" + entry[2] + "\t" + ref[j] + "\t" + alt[j] + "\t" + entry[5] + "\t" + entry[6] + "\t" + info + "\t" + entry[8] + "\t" + genotype + "\n")

				outputfile.write(entry[0] + "\t" + str(int(entry[1]) + len(ref) - 1) + "\t" + entry[2] + "\tX\tX" + alt[len(ref):len(alt)] + "\t" + entry[5] + "\t" + entry[6] + "\t" + info + "\t" + entry[8] + "\t" + genotype + "\n")


			elif (len(ref) == 1 and len(alt) == 1):						# SNP
				outputfile.write(entry[0] + "\t" + entry[1] + "\t" + entry[2] + "\t" + ref + "\t" + alt + "\t" + entry[5] + "\t" + entry[6] + "\t" + info + "\t" + entry[8] + "\t" + genotype + "\n")


			elif (ref[0] != alt[0] and len(ref) == len(alt)):				# multiple bases substitution -> transform to SNPs
				#print line
				for i in range(0,len(ref)):
	                        	outputfile.write(entry[0] + "\t" + str(int(entry[1]) + i) + "\t" + entry[2] + "\t" + ref[i] + "\t" + alt[i] + "\t" + entry[5] + "\t" + entry[6] + "\t" + info + "\t" + entry[8] + "\t" + genotype  + "\n")



			elif (ref[0] != alt[0] and len(ref) > len(alt)):				# ACAAG > TT  mixed SNP and deletion
				#print line
				for i in range(0,len(alt)):
					outputfile.write(entry[0] + "\t" + str(int(entry[1]) + i) + "\t" + entry[2] + "\t" + ref[i] + "\t" + alt[i] + "\t" + entry[5] + "\t" + entry[6] + "\t" + info + "\t" + entry[8] + "\t" + genotype + "\n")

				outputfile.write(entry[0] + "\t" + str(int(entry[1]) + len(alt) - 1) + "\t" + entry[2] + "\tX" + ref[len(alt):len(ref)] + "\tX\t" + entry[5] + "\t" + entry[6] + "\t" + info + "\t" + entry[8] + "\t" + genotype + "\n")


			elif (ref[0] != alt[0] and len(alt) > len(ref)):				# TT > ACAAG  mixed SNP and insertion
				#print line
				for j in range(0,len(ref)):
					outputfile.write(entry[0] + "\t" + str(int(entry[1]) + j) + "\t" + entry[2] + "\t" + ref[j] + "\t" + alt[j] + "\t" + entry[5] + "\t" + entry[6] + "\t" + info + "\t" + entry[8] + "\t" + genotype + "\n")

				outputfile.write(entry[0] + "\t" + str(int(entry[1]) + len(ref) - 1) + "\t" + entry[2] + "\tX\tX" + alt[len(ref):len(alt)] + "\t" + entry[5] + "\t" + entry[6] + "\t" + info + "\t" + entry[8] + "\t" + genotype  + "\n")

			else:
				print line

                else:
                        outputfile.write(line)


	file1.close()
	outputfile.close()



if __name__=='__main__':
        main(sys.argv[1:])

