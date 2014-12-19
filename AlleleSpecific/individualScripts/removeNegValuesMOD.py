#!/usr/bin/python

import sys, getopt

def main(argv):

        try:
               opts, args = getopt.getopt(argv,"hi:o:",["help","mpileupfile=","jfile=","snpfile=","ofile="])
        except getopt.GetoptError:
                print 'removeNegValuesMOD.py -i <infile> -o <output_file>'
                sys.exit(2)
        for opt, arg in opts:
                if opt == '-h':
                        print 'USAGE: removeNegValuesMOD.py -i <infile> -o <output_file>'
                        sys.exit()
                elif opt in ("-i", "--infile"):
                        infile = open( arg, 'r')
                elif opt in ("-o", "--ofile"):
                        outputfile = open( arg, "w")


	for line in infile:
		entry=line.rstrip().split("\t")
		if(len(entry) == 9):
			start = entry[3]
			end = entry[4]

			if(start[0] == "-" and end[0] == "-"):
				#outputfile.write(line.replace(start,str(int(start[1:len(start)])+1)).replace(end,str(int(end[1:len(end)])+1)))
				outputfile.write(entry[0] + "\t" + entry[1] + "\t" + entry[2] + "\t" + str(int(start[1:len(start)])+1) + "\t" + str(int(end[1:len(end)])+1) + "\t" + entry[5] + "\t" + entry[6] + "\t" + entry[7] + "\t" + entry[8] + "\n")
			elif (start[0] == "-"):
				#outputfile.write(line.replace(start,str(int(start[1:len(start)])+1)).replace(end,str(int(end)+1)))
				outputfile.write(entry[0] + "\t" + entry[1] + "\t" + entry[2] + "\t" + str(int(start[1:len(start)])+1) + "\t" + str(int(end)+1) + "\t" + entry[5] + "\t" + entry[6] + "\t" + entry[7] + "\t" + entry[8] + "\n")
			elif (end[0] == "-"):
				#outputfile.write(line.replace(end,str(int(end[1:len(end)])+1)).replace(start,str(int(start)+1)))
				outputfile.write(entry[0] + "\t" + entry[1] + "\t" + entry[2] + "\t" + str(int(start)+1) + "\t" + str(int(end[1:len(end)])+1) + "\t" + entry[5] + "\t" + entry[6] + "\t" + entry[7] + "\t" + entry[8] + "\n")
			else:
				#outputfile.write(line)
				outputfile.write(entry[0] + "\t" + entry[1] + "\t" + entry[2] + "\t" + str(int(start)+1) + "\t" + str(int(end)+1) + "\t" + entry[5] + "\t" + entry[6] + "\t" + entry[7] + "\t" + entry[8] + "\n")
		else:
			print line



	infile.close()
        outputfile.close()


if __name__=='__main__':
        main(sys.argv[1:])


# awk '{print $4,$5}' Drosophila_melanogaster.BDGP5.74.indv2.negRem.gtf | grep "-" | sort | uniq
#2048873	-2048887
#20801699	-20801709
