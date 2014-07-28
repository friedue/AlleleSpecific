#!/usr/bin/python

import sys
import optparse
import getopt

def main(argv):
        
	try:
               opts, args = getopt.getopt(argv,"hr:a:",["help","ofile1=","ofile2="])
        except getopt.GetoptError:
                print 'allelicFilterBAM.py -r <outfile1> -a <outfile2>' 
                sys.exit(2)
        for opt, arg in opts:
                if opt == '-h':
                        print 'This program splits BAM files processed by suspenders: reads are written into different output alignment files depending on their genomic origin (po:i tag). Reads which map equally well at the same position in both genomes are not reported (po:i:3). Reads reported by random filter are ignored!! \n\n /package/samtools/samtools view -h <bam file> | allelicFilterBAM.py -r <po:i:1_outfile> -a <po:i:2_outfile>'
                        sys.exit()
                elif opt in ("-r", "--ofile1"):
                        outputfile1 = open( arg, "w")
                elif opt in ("-a", "--ofile2"):
                        outputfile2 = open( arg, "w")
	
	origin = 0
	keep = True
	for line in sys.stdin:
	    
        	field = line.rstrip().split("\t")

       		if line[0] == '@':
			outputfile1.write(line)
			outputfile2.write(line)
    	      		continue


        	elif len(field) < 12: # checking the tag fields (after field 12)
			keep = false

            	for subfield in field[12:]:
                	(name,ftype,value) = subfield.split(":")
                	if name=='po' and int(value) == 1: # this will only work for two input files --> in theory, suspenders can compare more than 2 files (For example, given three inputs and an identical read in the second and third input, the flag would be 0b110 which would be stored as integer value 6.)
                		origin = 1
			elif name=='po' and int(value) == 2:
				origin = 2
			elif name=='ct' and value=='R': # ct:R = randomly assigned reads --> will be ignored
				keep = False

		if(origin == 1 and keep):
 	               	outputfile1.write(line)
		elif(origin == 2 and keep):
			outputfile2.write(line)

		origin = 0
		keep = True


    	outputfile1.close()
    	outputfile2.close()


if __name__ == '__main__':
    main(sys.argv[1:])



