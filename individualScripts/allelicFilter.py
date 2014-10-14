#!/usr/bin/env python

'''This program splits BAM files processed by suspenders: reads are written into different output alignment files (.bam) depending on their genomic origin (which is noted in the po:i tag).'''
'''If a third output file is indicated, reads which map equally well __at the same position_ in both genomes will also be reported (po:i:3).'''
'''Currently, this script ignores reads that were assigned by suspender's "RANDOM" filter'. '''
'''Possible improvements: let user choose which reads should be ignored (not just RANDOM reads, perhaps); allow SAM format for in- and output; multiprocessing'''

import sys
import argparse
import getopt
import pysam
import os

def get_args():
    parser=argparse.ArgumentParser(description='Split suspender BAM files')
    parser.add_argument('--BAMfile', '-in', type=str, required=True, help="suspender-generated BAM file")
    parser.add_argument('--outfile1', '-alt1', type=str, required=True, help = 'Name for output file with alternative genome 1')
    parser.add_argument('--outfile2', '-alt2', type=str, required=True, help = 'Name for output file with alternative genome 2')
    parser.add_argument('--outfile3', '-neither', type=str, help = 'Name for output file with reads that mapped equally well to both alternative genomes (optional).')
    parser.add_argument('--removeMultiMapped', '-vMulti', action='store_true', default=False, help = 'exclude reads that have more than one alignment based on the NH and/or XS tag')

    args=parser.parse_args()
    return args

    
def get_basenames(FileName):
    '''extracts the basename of the filenames for easy formating of the read count output'''
    Name = os.path.splitext(os.path.basename(FileName))[0]
    return(Name)


def remove_multiAlignments(ReadTagsEntry, KeepInfo):
    '''Checks for NH and XS tag to remove reads that aligned more than once'''    
    keepRead = KeepInfo
    if 'XS' in ReadTagsEntry:
        keepRead=False
    elif 'NH' in ReadTagsEntry and ReadTagsEntry[1] > 1:
        keepRead=False
    return(keepRead)


def main():

    args = get_args()

    infile = pysam.Samfile(args.BAMfile, "rb")
    out1 = pysam.Samfile(args.outfile1 , "wb", template = infile )
    out2 = pysam.Samfile(args.outfile2 , "wb", template = infile )
    if args.outfile3:
        out3 = pysam.Samfile(args.outfile3 , "wb", template = infile )

    keep = True
    origin = 0
    count1 = 0
    count2 = 0
    count3 = 0

    for DNAread in infile.fetch():
        intags = DNAread.tags

        for entry in intags:

            if 'po' in entry and entry[1] == 2:
                origin = 2
            elif 'po' in entry and entry[1] == 1:
                origin = 1
            
            if 'ct' in entry and entry[1]== 'R':
                keep = False

            if args.removeMultiMapped:
                keep = remove_multiAlignments(entry, keep)
        
        if(origin == 1 and keep):
            out1.write(DNAread)
            count1 += 1
        elif(origin == 2 and keep):
            out2.write(DNAread)
            count2 += 1
        elif keep:
            count3 += 1
            if args.outfile3:
                out3.write(DNAread)

        origin = 0
        keep = True

    out1.close()
    out2.close()

    CountOut="reads mapped to alternative genome 1 ({1}): {0}\nreads mapped to alternative genome 2 ({3}): {2}\nundistinguishable reads {4}".format(count1, get_basenames(args.outfile1), count2, get_basenames(args.outfile2), count3)
    print CountOut


if __name__ == '__main__':
    main()



