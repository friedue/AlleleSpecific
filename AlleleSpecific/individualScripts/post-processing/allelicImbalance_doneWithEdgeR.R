#! /package/R-3.1.0/bin/Rscript --vanilla
#########################################
### November 2014, FrieDue
### This script takes a composite feature count table
### (i.e. multiple BAM file counts for the same set of genome regions)
### and returns one bedgraph file for the indicated sample.
### It requires rather strict naming conventions (check the
### gsub routines here and in the functions; e.g. it can currently
### only handle 129S1 and CASTEiJ as the alleles).
### Since it makes use of the data.table package, even very
### large featureCount tables (> 1 mio rows) can be handled.
#########################################
### args: [1] path to feature count table, e.g. featCount_genomeCounts.txt; that file must contain 1 arbitratry line on top (e.g. beginning with #), 2nd line should contain the header
###       [2] sample # important for naming convention and file name of the output
###       [3] cell type # same as for [2]
######################################### 

source("f_edgeRForAllelicImbalance.R") ## can be found in the same github folder
library(edgeR) # must be installed first
library(data.table) # must be installed first

args <- commandArgs(TRUE)

# reading in
### assumes the following naming convention of the BAM files used with featureCounts that are present 
### in the header of the featureCount read count table:
### SAMPLEREPLICATE_CELL_129S1.ANYTHING.bam or SAMPLE_CELL_CASTEiJ.ANYTHING.bam, e.g.
### InputA_myCell_129S1.awesomeReads.bam
In =  ReadingInFeatureCountData.edgeR.genome(FileName=args[1], Sample = args[2], Cell = args[3])
print("Running EdgeR now")

# pairwise test
ptm <- proc.time()
In.edgeR = RunningEdgeR(calculateSizeFactor=FALSE, Counts.DF=In$forEdgeR,
                            Group = gsub("(.*)_([0-9AB]*)_(.*)_(.*)", "\\4", names(In$forEdgeR)), # this assumes that the sample names could be extracted from the featureCount table into the following format: SAMPLE_REPLICATE_CELL_GENOTYPE
                            PairwiseTest=TRUE, ExpDesign="Genotype")
proc.time()-ptm
print("Done with EdgeR; calculating allelic imbalance score from edgeR'slogFC")

# calculating AI
TT = topTags(In.edgeR$ExactTest, n=Inf, sort.by = "none")
TT.ai = AIFromlog2FC.df(as.data.frame(TT))

print("formatting everything to match bedgraph specifications")

# obtaining chr, start, end for bedgraph
Bins = as.data.table(In$Bins[which(In.edgeR$kept),])

# merging bedgraph info and AI score
TT.ai$Bin = row.names(TT.ai)
TT.ai = as.data.table(TT.ai)
setkey(Bins, Geneid)
setkey(TT.ai, Bin)
TT.out = TT.ai[Bins, nomatch=FALSE]

# writing out a bedgraph file
write.table(TT.out[,c("Chr","Start","End","AI"), with=FALSE], file = paste(args[2], args[3], "AI_edgeR.bedgraph", sep ="_"),
            sep = "\t", quote=F, col.names=F, row.names=F)
