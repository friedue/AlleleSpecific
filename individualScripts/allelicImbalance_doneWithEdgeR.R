#! /package/R-3.1.0/bin/Rscript --vanilla
#########################################
### This script takes a composite feature count table
### (i.e. multiple BAM file counts for the same set of genome regions)
### and returns one bedgraph file for the indicated sample.
### It requires rather strict naming conventions (check the
### gsub routines here and in the functions).
### Since it makes use of the data.table package, even very
### large featureCount tables (> 1 mio rows) can be handled.
#########################################
### args: [1] path to feature count table, e.g. featCount_genomeCounts.txt; must contain 1 arbitratry line on top, 2nd line should contain the header
###       [2] sample
###       [3] cell type
#########################################

source("/data/projects_2/muehlpfordt/2014_AlleleSpecificMapping/2014-10-28_edgeR_for_allelicImbalance/f_edgeRForAllelicImbalance.R")
library(edgeR)
library(data.table)

args <- commandArgs(TRUE)

# reading in
In =  ReadingInFeatureCountData.edgeR.genome(FileName=args[1], Sample = args[2], Cell = args[3])
print("Running EdgeR now")

# pairwise test
ptm <- proc.time()
In.edgeR = RunningEdgeR(calculateSizeFactor=FALSE, Counts.DF=In$forEdgeR,
                            Group = gsub("(.*)_([0-9AB]*)_(.*)_(.*)", "\\4", names(In$forEdgeR)),
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
