Once the pipeline as described in [AlleleSpec_Pipeline.md](https://github.com/friedue/AlleleSpecific/blob/master/Notes/AlleleSpec_Pipeline.md) has been completed, I generate several files for different aspects of downstream analyses.

1. [__bigWigs of "_optimally aligned reads_"__](#optimal)
2. [__bigWigs of _maternal_ or _paternal_ reads__](#MatPat)
3. [__bigWigs of allelic imbalances__](#AI)
4. [__count tables of maternal and paternal reads for determining regions of significant allelic differences__](#countTables)

===============================

<a name="optimal"></a>
## 1. Optimally aligned reads

* bamCoverage and/or bamCompare

<a name="MatPat"></a>
## 2. Maternal/paternal reads

* bamCoverage and/or bamCompare

<a name="AI"></a>
## 3. Allelic imbalance

I define allelic imbalance based on the ratio of read counts for the maternal vs. the paternal allele. 

```
AI = (Nmat - Npat) / (Nmat + Npat)
```

* -1 = paternal origin
* 0 = equal number of reads from mat and pat
* 1 = maternal origin

* I tried two ways for the genome-wide calculation: a [modified version of bamCoverage](#AIbamCov) from deepTools and [edgeR](#AIedgeR)

<a name="AIbamCov"></a>
### allelic imbalance with bamCoverage

<a name="AIedgeR"></a>
### allelic imbalance with edgeR

I wanted to make use of edgeR's more sophisticated way of calculating logFC between the maternal and paternal alleles. Since edgeR's input is a read count matrix, I first separated the genome into 100 bp bins, used featureCounts for the read counts, and let edgeR output a bedgraph file of logFC per bin.

1. __genome tiling__

--> [tiling.pl](https://github.com/friedue/randomScripts/blob/master/tiling.pl "genome tiling script")

```
tiling.pl -length 100 -shift 101 mm9.fa.fai > mm9.tiled100bp.bed
awk '{OFS="\t"; print "bin"NR, $1,$2,$3,$4,"+"}' mm9.tiled100bp.bed > mm9.tiled100bp.saf # saf format for featureCounts
```

2. __featureCounts__

--> [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) from the subread package 

```
subread-1.4.5-p1-source/bin/featureCounts --ignoreDup -f -O -Q 20 -p -T 35 -a mm9.tiled100bp.saf -F SAF \
                                          -o featCount_${OUT}.txt \
                                          ${BAMFOLDER}/*markedDupl.bam \ # this indicates that all files ending with .markedDupl.bam will be used by featureCount
                                          2> featCount_${OUT}.log
```

3. __edgeR__

The featureCount table will be fed to edgeR which will assess the read counts of the replicates for each sample and will compare them for the different alleles. edgeR will calculate a logFC between both alleles (using CASTEiJ as the reference "condition") which I will turn into the allelic imbalance score described above.

--> edgeR must be installed in R
--> [allelicImbalance_doneWithEdgeR.R](https://github.com/friedue/AlleleSpecific/blob/master/individualScripts/post-processing/allelicImbalance_doneWithEdgeR.R) script for optimized reading in and handling of large featureCount tables and the corresponding [R functions](https://github.com/friedue/AlleleSpecific/blob/master/individualScripts/post-processing/f_edgeRForAllelicImbalance.R)

```
./allelicImbalance_doneWithEdgeR.R featCount_genomeCounts.txt ${SAMPLE} ${CELL}  # running edgeR, getting bedgraph output
```

the output is a bedGraph-style file 

4. __bedGraph to bigWig__

--> UCSCtools: bedGraphToBigWig

```
# sorting by coordinates
for FILE in *AI_edgeR.bedgraph
do
  sort -k1,1 -k2,2n ${FILE} > ${FILE}.sort &
done
wait

# bedGraph to bigWig
for FILE in *AI_edgeR.bedgraph.sort
do
  /package/UCSCtools/bedGraphToBigWig ${FILE} mm9.fa.fai ${FILE}.bw &
done
```

<a name="countTables"></a>
## 4. Count tables

* featureCounts
* gene region file (SAF)
* for use with edgeR or DESeq
