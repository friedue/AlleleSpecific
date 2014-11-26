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

* either with modified deepTools code or with edgeR

<a name="countTables"></a>
## 4. Count tables

* featureCounts
* gene region file (SAF)
* for use with edgeR or DESeq
