Given the presence of fasta files with the information about the pseudogenomes, the pipeline is as follows:

## 1. MAPPING to pseudogenomes

||Bowtie|TopHat|
|-----|------|-----|
|**In** | indeces for pseudo-genome FASTA file | indeces for pseudo-genome FASTA file, GTF file with coordinates matching the respective pseudo-genome |
|**Out** | coordinate-sorted BAM | coordinate-sorted BAM|

## 2. LAPELS

pseudogenome coordinates --> reference genome coordinates

- **In**: MOD file, 1 coordinate-sorted BAM file
- **Out**: 1 coordinate-sorted BAM file with reference genome coordinates

## 3. SUSPENDERS

merging the two different mappings, identifying most probable origin for each read

- **In**: 2 read-name-sorted BAM files
- **Out**: 1 read-name-sorted BAM file with new flags indicating the origin of each read

## 4. SPLITTING

- **In**: 1 (preferably read-name-sorted) BAM file with flags set by suspenders, bowtie and/or tophat
- **Out**: 2-3 (if wanted coordinate-sorted) BAM files - one for each origin (maternal, paternal, ambigious)

=====================================================================
__see this [script on the mapping pipeline](https://github.com/friedue/AlleleSpecific/blob/master/Notes/AlleleSpec_Pipeline.md) for the details__

=====================================================================

## 5. POST-PROCESSING

e.g. count tables of reads per genomic regions, bigWig files etc.

see my notes on [BAM Post-Processing](https://github.com/friedue/AlleleSpecific/blob/master/Notes/BAMProcessing.md)
