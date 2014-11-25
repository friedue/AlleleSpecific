Given a MOD file, the pipeline is as follows:

## 1. MAPPING to pseudogenomes

||Bowtie|TopHat|
|-----|------|-----|
|**Input** | indeces, pseudo-genome | indeces, pseudo-genome, GTF file with coordinates matching pseudo-genome |
|**Output** | coordinate-sorted BAM | coordinate-sorted BAM|

## 2. LAPELS

pseudogenome coordinates --> reference genome coordinates

- **In**: MOD file, 1 coordinate-sorted BAM file
- **Out**: 1 coordinate-sorted BAM file with reference genome coordinates

## 3. SUSPENDERS

- **In**: 2 read-name-sorted BAM files
- **Out**: 1 read-name-sorted BAM file with new flags indicating the origin of each read

## 4. SPLITTING

- **In**: 1 (preferably read-name-sorted) BAM file with flags set by suspenders, bowtie and/or tophat
- **Out**: 2-3 (if wanted coordinate-sorted) BAM files - one for each origin (maternal, paternal, ambigious)
