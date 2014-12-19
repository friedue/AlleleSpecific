Allele-specific mapping of short reads
======================================

Based on [lapels](https://code.google.com/p/lapels/) and [suspenders](https://pypi.python.org/pypi/suspenders/), we want to establish a Galaxy-based pipeline for allele-specific mapping of RNA-seq and DNA-seq applications.

Please see those publications for details on lapels, suspenders and MOD files:

[Holt, J., Huang, S., McMillan, L., & Wang, W. (2013). Read Annotation Pipeline for High-Throughput Sequencing Data. In J. Gao (Ed.), ACM Conference on Bioinformatics, Computational Biology and Biomedical Informatics. Washington, DC, USA.](http://csbio.unc.edu/CCstatus/Media/suspenders_acmbcb2013.pdf)

[Huang, S., Kao, C.-Y., McMillan, L., & Wang, W. (2013). Transforming Genomes Using MOD Files with Applications. In J. Gao (Ed.), ACM Conference on Bioinformatics, Computational Biology and Biomedical Informatics. Washington, DC, USA.](http://csbio.unc.edu/CCstatus/Media/mod_acmbcb2013.pdf)

[Huang, S., Holt, J., Kao, C.-Y., McMillan, L., & Wang, W. (2014). A novel multi-alignment pipeline for high-throughput sequencing data. Database : The Journal of Biological Databases and Curation, 2014(0).](http://dx.doi.org/10.1093/database/bau057)

----------------------------------------------
### Input

##### files 

* FASTQ files of raw reads

```
    FASTQ_FOLDER=(fastqs/ ) # assumes that fastq files end with .fastq.gz and paired-end sequencing, (_R1, _R2) 
```

* MOD file -- if the MOD file is not available, it must be created or downloaded
  * check whether your MOD file is available [here](http://www.csbio.unc.edu/CCstatus/index.py?run=Pseudo "MOD files supplied by Huang et al.") 
  * otherwise, you will absolutely need VCF files with the information for your strain(s) of interest and their differences to the reference genome and the reference genome FASTA file:

```
    REF_FASTA=ce10.fa # fasta file for reference genome, e.g. mm9, dm3, ce10....
    VCF_SNP=snps.vcf # must contain the information about genetic variants for at least both strains of interest
    VCF_INDELS=indels.vcf # like snps.vcf, must contain information about the genetic variations
```

* for RNA-seq: GTF file is needed, too
    
##### information

these variables should be set regardless whether a MOD file is available already or not

```
    REF_GENOME=ce10 # reference genome
    MAT_STRAIN=A # maternal strain --> must match entry from to vcf column
    PAT_STRAIN=B # paternal strain --> must match entry from vcf column
    CHROMOSOMES=1,2,3,X --> in case, one is interested in a couple of chromosomes only
```

##### programs

###### MOD file generation

lapels comes with a set of supplementary scripts that can be used to generate the MOD file if needed:

1. get_refmeta
2. vcf2mod
3. insilico

* custom-made script: changeChrNamesInMODFile.awk 
* additional tool: tabix

###### After MOD file generation

These scripts contain the core of the lapels/suspenders package:

4. __modmap__: lifting over coordinates of _annotation_ files (e.g. BED, GTF)
5. __lapels__: lifting over the coordinates of _read_ files (BAM)
6. __suspenders__: merging of 2 BAM files that were aligned to different "pseudogenomes"

* additional tools: bowtie2/tophat - for read alignment
* custom-made scripts:
 * allelicFilter.py
 * additional script for insert-size determination (needed only if TopHat is used for alignment)

![Overview](https://raw.githubusercontent.com/friedue/AlleleSpecific/master/images/pipelineOverview.png)

Image modified from [Huang et al., 2014](http://dx.doi.org/10.1093/database/bau057) to illustrate the different steps of our pipeline

-----------------------------------------------

Part I: Preparing pseudogenomes - this will only need to be done once per genome version!
---------------

###  I.1. Meta-data for reference genome

-->  lapels/bin/get_refmeta

      lapels-1.0.5/bin/get_refmeta -o ${REF_GENOME}.meta ${REF_GENOME} ${REF_FASTA}


### I.2. generating MOD files

##### a) indexing of VCF files
   
	for vcf in ${VCF_INDELS} ${VCF_SNP}
	do
          tabix-0.2.6/bgzip -c ${vcf} > ${vcf}.gz
    	  tabix-0.2.6/tabix -p vcf ${vcf}.gz
    	done

##### b) VCF to MOD format

--> lapels/bin/vcf2mod

1 MOD file per VCF file (SNPs, INDELS) and genotype (maternal, paternal)
SNPs with bad quality (FI tag = 0) will be discarded
    
    lapels-1.0.5/bin/vcf2mod -c ${CHROMOSOMES} -f \
		-o ${REF_GENOME}_SNPs_${genotype}.mod \
		${REF_GENOME} ${REF_GENOME}.meta \
		${genotype} ${VCF_SNP} 2> ${genotype}.vcf2mod_snps_`date +"%Y%m%d%H%M%S"`.log &

##### c) merge MOD files for SNPs and indels per genotype

    cat ${REF_GENOME}_SNPs_${genotype}.mod ${REF_GENOME}_indels_${genotype}.mod |\
	sort -k2,2n -k3,3n | uniq | awk -f changeChrNamesInMODFile.awk - > ${REF_GENOME}_indels_SNPs_${genotype}_changedChr.mod 

### I.3. Generating pseudogenomes

--> lapels/bin/insilico

CAVE: after this step, the MOD file will be gzipped (without any indication in the file name)
    
	lapels-1.0.5/bin/insilico \
		${REF_GENOME}_indels_SNPs_${genotype}_changedChr.mod \
		${REF_FASTA} -v -f \
		-o pseudogenome_${REF_GENOME}_${genotype}.fa >${genotype}.insilico_`date +"%Y%m%d%H%M%S"`.log 

---------------------------------------

Part II: Mapping - this needs to be done for every sample
------------------

### II.1. Bowtie/Tophat: mapping to each pseudogenome

In principle, this step can be done with any read alignment programme. Here I show the routine for genome and transcriptome data that is aligned with either bowtie2 or tophat.

#### Non-RNA-seq experiment

--> bowtie2

**building bowtie indeces**

	bowtie2-2.1.0/bowtie2-build pseudogenome_${REF_GENOME}_${genotype}.fa pseudogenome_${REF_GENOME}_${genotype} > ${genotype}.bowtieIndex_`date +"%Y%m%d%H%M%S"`.log 2> ${genotype}.bowtieIndex_Error_`date +"%Y%m%d%H%M%S"`.log

**mapping with bowtie2** (example is shown for paired-end sequencing experiment (-1, -2) with maximum distance between the mates of 1 kb (-X)

	(/bowtie2-2.1.0/bowtie2 -x pseudogenome_${REF_GENOME}_${genotype} \
		-1 ${fastq_folder}${sample}_R1.fastq.gz -2 ${fastq_folder}${sample}_R2.fastq.gz \
		-X 1000 -p 40 \
		--rg-id mpi-ie --rg CN:deep_sequencing_unit --rg PL:illumina | \
	samtools view -Sb - |\
	samtools sort - ${genotype}_${sample}.sorted) 2>&1 > MappingTo${genotype}_${sample}.log

    samtools index ${genotype}_${sample}.sorted.bam &
  
#### RNA-seq: tophat

In contrast to genome alignment, transcriptome alignment programs often benefit from supplying a file containing known gene-coordinates (e.g. a GTF file). This file must be present for each pseudogenome.

* a) generate GTFs for maternal and paternal versions: **modmap**
* b) transcriptome build for both genotypes using 1 set of FASTQs with **TopHat**
* c) run TopHat on remaining samples

###### a) GTF for pseudogenome

--> lapels/bin/modmap

modmap expects 0-based positions, since gtf files are 1-based, I first convert them into 0-based, use modmap to change the annotation to match the individual genomes, and turn the resulting 0-based gtf file back into 1-based (additionally, I remove the negative numbers introduced by modmap for regions that fall into deletions)

    awk -F "\t" '{OFS="\t"; print $1,$2,$3,$4-1, $5-1, $6, $7, $8,$9}' ${REF_GENOME}.gtf > ${REF_GENOME}_0based.gtf
    
    lapels-1.0.5/bin/modmap -d '\t' \
		-f ${REF_GENOME}_indels_SNPs_${genotype}_changedChr.mod	\
		${REF_GENOME}_0based.gtf ${REF_GENOME}_${genotype}_0based.gtf 1,4 1,5

the awk magic here is neccessary to turn the negative values into positive ones so that the GTF file does not contain negative chromosome coordinates 

     awk -F "\t" '{OFS="\t";print $1,$2,$3, $4<0?$4*-1:$4+1,$5<0?$5*-1:$5+1,$6,$7,$8,$9}' ${REF_GENOME}_${genotype}_0based.gtf > ${REF_GENOME}_${genotype}_1based.gtf
    

###### b) building the transcriptome index for both genotypes using TopHat

* you'll need the insert sizes for the reads from RNA-seq experiments --> I used a script from Andreas that generated txt files named "insert_stats_RNA..." using Picard

extracting the values for standard deviation and inner distance from the text files:

    STD_DEV=`sed '1d' fastq/insert_stats/insert_stats_RNA__1.txt | sed -n '3p' | awk '{print sprintf("%.0f",$2)}'`
    INNER_DIST=`sed '1d' fastq/insert_stats/insert_stats_RNA_WT_1.txt | sed -n '3p' | awk '{print sprintf("%.0f",$1)}'`

running top hat _once per pseudogenome_ to obtain the transcriptome indeces:      

    tophat2 -G ${REF_GENOME}_${genotype}_1based.gtf \
	--transcriptome-index transcriptome_data/${REF_GENOME}_${genotype}_transcriptomeIndex \
    	--no-coverage-search --library-type fr-firststrand --mate-inner-dist ${INNER_DIST} \
    	--mate-std-dev ${STD_DEV} \
		pseudogenome_${REF_GENOME}_${genotype} \
    	-o RNA_WT_${genotype} \
	    SAMPLE_R1.fastq.gz \
    	    SAMPLE_R2.fastq.gz
    
###### c) looping over remaining files fastq files (if any)

        for fastq in ${FASTQ_FOLDER}*_R1.fastq.gz
        do
        
        # extracting std dev and inner distance from the output of Andreas' script
        forStats=$(basename "$fastq" .fastq.gz | sed 's/_R[1-2]//' | awk '{print substr($0,index($0,"_RNA")+1)}')
        STD_DEV=`sed '1d' fastq/insert_stats/insert_stats_${forStats}.txt | sed -n '3p' | awk '{print sprintf("%.0f",$2)}'`
        INNER_DIST=`sed '1d' fastq/insert_stats/insert_stats_${forStats}.txt | sed -n '3p' | awk '{print sprintf("%.0f",$1)}'`

	# extracting sample name        	
        sample=$(basename "$fastq" .fastq.gz | sed 's/_R[1-2]//')
        
        # aligning reads with tophat
        tophat2 -o ${sample}_${genotype} -p 20 \
        	-G ${REF_GENOME}_${genotype}_1based.gtf \
        	--transcriptome-index transcriptome_data/${genotype}_transcriptomeIndex \
        	--no-coverage-search --library-type fr-firststrand \
        	--mate-inner-dist ${INNER_DIST} --mate-std-dev ${STD_DEV} \
        	pseudogenome_${REF_GENOME}_${genotype} \
        	${FASTQ_FOLDER}/${sample}_R1.fastq.gz ${FASTQ_FOLDER}/${sample}_R2.fastq.gz
        
        # indexing aligned reads
        samtools index ${FOLDER}/accepted_hits.bam &
        done
  
### II.2. Translating the different pseudogenome mappings back to reference genome

This is the core functionality of the lapels package.

--> lapels/bin/pylapels

    lapels-1.0.5/bin/pylapels -p 50 -f \
    	-o ${bam}_mappedBackTo_${REF_GENOME}.bam \
    	${REF_GENOME}_indels_SNPs_${genotype}_changedChr.mod ${BAM} > ${bam}_stdout.txt 2> ${bam}_stderr.txt
    
 

### II.3 Merging pseudogenome mappings

--> suspenders/bin/pysuspenders

suspenders will determine the best fit for each read

input: 2 BAMs, one for each pseudogenome mapping (parental, maternal)

output: 1 BAM file with reads where flags indicate maternal or paternal origin

example for suspenders on RNA-seq data aligned with tophat:

     suspenders-0.2.4/bin/pysuspenders \
     	--lapels --quality -p 15 \
     	mergedAlignment_${SAMPLE}.bam \
     	tophat_${SAMPLE}_${MAT_STRAIN}/accepted_hits_pylapels.bam \
     	tophat_${SAMPLE}_${PAT_STRAIN}/accepted_hits_pylapels.bam

-----------------------------------------

Part III: Counting allele-specific reads
-----------------------------------------

### III BAM file filter, merge, sort, index

Suspender's output needs to be modified for most down-stream analyses:

   * __filtering__ for multiple alignments (at read alignment step), randomly assigned reads (by suspenders), chrM, random chromosomes
   * __sorting__ - suspenders outputs a BAM file sorted by read name, this cannot be indexed and is thus unsuitable for most down-stream analyses --> sorting by coordinates
   * __splitting__ into separate files for maternal and paternal reads

---> allelicFilterOct2014.py does all of that 

     for BAM in mergedAlignment*.bam
     do
     	SAMPLE=$(basename "$BAM" .bam)
     	echo ${SAMPLE}
     	python allelicFilterOct2014.py -in ${BAM} \
        	--outfile1 ${SAMPLE}_maternal --outfile2 ${SAMPLE}_paternal \ # optionally, a third output file for those reads that cannot be tracked back to either pseudogenome can also be indicated
     		-vMulti \ # exclusion of multiple alignments determined by XS and NH tags
     		--coordinateSorting > counts_${SAMPLE}.txt & \ # counts*.txt will contain statistics about the output files
     done

indexing of coordinate-sorted BAM files

        for BAM in *.bam
        do
         samtools index ${BAM}
        done   
