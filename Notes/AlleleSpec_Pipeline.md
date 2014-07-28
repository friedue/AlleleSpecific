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

    REF_FASTA=mm9.fa
    VCF_SNP=mgp.v2.snps.annot.reformat.vcf
    VCF_INDELS=mgp.v2.indels.annot.reformat.vcf
    FASTQ_FOLDER=(fastqs/ ) # assumes that fastq files end with .fastq.gz and paired-end sequencing, (_R1, _R2)

for RNA-seq: GTF file is needed, too
    
##### information

    REF_GENOME=mm9 # reference genome
    MAT_STRAIN=129S1 # refers to vcf column
    PAT_STRAIN=CASTEiJ # refers to vcf column
    CHROMOSOMES=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X

##### programs

###### lapels/suspenders tools

1. get_refmeta
2. vcf2mod
3. insilico
4. modmap (gene annotation)

###### additional tools

* tabix
* bowtie2/tophat

###### custom-made scripts

* changeChrNamesInMODFile.awk
* samFilter_v2_suspenders.py

![Overview](https://raw.githubusercontent.com/friedue/AlleleSpecific/master/images/pipelineOverview.png)

Image modified from [Huang et al., 2014](http://dx.doi.org/10.1093/database/bau057) to illustrate the different steps of our pipeline

-----------------------------------------------

Part I: Preparing pseudogenomes
---------------

###  I.1. lapels/bin/get_refmeta 

Meta-data for reference genome

      /package/lapels-1.0.5/bin/get_refmeta -o ${REF_GENOME}.meta ${REF_GENOME} ${REF_FASTA}


### I.2. generating MOD files

##### a) indexing of vcf files
   
	for vcf in ${VCF_INDELS} ${VCF_SNP}
    do
    	/package/tabix-0.2.6/bgzip -c ${vcf} > ${vcf}.gz
    	/package/tabix-0.2.6/tabix -p vcf ${vcf}.gz
    done

##### b) lapels/bin/vcf2mod

1 MOD file per VCF file (SNPs, INDELS) and genotype (maternal, paternal)
SNPs with bad quality (FI tag = 0) will be discarded
    
    /package/lapels-1.0.5/bin/vcf2mod -c ${CHROMOSOMES} -f \
		-o ${REF_GENOME}_SNPs_${genotype}.mod \
		${REF_GENOME} ${REF_GENOME}.meta \
		${genotype} ${VCF_SNP} 2> ${genotype}.vcf2mod_snps_`date +"%Y%m%d%H%M%S"`.log &

##### c) merge MOD files for SNPs and indels per genotype

    cat ${REF_GENOME}_SNPs_${genotype}.mod ${REF_GENOME}_indels_${genotype}.mod |\
	sort -k2,2n -k3,3n | uniq | awk -f changeChrNamesInMODFile.awk - > ${REF_GENOME}_indels_SNPs_${genotype}_changedChr.mod 

### I.3. lapels/bin/insilico

generating pseudogenomes

CAVE: after this step, the MOD file will be gzipped (without any indication in the file name)
    
	/package/lapels-1.0.5/bin/insilico \
		${REF_GENOME}_indels_SNPs_${genotype}_changedChr.mod \
		${REF_FASTA} -v -f \
		-o pseudogenome_${REF_GENOME}_${genotype}.fa >${genotype}.insilico_`date +"%Y%m%d%H%M%S"`.log 

---------------------------------------

Part II: Mapping
------------------

### II.1. Bowtie/Tophat: mapping to each pseudogenome

bowtie/tophat indeces

	/package/bowtie2-2.1.0/bowtie2-build pseudogenome_${REF_GENOME}_${genotype}.fa pseudogenome_${REF_GENOME}_${genotype} > ${genotype}.bowtieIndex_`date +"%Y%m%d%H%M%S"`.log 2> ${genotype}.bowtieIndex_Error_`date +"%Y%m%d%H%M%S"`.log

##### bowtie2

	(/package/bowtie2-2.1.0/bowtie2 -x pseudogenome_${REF_GENOME}_${genotype} \
		-1 ${fastq_folder}${sample}_R1.fastq.gz -2 ${fastq_folder}${sample}_R2.fastq.gz \
		-X 1000 -p 40 \
		--rg-id mpi-ie --rg CN:deep_sequencing_unit --rg PL:illumina | \
	/package/samtools/samtools view -Sb - |\
	/package/samtools/samtools sort - ${genotype}_${sample}.sorted) 2>&1 > MappingTo${genotype}_${sample}.log

    /package/samtools/samtools index ${genotype}_${sample}.sorted.bam &
  
##### RNA-seq: tophat

1. generate GTFs for S129 and CASTEiJ: **modmap**
2. transcriptome build for S129I and CASTEiJ genotypes using 1 set of FASTQs with **TopHat**
3. run TopHat on remaining samples

###### modmap

modmap expects 0-based positions, since gtf files are 1-based, I first convert them into 0-based, use modmap to change the annotation to match the individual genomes, and turn the resulting 0-based gtf file back into 1-based (additionally, I remove the negative numbers introduced by modmap for regions that fall into deletions)
    
    awk -F "\t" '{OFS="\t"; print $1,$2,$3,$4-1, $5-1, $6, $7, $8,$9}' ${REF_GENOME}.gtf > ${REF_GENOME}_0based.gtf
    
    /package/lapels-1.0.5/bin/modmap -d '\t' \
		-f ${REF_GENOME}_indels_SNPs_${genotype}_changedChr.mod	\
		${REF_GENOME}_0based.gtf ${REF_GENOME}_${genotype}_0based.gtf 1,4 1,5
    
     awk -F "\t" '{OFS="\t";print $1,$2,$3, $4<0?$4*-1:$4+1,$5<0?$5*-1:$5+1,$6,$7,$8,$9}' ${REF_GENOME}_${genotype}_0based.gtf > ${REF_GENOME}_${genotype}_1based.gtf
    

###### TopHat

 
**building the transcriptome index for 129S1 and CASTEiJ**

    STD_DEV=`sed '1d' fastq/insert_stats/insert_stats_RNA__1.txt | sed -n '3p' | awk '{print sprintf("%.0f",$2)}'`
    INNER_DIST=`sed '1d' fastq/insert_stats/insert_stats_RNA_WT_1.txt | sed -n '3p' | awk '{print sprintf("%.0f",$1)}'`
      
    /package/tophat-2.0.8b.Linux_x86_64/tophat2 -G ${REF_GENOME}_${genotype}_1based.gtf \
		--transcriptome-index transcriptome_data/${REF_GENOME}_${genotype}_transcriptomeIndex \
    	--no-coverage-search --library-type fr-firststrand --mate-inner-dist ${INNER_DIST} \
    	--mate-std-dev ${STD_DEV} \
		pseudogenome_${REF_GENOME}_${genotype} \
    	-o RNA_WT_${genotype} \
	    SAMPLE_R1.fastq.gz \
    	SAMPLE_R2.fastq.gz
    
**looping over remaining files (for simplicity, I moved the fastq files)**

        for fastq in ${FASTQ_FOLDER}*_R1.fastq.gz
        
		forStats=$(basename "$fastq" .fastq.gz | sed 's/_R[1-2]//' |  awk '{print substr($0,index($0,"_RNA")+1)}')
        
		STD_DEV=`sed '1d' fastq/insert_stats/insert_stats_${forStats}.txt | sed -n '3p' | awk '{print sprintf("%.0f",$2)}'`
        
		INNER_DIST=`sed '1d' fastq/insert_stats/insert_stats_${forStats}.txt | sed -n '3p' | awk '{print sprintf("%.0f",$1)}'`
        
        sample=$(basename "$fastq" .fastq.gz | sed 's/_R[1-2]//') # extracting sample name
        
        /package/tophat-2.0.8b.Linux_x86_64/tophat2 \
			-o ${sample}_${genotype} -p 20 \
			-G ${REF_GENOME}_${genotype}_1based.gtf \
			--transcriptome-index transcriptome_data/Mus_musculus_${genotype}_transcriptomeIndex \
			--no-coverage-search --library-type fr-firststrand \
			--mate-inner-dist ${INNER_DIST} --mate-std-dev ${STD_DEV} \
			pseudogenome_${REF_GENOME}_${genotype} \
			${FASTQ_FOLDER}/${sample}_R1.fastq.gz ${FASTQ_FOLDER}/${sample}_R2.fastq.gz

**indexing**
   
		/package/samtools/samtools index ${FOLDER}/accepted_hits.bam &
  

### II.2. lapels/bin/pylapels

**translate different pseudogenome mappings back to reference genome**


    /package/lapels-1.0.5/bin/pylapels -p 50 -f \
		-o ${bam}_mappedBackTo_${REF_GENOME}.bam \
		${REF_GENOME}_indels_SNPs_${genotype}_changedChr.mod ${BAM} > ${bam}_stdout.txt 2> ${bam}_stderr.txt
    
 

### II.3 suspenders/bin/pysuspenders

* merging the parental and maternal mapping, determining best fit for each read --> 1 BAM file with reads where flags indicate maternal or paternal origin

		/package/suspenders-0.2.4/bin/pysuspenders #
			--lapels --quality -p 15 mergedAlignment_${SAMPLE}.bam #
			tophat_${SAMPLE}_${MAT_STRAIN}/accepted_hits_pylapels.bam #
			tophat_${SAMPLE}_${PAT_STRAIN}/accepted_hits_pylapels.bam

-----------------------------------------

Part III: Counting allele-specific reads
-----------------------------------------
### III BAM file filter, merge, sort, index

Suspender files should be filtered:

   * multiple alignments
   * randomly assigned reads 
   * chrM, random chromosomes (/data/projects/scripts/samFilter_v2.py)

     /package/samtools/samtools view -h -F 4 ${BAM}|\
	./samFilter_v2_suspenders.py \
		--filter_out_from_BED conspicuousEnrichments.forDelete.bed \
	    --random --chrM --multiple --mismatch --lowqual --suspendersRandom |\
      /package/samtools/samtools view -Sb - > ${out}_filtered.rnsort.bam &
      
coordinate sorting (most tools want that) + indexing

        for BAM in *rnsort*bam
        do
        out=`basename "$BAM" .bam`
        /package/samtools/samtools sort -m 4000000000 ${BAM} ${out}.sort 
        done   
        
        for BAM in *.sort.bam
        do
        /package/samtools/samtools index ${BAM} &
        done
        wait
