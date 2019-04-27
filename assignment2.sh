###########################################
#  - NGS2 Course - Assignment             #
#  - Bash Script                          #
#  - April 27,2019                        #
#  - Copyright: Mohamed AboelEla          #
#  - Nile University                      #
###########################################
#!/bin/bash

### Downloading Reference Genome ###
mkdir -p ../ref/
wget -P ../ref/ http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa

### Changing Data Files Name ###
mv ../data/SRR8797509_1.part_001.part_001.fastq.gz ../data/S1_L001_R1_001.fastq.gz
mv ../data/SRR8797509_2.part_001.part_001.fastq.gz ../data/S1_L001_R2_001.fastq.gz
mv ../data/shuffled_SRR8797509_1.part_001.part_001.fastq.gz ../data/S2_L001_R1_001.fastq.gz
mv ../data/shuffled_SRR8797509_2.part_001.part_001.fastq.gz ../data/S2_L001_R2_001.fastq.gz

### Unzipping Data Files ###
cd ../data/
gunzip -k *.fastq.gz
cd -

### 01. 2-PASS Alignment ###
### a. 1st Indexing ###
GENOME_DIR="../ref/star-index/"
mkdir -p $GENOME_DIR && cd $GENOME_DIR
ln -s ../chr22_with_ERCC92.fa .
cd -
STAR --runMode genomeGenerate --genomeDir $GENOME_DIR --genomeFastaFiles ../ref/chr22_with_ERCC92.fa

### b. 1st Alignment ###
PASS_1_DIR=../GATK_results/star_res/1pass
mkdir -p $PASS_1_DIR
for R1 in ../data/*_R1_001.fastq;
do
	mkdir -p $PASS_1_DIR/$(basename $R1 _R1_001.fastq)
	cd $PASS_1_DIR/$(basename $R1 _R1_001.fastq)
	R2=$(echo $R1 | sed 's/_R1_/_R2_/')
	echo $R1 $R2
	STAR --genomeDir ../../../../ref/star-index/ --readFilesIn ../../../$R1 ../../../$R2
	cd -
done

### c. 2nd Indexing ###
for S in {1..2};
do
	GENOME_DIR="../ref/star-index-2/S"$S
	mkdir -p $GENOME_DIR && cd $GENOME_DIR
	ln -s ../../chr22_with_ERCC92.fa .\
	cd -
	STAR --runMode genomeGenerate --genomeDir $GENOME_DIR --genomeFastaFiles ../ref/chr22_with_ERCC92.fa --sjdbFileChrStartEnd ../GATK_results/star_res/1pass/S${S}_L001/SJ.out.tab --sjdbOverhang 75 --runThreadN 4
 
done

### d. 2nd Alignment ###
PASS_2_DIR=../GATK_results/star_res/2pass
mkdir -p $PASS_2_DIR
for R1 in ../data/*_R1_001.fastq;
do
	mkdir -p $PASS_2_DIR/$(basename $R1 _R1_001.fastq)
	cd $PASS_2_DIR/$(basename $R1 _R1_001.fastq)
	R2=$(echo $R1 | sed 's/_R1_/_R2_/')
	echo $R1 $R2
	echo ../../../../ref/star-index-2/$(basename $R1 | cut -d"_" -f1)
	STAR --genomeDir ../../../../ref/star-index-2/$(basename $R1 | cut -d"_" -f1) --readFilesIn ../../../$R1 ../../../$R2
	cd -
done

### 02. Adding Read Groups, Sortting, Marking Duplicates, and Creating Index ###
mkdir -p ../GATK_results/picard_res
for i in S1_L001 S2_L001;
do
	SM=$(basename $i | cut -d"_" -f1)      
	LB=$i                             
	PL="Illumina"                     
	RGID=$(cat ../data/${i}_R1_001.fastq | head -n1 | sed 's/ /_/g' | cut -d "_" -f1) 
	PU=$RGID.$LB                                                     
        echo -e "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU"
	
	# Adding reading group information and sorting
	picard AddOrReplaceReadGroups I=../GATK_results/star_res/2pass/$i/Aligned.out.sam O=../GATK_results/picard_res/${SM}_rg_added_sorted.bam SO=coordinate RGID=$RGID RGLB=$LB RGPL=$PL RGPU=$PU RGSM=$SM 
	
	# Marking duplicates and creating index
	picard MarkDuplicates I=../GATK_results/picard_res/${SM}_rg_added_sorted.bam O=../GATK_results/picard_res/${SM}_dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=../GATK_results/picard_res/${SM}_output.metrics 
done

### 03. Splitting 'N'Trim and reassign mapping qualities ###
#Creating dictionary and fastq index files for reference file
cd ../ref
samtools fqidx chr22_with_ERCC92.fa 
gatk CreateSequenceDictionary -R chr22_with_ERCC92.fa -O chr22_with_ERCC92.dict
cd -

#Splitting reads into exon and hard-clipping sequences overhanging into the intronic regions
mkdir -p ../GATK_results/gatk_res
for i in S1 S2;
do
	gatk SplitNCigarReads -R ../ref/chr22_with_ERCC92.fa -I ../GATK_results/picard_res/${i}_dedupped.bam -O ../GATK_results/gatk_res/${i}_split.bam
done 

### 04. Base Recalibration ###
#Downloading known variant for Chr22
wget -P ../ref ftp://ftp.ensembl.org/pub/release-96/variation/vcf/homo_sapiens/homo_sapiens-chr22.vcf.gz

#Correctting chromosome name
cd ../ref
gunzip -k homo_sapiens-chr22.vcf.gz 
grep "^22" homo_sapiens-chr22.vcf | sed 's/^22/chr22/' > chr22.vcf
gatk IndexFeatureFile -F chr22.vcf
cd -

#Recalibrating bases using GATK BaseRecalibrator
for sample in S1 S2;
do
	name=${sample%.dedup.bam}
	gatk --java-options "-Xmx2G" BaseRecalibrator \
	-R ../refdog_chr5.fa -I $sample --known-sites ../ref/chr22.vcf \
	-O $name.report
	gatk --java-options "-Xmx2G" ApplyBQSR \
	-R dog_chr5.fa -I $sample -bqsr $name.report \
	-O $name.bqsr.bam --add-output-sam-program-record --emit-original-quals
done

### 05. Variant Calling ###

gatk HaplotypeCaller -R ../ref/chr22_with_ERCC92.fa -I ../GATK_results/picard_res/S1_dedupped.bam -O ../GATK_results/gatk_res/S1_dedupped.vcf
gatk HaplotypeCaller -R ../ref/chr22_with_ERCC92.fa -I ../GATK_results/picard_res/S2_dedupped.bam -O ../GATK_results/gatk_res/S2_dedupped.vcf

### 06. Variant Filteration ###
gatk VariantFiltration -R ../ref/chr22_with_ERCC92.fa -V ../GATK_results/gatk_res/S1_dedupped.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ../GATK_results/gatk_res/S1.vcf
gatk VariantFiltration -R ../ref/chr22_with_ERCC92.fa -V ../GATK_results/gatk_res/S1_dedupped.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ../GATK_results/gatk_res/S1.vcf

### THE END ###
