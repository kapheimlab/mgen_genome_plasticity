# Halictid genome resequencing project - _Megalopta genalis_
## SNP calling pipeline
## Karen M. Kapheim
## began 28 April 2016

scripts are in: `/uufs/chpc.utah.edu/common/home/u6001661/scripts`
(now moved to `/uufs/chpc.utah.edu/common/home/kapheim-group1/scripts`)

log files are in: `/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen/jobreports`

Following GATK best practices (https://software.broadinstitute.org/gatk/best-practices/).

## Pre-processing Part I

#### 1. Ran all of the following in a single script:
  * Quality trimming (sickle)
  * Make UBAM
  * Mark adapters
  * Align with bwa_mem
  * Mark duplicates

Ran separate script for each sample.
example script: `GATKpreprocessing_2612_06.slurm`
log files: `/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen/jobreports/preprocessing`


```bash
#!/bin/bash
#
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=kapheim
#SBATCH --partition=kingspeak
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load sickle
module load picard/2.1.1
module load bwa/0.7.10
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
FILE1=2612_06_CAACCTAACA_L006_R1_001.fastq
FILE2=2612_06_CAACCTAACA_L006_R2_001.fastq
SAMPLE=2612_06
LANE=L006
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/test
#
#INDEX GENOME
#bwa index $GENOME
#java -jar /uufs/chpc.utah.edu/sys/pkg/picard/1.55/CreateSequenceDictionary.jar R=$GENOME O=Mgen.soapdenovo.dict
#
cd $DATADIR
#
#QUALITY TRIMMING
sickle pe -f $FILE1 -r $FILE2 -t sanger -o ${SAMPLE}_R1.qc.fastq -p ${SAMPLE}_R2.qc.fastq -s ${SAMPLE}_singles.qc.fastq
#
#MAKE UBAM FILE
java -Xmx8G -jar $PICARDPATH/picard.jar FastqToSam FASTQ=${SAMPLE}_R1.qc.fastq FASTQ2=${SAMPLE}_R2.qc.fastq OUTPUT=${SAMPLE}_fastq2bam.bam READ_GROUP_NAME=$LANE SAMPLE_NAME=$SAMPLE LIBRARY_NAME=hyper_kapa PLATFORM=Illumina_HiSeq4000
#
#MARK ADAPTERS
java -Xmx8G -jar $PICARDPATH/picard.jar MarkIlluminaAdapters I=${SAMPLE}_fastq2bam.bam O=${SAMPLE}_markadapters.bam M=${SAMPLE}_markadapters_metrics.txt TMP_DIR=$SCRATCH
#
#ALIGNMENT
set -o pipefail
java -Xmx8G -jar $PICARDPATH/picard.jar SamToFastq I=${SAMPLE}_markadapters.bam FASTQ=${SAMPLE}_sam2fq.fastq CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true TMP_DIR=$SCRATCH
bwa mem -M -t 16 -p $GENOME ${SAMPLE}_sam2fq.fastq > ${SAMPLE}_bwa_mem.sam
java -Xmx16G -jar $PICARDPATH/picard.jar MergeBamAlignment R=$GENOME UNMAPPED_BAM=${SAMPLE}_fastq2bam.bam ALIGNED_BAM=${SAMPLE}_bwa_mem.sam O=${SAMPLE}_mergebamaln.bam CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true \
MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS TMP_DIR=$SCRATCH
#
#MARK DUPS
java -Xmx32G -jar $PICARDPATH/picard.jar MarkDuplicatesWithMateCigar INPUT=${SAMPLE}_mergebamaln.bam OUTPUT=${SAMPLE}_mdups_MC.bam METRICS_FILE=${SAMPLE}_mdups_mc_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=$SCRATCH
#
echo "complete"
[u6001661@kingspeak1 scripts]$ cat GATKpreprocessing_2612_06.slurm
#!/bin/bash
#
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=kapheim
#SBATCH --partition=kingspeak
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load sickle
module load picard/2.1.1
module load bwa/0.7.10
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
FILE1=2612_06_CAACCTAACA_L006_R1_001.fastq
FILE2=2612_06_CAACCTAACA_L006_R2_001.fastq
SAMPLE=2612_06
LANE=L006
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/test
#
#INDEX GENOME
#bwa index $GENOME
#java -jar /uufs/chpc.utah.edu/sys/pkg/picard/1.55/CreateSequenceDictionary.jar R=$GENOME O=Mgen.soapdenovo.dict
#
cd $DATADIR
#
#QUALITY TRIMMING
sickle pe -f $FILE1 -r $FILE2 -t sanger -o ${SAMPLE}_R1.qc.fastq -p ${SAMPLE}_R2.qc.fastq -s ${SAMPLE}_singles.qc.fastq
#
#MAKE UBAM FILE
java -Xmx8G -jar $PICARDPATH/picard.jar FastqToSam FASTQ=${SAMPLE}_R1.qc.fastq FASTQ2=${SAMPLE}_R2.qc.fastq OUTPUT=${SAMPLE}_fastq2bam.bam READ_GROUP_NAME=$LANE SAMPLE_NAME=$SAMPLE LIBRARY_NAME=hyper_kapa PLATFORM=Illumina_HiSeq4000
#
#MARK ADAPTERS
java -Xmx8G -jar $PICARDPATH/picard.jar MarkIlluminaAdapters I=${SAMPLE}_fastq2bam.bam O=${SAMPLE}_markadapters.bam M=${SAMPLE}_markadapters_metrics.txt TMP_DIR=$SCRATCH
#
#ALIGNMENT
set -o pipefail
java -Xmx8G -jar $PICARDPATH/picard.jar SamToFastq I=${SAMPLE}_markadapters.bam FASTQ=${SAMPLE}_sam2fq.fastq CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true TMP_DIR=$SCRATCH
bwa mem -M -t 16 -p $GENOME ${SAMPLE}_sam2fq.fastq > ${SAMPLE}_bwa_mem.sam
java -Xmx16G -jar $PICARDPATH/picard.jar MergeBamAlignment R=$GENOME UNMAPPED_BAM=${SAMPLE}_fastq2bam.bam ALIGNED_BAM=${SAMPLE}_bwa_mem.sam O=${SAMPLE}_mergebamaln.bam CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true \
MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS TMP_DIR=$SCRATCH
#
#MARK DUPS
java -Xmx32G -jar $PICARDPATH/picard.jar MarkDuplicatesWithMateCigar INPUT=${SAMPLE}_mergebamaln.bam OUTPUT=${SAMPLE}_mdups_MC.bam METRICS_FILE=${SAMPLE}_mdups_mc_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=$SCRATCH
#
echo "complete"
```

#### 2. Ran stats

  * Sequencing depth: `seqdepth.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=kapheim
#SBATCH --partition=kingspeak
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load samtools
#
#SET VARS
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
cd $DATADIR
#
#CALCULATE SEQ DEPTH
for file in *_mergebamaln.bam; do  echo "$file"; samtools depth $file  |  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}'; done
#
echo "complete"
```

  * Flagstats run on `SAMPLE_mergebamaln.bam`: `flgstats.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=kapheim
#SBATCH --partition=kingspeak
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load samtools
#
#SET VARS
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
cd $DATADIR
#
#CALCULATE SEQ DEPTH
for file in *_mergebamaln.bam; do  echo "$file"; samtools depth $file  |  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}'; done
#
echo "complete"
[u6001661@kingspeak1 scripts]$ cat flgstats.slurm
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=kapheim
#SBATCH --partition=kingspeak
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load samtools
#
#SET VARS
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
cd $DATADIR
#
#CALCULATE SEQ DEPTH
for file in *_mergebamaln.bam; do  echo "$file"; samtools flagstat $file > ${file}.flagstats.txt; done
#
echo "complete"
```

jobreport: `slurm-1157641.out-kp160`

  * samtools on 'SAMPLE_mergebamaln.bam': `samtools_uniqmaps.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=kapheim
#SBATCH --partition=kingspeak
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load samtools
#
#SET VARS
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
cd $DATADIR
#
#CALCULATE SEQ DEPTH
for file in *_mergebamaln.bam; do  echo "$file"; samtools flagstat $file > ${file}.flagstats.txt; done
#
echo "complete"
[u6001661@kingspeak1 scripts]$ cat samtools_uniqmaps.slurm
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim
#SBATCH --partition=kingspeak
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load samtools
#
#SET VARS
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
cd $DATADIR
#
#CALCULATE SEQ DEPTH
for file in *_mergebamaln.bam; do  echo "$file";  samtools view $file | cut -f1 | sort | uniq | wc -l > ${file}.uniq_mapped_reads.txt; done
#
echo "complete"
```

job report: `slurm-1157641.out-kp160`

#### 3. INDEL identification & realignment

Ran  separate scripts for each individual sample.
Skipped BSQR step with the intention of calling high-confidence variants and then using these as a training set.

example script: `GATKpreprocessing_part2_2612_06.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SAMPLE=2612_06
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "begin indel realignment & base quality recalibration starting with the marked duplicates for sample ${SAMPLE}"
#IDENTIFY INDELS
java -Xmx32G -jar $GATKPATH/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $GENOME -I ${SAMPLE}_mdups_MC.bam -o ${SAMPLE}_mdups_MC_forInDelRealn.intervals
#REALIGN
java -Xmx32G -jar $GATKPATH/GenomeAnalysisTK.jar -T IndelRealigner -R $GENOME -I ${SAMPLE}_mdups_MC.bam -targetIntervals ${SAMPLE}_mdups_MC_forInDelRealn.intervals -o ${SAMPLE}_mdups_MC_realn.bam
#
echo "complete"
```

log files:
`slurm-1180762.err-kp292`
`slurm-41210.err-lp013`
`slurm-41211.err-lp010`
`slurm-41212.err-lp012`
`slurm-41213.err-lp009`
`slurm-41214.err-lp015`
`slurm-41215.err-lp016`
`slurm-41216.err-lp011`
`slurm-41217.err-lp014`
`slurm-41219.err-lp013`
`slurm-41220.err-lp014`
`slurm-41221.err-lp010`
`slurm-41222.err-lp015`
`slurm-41223.err-lp012`
`slurm-41224.err-lp016`
`slurm-41225.err-lp009`
`slurm-41226.err-lp011`
`slurm-41227.err-lp013`
`slurm-41228.err-lp016`
`slurm-41229.err-lp010`

#### 4. Run stats
Should rerun stats on `SAMPLE_mdups_MC_realn.bam` files in case these last preprocessing steps made a difference, but will wait until after BSQR is completed.

 ## Variant Calling - round 1
The goal here is to get a high confidence SNP set to use in BQSR.

#### 1. Haplotype caller on each individual sample

This timed out several times on lonepeak, so ran it as a bundle on my kingspeak node.

the file: `multiprog_GATK_step1.slurm`
calls the config file: `GATK_step1.conf`

example script: `GATK_haplocall_step1_2612_06.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SAMPLE=2612_06
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "begin variant calling for sample ${SAMPLE}"
echo "from GATK best practices article 3893: 'Note that versions older than 3.4 require passing the options --variant_index_type LINEAR --variant_index_parameter 128000 to set the correct index strategy for the output gVCF'"
echo "this is step 1 on best practices in website 'http://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode'"
#
#HAPLOTYPE CALLER
java -Xmx32G -jar $GATKPATH/GenomeAnalysisTK.jar -T HaplotypeCaller -R $GENOME -I ${SAMPLE}_mdups_MC_realn.bam --variant_index_type LINEAR --variant_index_parameter 128000 -ERC GVCF -o ${SAMPLE}_rawvars_v1.g.vcf
#
echo "complete"
```

log files:
`slurm-1236605.err-kp292`
`slurm-1237584.err-kp292`
`slurm-1238798.err-kp292`
`slurm-1188900.err-kp292`
`slurm-42166.err-lp013`
`slurm-42168.err-lp011`
`slurm-42169.err-lp012`
`slurm-42171.err-lp014`
`slurm-1251522.err-kp292`

#### 2. Joint genotyping for just the _M. genalis_ females
Using GenotypeGVCFs
One script for all 18 females: `GATK_haplocall_step2.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SUFFIX=_rawvars_v1.g.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "this is second step of first round of variant calling; including all mgen females"
echo "from GATK best practices article 3893: 'Note that versions older than 3.4 require passing the options --variant_index_type LINEAR --variant_index_parameter 128000 to set the correct index strategy for the output gVCF'"
echo "this is step 2 on best practices in website 'http://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode'"
#
#HAPLOTYPE CALLER
java -Xmx256G -jar $GATKPATH/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $GENOME \
-V 2612_06$SUFFIX \
-V 2622_6$SUFFIX \
-V C12_4$SUFFIX \
-V C14_4$SUFFIX \
-V C1_9$SUFFIX \
-V C20_4$SUFFIX \
-V C33_6$SUFFIX \
-V C36_6$SUFFIX \
-V C43_2$SUFFIX \
-V C49_9$SUFFIX \
-V C73_7$SUFFIX \
-V C8_1$SUFFIX \
-V GT118_2$SUFFIX \
-V GT119_4$SUFFIX \
-V GT253_4$SUFFIX \
-V GT27_8$SUFFIX \
-V GT51_2$SUFFIX \
-V GT83_6$SUFFIX \
-nt 4 \
-log haplocall_step2_v1.log \
-o haplocall_v1.vcf
#
echo "complete"
```

job report: `slurm-1283084.err-kp292`

#### 3. Genotype just the _M. genalis_ male
Using GenotypeGVCFs
Haplotype caller was set to ploidy level = 2n, so we will use the genotyping data to identify heterozygous SNPs.
Since males should not have any heterozygous SNPs, these will serve as a filtering tool for low-confidence or spurious SNPs.

script: `GATK_haplocall_step2_male2N.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SUFFIX=_rawvars_v1.g.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "this is second step of first round of variant calling; for the male Mgen only"
echo "from GATK best practices article 3893: 'Note that versions older than 3.4 require passing the options --variant_index_type LINEAR --variant_index_parameter 128000 to set the correct index strategy for the output gVCF'"
echo "this is step 2 on best practices in website 'http://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode'"
#
#HAPLOTYPE CALLER
java -Xmx256G -jar $GATKPATH/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $GENOME \
-V GT118_2_1$SUFFIX \
-nt 4 \
-log haplocall_step2_v1.log \
-o haplocall_male2N_v1.vcf
#
echo "complete"
```

job report: `slurm-1285681.err-kp292`

#### 4. Genotype just the _M. centralis_ female

Using GenotypeGVCFs
script: `GATK_haplocall_step2_MC1.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SUFFIX=_rawvars_v1.g.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "this is second step of first round of variant calling; for the M. centralis female"
echo "from GATK best practices article 3893: 'Note that versions older than 3.4 require passing the options --variant_index_type LINEAR --variant_index_parameter 128000 to set the correct index strategy for the output gVCF'"
echo "this is step 2 on best practices in website 'http://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode'"
#
#HAPLOTYPE CALLER
java -Xmx256G -jar $GATKPATH/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $GENOME \
-V MC1$SUFFIX \
-nt 4 \
-log MC1_haplocall_step2_v1.log \
-o haplocall_MC1_v1.vcf
#
echo "complete"
```

job report: `slurm-1868753.err-kp292`

## Variant Filtering - round 1
Goal here is to get a high quality SNP set to use in BQSR

#### 1. Extract SNPs
Extract SNPs from `haplocall_v1.vcf`
script for _M. genalis_ females: `GATK_extractSNPs_v1.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=owner-guest
#SBATCH --partition=lonepeak-guest
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
VARS=haplocall_v1.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "begin hard filtering - extract SNPs"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#EXTRACT SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $GENOME -V $VARS -selectType SNP -o ${VARS}_rawSNPs.vcf
#
echo "complete"
```
job report: `slurm-47958.err-lp070`

script for _M. centralis_ female: `GATK_extractSNPs_v1_MC1.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
VARS=haplocall_MC1_v1.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "begin hard filtering - extract SNPs"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#EXTRACT SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $GENOME -V $VARS -selectType SNP -o ${VARS}_rawSNPs.vcf
#
echo "complete"
```

job report: `slurm-1874106.err-kp292`

#### 2. Filter based on GATK generic recommendations:
https://www.broadinstitute.org/gatk/guide/article?id=3225
https://www.broadinstitute.org/gatk/guide/article?id=6925
https://software.broadinstitute.org/gatk/guide/article?id=2806

###### Tag SNPs to be filtered

Script for _M. genalis_ females: `GATK_filterSNPs_r1_v1.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=owner-guest
#SBATCH --partition=lonepeak-guest
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
VARS=haplocall_v1.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "hard filtering: apply filter to extracted SNPs"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#FILTER SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T VariantFiltration -R $GENOME -V ${VARS}_rawSNPs.vcf \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filterName "GATK_snp_filter" \
-o filtered_snps_v1.vcf
#
echo "complete"
```
NOTE: The output file from this script was overwritten 28nov2016 when processing the _M. centralis_ sample.
The script will need to be rerun in order to re-generate the file.
NOTE: Repeatedly get an error that MQRankSum is an undefined variable. So this part was not included in the filtering.

job report: `slurm-47960.err-lp070`

Script for _M. centralis_ female: `GATK_filterSNPs_r1_v1_MC1.slurm`

  ```bash
  #!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
VARS=haplocall_MC1_v1.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "hard filtering: apply filter to extracted SNPs"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#FILTER SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T VariantFiltration -R $GENOME -V ${VARS}_rawSNPs.vcf \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filterName "GATK_snp_filter" \
-o MC1_filtered_snps_v1.vcf
#
echo "complete"
```

job report: `slurm-1874200.err-kp292`

###### Exclude tagged SNPs
script for _M. genalis_ females: `GATK_excl_filterSNPs_r1_v1.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=owner-guest
#SBATCH --partition=lonepeak-guest
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SNPs=filtered_snps_v1.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "hard filtering: extract SNPs filered after round 1"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#FILTER SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $GENOME -V $SNPs \
--excludeFiltered \
-o filtered_excl_snps_r1_v1.vcf
#
echo "complete"
```

job report: `slurm-48232.err-lp142`

Script for _M. centralis_ female: `GATK_excl_filterSNPs_r1_v1_MC1.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SNPs=MC1_filtered_snps_v1.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "hard filtering: extract SNPs filered after round 1"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#FILTER SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $GENOME -V $SNPs \
--excludeFiltered \
-o MC1_filtered_excl_snps_r1_v1.vcf
#
echo "complete"
```

job report: `slurm-1874271.err-kp292`

#### 3. Filter based on allele type
Keep only bi-allelic SNPs

Script for _M. genalis_ females: `GATK_excl_filterSNPs_r2_v1.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=owner-guest
#SBATCH --partition=lonepeak-guest
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SNPs=filtered_excl_snps_r1_v1.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "hard filtering: keep only biallelic snps"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#FILTER SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $GENOME -V $SNPs \
--restrictAllelesTo BIALLELIC \
-o filtered_excl_snps_r2_v1.vcf
#
echo "complete"
```

job report: `slurm-48240.err-lp066`

Script for _M. centralis_ female: `GATK_excl_filterSNPs_r2_v1_MC1.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SNPs=MC1_filtered_excl_snps_r1_v1.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "hard filtering: keep only biallelic snps"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#FILTER SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $GENOME -V $SNPs \
--restrictAllelesTo BIALLELIC \
-o MC1_filtered_excl_snps_r2_v1.vcf
#
echo "complete"
```

job report: `slurm-1874288.err-kp292`

#### 4. Filter based on male heterozygosity

###### Extract SNPs
script: `GATK_extractSNPs_male2N.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=owner-guest
#SBATCH --partition=lonepeak-guest
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
VARS=haplocall_male2N_v1.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "begin filtering male heterozygous snps - extract SNPs"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#EXTRACT SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $GENOME -V $VARS -selectType SNP -o ${VARS}_rawSNPs.vcf
#
echo "complete"
```

job report: `slurm-48236.err-lp142`

###### Extract heterozygous SNPs
Pull out list of SNPs that were heterozygous in the male 2N genotyping output

script:   `GATK_filterSNPs_male2NHet.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=owner-guest
#SBATCH --partition=lonepeak-guest
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
VARS=haplocall_male2N_v1.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "hard filtering: apply filter to extracted SNPs"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#FILTER SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $GENOME -V ${VARS}_rawSNPs.vcf \
-select 'vc.getGenotype("GT118_2_1").isHet()' \
-o filtered_snps_male2NHet.vcf
#
echo "complete"
```

job report: `slurm-48281.err-lp106`

###### Exclude male heterozygous SNPs from female SNP set
This is only for _M. genalis_ females, because we do not have any male data for _M. centralis_.

script: `GATK_excl_filterSNPs_r3_v1.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=owner-guest
#SBATCH --partition=lonepeak-guest
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SNPs=filtered_excl_snps_r2_v1.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "hard filtering: get rid of snps that are heterozygous in males"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#FILTER SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $GENOME -V $SNPs \
--discordance filtered_snps_male2NHet.vcf \
-o filtered_excl_snps_r3_v1.vcf
#
echo "complete"
```

job report: `slurm-48282.err-lp121`

#### 5. Exclude variants with missing genotypes
We are still just daeling with the _M. genalis_ females here, because there is only 1 _M. centralis_ female.

script: `GATK_excl_filterSNPs_r4_v1.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim
#SBATCH --partition=kingspeak
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SNPs=filtered_excl_snps_r3_v1.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "hard filtering: get rid of SNPs that are missing genotypes in any samples"
echo "using VCFtools"
#
#FILTER SNPs
/uufs/chpc.utah.edu/sys/pkg/vcftools/vcftools_0.1.11/bin/vcftools --vcf $SNPs --max-missing-count 0 --recode --out filtered_excl_snps_r4_v1.vcf
#
echo "complete"
```

job report: `slurm-1348093.out-kp026`

We are now left with 18 individuals and 6,656,426 SNPs.
We began this step with 9,958,669 SNPs.

## Pre-processing Part II (BQSR)

Following the BQSR tutorial (https://software.broadinstitute.org/gatk/documentation/article?id=2801)
Using the high-confidence, stringent SNP set that we just made as a training set.
Start with bam files after duplicate marking and indel realignment.
Run a separate script for each sample.

example script: `GATKpreprocessing_BQSRp1_2612_06.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --account=kapheim
#SBATCH --partition=kingspeak
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
GOODSNPS=filtered_excl_snps_r4_v1.vcf
SAMPLE=2612_06
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "begin base quality recalibration starting with the samples after dups are marked and indels are realigned for sample ${SAMPLE}"
#produce covariation data tables
java -Xmx32G -jar $GATKPATH/GenomeAnalysisTK.jar -T BaseRecalibrator -R $GENOME -I ${SAMPLE}_mdups_MC_realn.bam -knownSites $GOODSNPS -o ${SAMPLE}_recal_data.table
#analyze remaining covariation after recalibration
java -Xmx32G -jar $GATKPATH/GenomeAnalysisTK.jar -T BaseRecalibrator -R $GENOME -I ${SAMPLE}_mdups_MC_realn.bam -knownSites $GOODSNPS -BQSR ${SAMPLE}_recal_data.table  -o ${SAMPLE}_post_recal_data.table
#make before/after plots
java -Xmx32G -jar $GATKPATH/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $GENOME -before ${SAMPLE}_recal_data.table -after ${SAMPLE}_post_recal_data.table -plots ${SAMPLE}_recalibration_plots.pdf
#apply the recalibration to the data
java -Xmx32G -jar $GATKPATH/GenomeAnalysisTK.jar -T PrintReads -R $GENOME -I ${SAMPLE}_mdups_MC_realn.bam -BQSR ${SAMPLE}_recal_data.table -o ${SAMPLE}_mdups_MC_realn_recal.bam
#
echo "complete"
```

Job reports:
`slurm-1351067.err-kp002`
`slurm-1351068.err-kp029`
`slurm-1351069.err-kp030`
`slurm-1351070.err-kp031`
`slurm-1351071.err-kp018`
`slurm-1351072.err-kp004`
`slurm-1351073.err-kp111`
`slurm-1351074.err-kp160`
`slurm-1351075.err-kp027`
`slurm-1351076.err-kp028`
`slurm-1351077.err-kp007`
`slurm-1351078.err-kp013`
`slurm-1351079.err-kp032`
`slurm-1351080.err-kp165`
`slurm-1351081.err-kp198`
`slurm-1351082.err-kp005`
`slurm-1351083.err-kp006`
`slurm-1351086.err-kp031`

## Variant Calling - round 2
This is the final version, now that we have completed the pre-processing.

#### 1. Haplotype caller
Run separately for each sample
example script: `GATK_haplocall_step1_v2_2612_06.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --account=kapheim
#SBATCH --partition=kingspeak
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SAMPLE=2612_06
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "second round of variant calling - after BQSR - this one is for real"
echo "begin variant calling for sample ${SAMPLE}"
echo "from GATK best practices article 3893: 'Note that versions older than 3.4 require passing the options --variant_index_type LINEAR --variant_index_parameter 128000 to set the correct index strategy for the output gVCF'"
echo "this is step 1 on best practices in website 'http://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode'"
#
#HAPLOTYPE CALLER
java -Xmx32G -jar $GATKPATH/GenomeAnalysisTK.jar -T HaplotypeCaller -R $GENOME -I ${SAMPLE}_mdups_MC_realn_recal.bam --variant_index_type LINEAR --variant_index_parameter 128000 -ERC GVCF -o ${SAMPLE}_rawvars_v2.g.vcf
#
echo "complete"
```

Haplocaller mapping stats:
| Sample | Total reads | Reads filtered out | failing DuplicateReadFilter |	failing FailsVendorQualityCheckFilter | failing HCMappingQualityFilter | failing MalformedReadFilter | failing MappingQualityUnavailableFilter | failing NotPrimaryAlignmentFilter | failing UnmappedReadFilter | jobreport |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| 2612_06	| 35818824 | 15799637	| 2996542	| 0	| 12463314 | 0 |	0 |	339781 | 0 | slurm-1353535.err-kp012
| 2622_6	| 43991580 | 16710718 | 3436258	| 0	| 12842118 | 0 | 0	| 432342 | 0 | slurm-1353536.err-kp019
| C12_4	| 41146141 | 16284147 | 3170558	| 0 |	12703479 | 0 | 0	| 410110 | 0 | slurm-1353537.err-kp023
| C14_4	| 49120729 | 21173897	| 4123624	| 0	| 16552120 | 0 | 0	| 498153 | 0 | slurm-1353538.err-kp026
| C1_9	| 36152862 | 14975485	| 2807818	| 0 |	11803520 | 0 | 0	| 364147 | 0 | slurm-1353539.err-kp014
| C20_4	| 44488320 | 3780318	| 3780318	| 0 |	15292106 | 0 | 0	| 421910 | 0 | slurm-1353540.err-kp015
| C33_6	| 39707541 | 16422895	| 3082815	| 0	| 12929530 | 0 | 0	| 410550 | 0 | slurm-1353541.err-kp016
| C36_6	| 47348088 | 18666178	| 3534707	| 0	| 14645049 | 0 | 0	| 486422 | 0 | slurm-1353542.err-kp001
| C43_2	| 42498009 | 17926717	| 3297700	| 0	| 14216320 | 0 | 0	| 412697 | 0 | slurm-1353543.err-kp008
| C49_9	| 39171564 | 15756379	| 3118202	| 0	| 12247175 | 0 | 0	| 391002 | 0 | slurm-1353544.err-kp014
| C73_7	| 46320346 | 18950283	| 3390982	| 0	| 15085436 | 0 | 0	| 473865 | 0 | slurm-1353545.err-kp016
| GT118_2	| 36005377 | 15370162	| 2913292	| 0	| 12091920 | 0 | 0	| 364950 | 0 | slurm-1353547.err-kp029
| GT119_4	| 32714835 | 13836556	| 2493191	| 0	| 11023767 | 0 | 0 |	319598 | 0 | slurm-1353548.err-kp165
| GT253_4	| 33284605 | 13835600	| 2593799	| 0	| 10917236 | 0 | 0	| 324565 | 0 | slurm-1353549.err-kp015
| GT27_8	| 39920634 | 15829169	| 2738980	| 0	| 12688009 | 0 | 0	| 402180 | 0 | slurm-1353550.err-kp023
| GT51_2	| 44097262 | 18228880	| 3331421	| 0	| 14459871 | 0 | 0	| 437588 | 0 | slurm-1353551.err-kp111
| GT83_6	| 44957338 | 17027572	| 3301459	| 0	| 13262521 | 0 | 0	| 463592 | 0 | slurm-1353552.err-kp026

#### 2. Joint genotyping
Just _M. genalis_ females
Using GenotypeGVCFs

script: `GATK_haplocall_step2_v2.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SUFFIX=_rawvars_v2.g.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "this is second step of 2nd round of variant calling; including all mgen females"
echo "from GATK best practices article 3893: 'Note that versions older than 3.4 require passing the options --variant_index_type LINEAR --variant_index_parameter 128000 to set the correct index strategy for the output gVCF'"
echo "this is step 2 on best practices in website 'http://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode'"
#
#HAPLOTYPE CALLER
java -Xmx256G -jar $GATKPATH/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $GENOME \
-V 2612_06$SUFFIX \
-V 2622_6$SUFFIX \
-V C12_4$SUFFIX \
-V C14_4$SUFFIX \
-V C1_9$SUFFIX \
-V C20_4$SUFFIX \
-V C33_6$SUFFIX \
-V C36_6$SUFFIX \
-V C43_2$SUFFIX \
-V C49_9$SUFFIX \
-V C73_7$SUFFIX \
-V C8_1$SUFFIX \
-V GT118_2$SUFFIX \
-V GT119_4$SUFFIX \
-V GT253_4$SUFFIX \
-V GT27_8$SUFFIX \
-V GT51_2$SUFFIX \
-V GT83_6$SUFFIX \
-nt 4 \
-log haplocall_step2_v2.log \
-o haplocall_v2.vcf
#
echo "complete"
```

job report: `slurm-1358908.out-kp292`

## Variant Filtering - round 2
This is the final round of variant filtering.

#### 1. Extract SNPs
Extract SNPs from `haplocall_v2.vcf`.

script: `GATK_extractSNPs_v2.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
VARS=haplocall_v2.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "begin hard filtering - extract SNPs"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#EXTRACT SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $GENOME -V $VARS -selectType SNP -o ${VARS}_rawSNPs.vcf
#
echo "complete"
```

job report: `slurm-1360193.err-kp292`

#### 2. Filter SNPs based on GATK recommendations
Filter based on GATK generic recommendations
https://www.broadinstitute.org/gatk/guide/article?id=3225
https://www.broadinstitute.org/gatk/guide/article?id=6925
https://software.broadinstitute.org/gatk/guide/article?id=2806

###### Tag SNPs to be filtered

script: `GATK_filterSNPs_r1_v2.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
VARS=haplocall_v2.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "hard filtering: apply filter to extracted SNPs"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#FILTER SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T VariantFiltration -R $GENOME -V ${VARS}_rawSNPs.vcf \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filterName "GATK_snp_filter" \
-o filtered_snps_v2.vcf
#
echo "complete"
```

job report: `slurm-1360203.err-kp292`

###### Exclude tagged SNPs

script: `GATK_excl_filterSNPs_r1_v2.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SNPs=filtered_snps_v2.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "hard filtering: extract SNPs filered after round 1"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#FILTER SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $GENOME -V $SNPs \
--excludeFiltered \
-o filtered_excl_snps_r1_v2.vcf
#
echo "complete"
```

job report: `slurm-1360220.err-kp292`

#### 3. Filter based on allele type
Keep only bi-allelic SNPs

Script: `GATK_excl_filterSNPs_r2_v2.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SNPs=filtered_excl_snps_r1_v2.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "hard filtering: keep only biallelic snps"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#FILTER SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $GENOME -V $SNPs \
--restrictAllelesTo BIALLELIC \
-o filtered_excl_snps_r2_v2.vcf
#
echo "complete"
```

job report: `slurm-1360225.err-kp292`

#### 4. Filter based on male heterozygosity

Script: `GATK_excl_filterSNPs_r3_v2.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SNPs=filtered_excl_snps_r2_v2.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "hard filtering: get rid of snps that are heterozygous in males"
echo "from GATK best practices article 2806:https://www.broadinstitute.org/gatk/guide/article?id=2806"
#
#FILTER SNPs
java -Xmx16G -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $GENOME -V $SNPs \
--discordance filtered_snps_male2NHet.vcf \
-o filtered_excl_snps_r3_v2.vcf
#
echo "complete"
```

job report: `slurm-1360230.err-kp292`

#### 5. Filter based on missingness
Exclude SNPs for which there are more than 8 missing genotypes.

Script: `GATK_excl_filterSNPs_r4_v2.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module purge
module load picard/2.1.1
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SNPs=filtered_excl_snps_r3_v2.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "hard filtering: get rid of SNPs that are missing genotypes in 8 or more samples"
echo "using VCFtools"
#
#FILTER SNPs
/uufs/chpc.utah.edu/sys/pkg/vcftools/vcftools_0.1.11/bin/vcftools --vcf $SNPs --max-missing-count 8 --recode --recode-INFO-all --out filtered_excl_snps_r4_v2.vcf
#
echo "complete"
```

job report: `slurm-1453222.out-kp292`

#### 6. Summary stats on final SNP set
Calculate # SNPs in each step wtih the command `egrep -v "^#" file.vcf | wc -l`.

**SNP filtering results**

| script | file | # snps | % filtered |
| --- | --- | :---: | :---: |
| GATK_extractSNPs_v1.slurm	| haplocall_v1.vcf_rawSNPs.vcf | 10,566,313
| GATK_filterSNPs_r1_v1.slurm	| filtered_snps_v1.vcf | 10,566,313
| GATK_excl_filterSNPs_r1_v1.slurm	| filtered_excl_snps_r1_v1.vcf	| 10,189,445 |	3.57 |
| GATK_excl_filterSNPs_r2_v1.slurm	| filtered_excl_snps_r2_v1.vcf	| 10,075,428 |	4.65 |
| GATK_extractSNPs_male2N.slurm	| haplocall_male2N_v1.vcf_rawSNPs.vcf	| 2,735,619 | |
| GATK_filterSNPs_male2NHet.slurm	| filtered_snps_male2NHet.vcf	 | 165,918 |	6.07 |
| GATK_excl_filterSNPs_r3_v1.slurm	| filtered_excl_snps_r3_v1.vcf	|  9,958,669 |	5.75 |
| GATK_excl_filterSNPs_r4_v1.slurm	| filtered_excl_snps_r4_v1.vcf	|6,656,426 |	37.00 |
| GATK_extractSNPs_v2.slurm	| haplocall_v2.vcf_rawSNPs.vcf 	|  6,351,929 | |
| GATK_filterSNPs_r1_v2.slurm	| filtered_snps_v2.vcf	| 6,351,929 | |
| GATK_excl_filterSNPs_r1_v2.slurm	| filtered_excl_snps_r1_v2.vcf	|  6,167,737 | 2.90
| GATK_excl_filterSNPs_r2_v2.slurm	| filtered_excl_snps_r2_v2.vcf	| 6,133,188 | 0.54
| GATK_excl_filterSNPs_r3_v2.slurm	| filtered_excl_snps_r3_v2.vcf	| 6,032,885 | 1.58
| GATK_excl_filterSNPs_r4_v2.slurm	| filtered_excl_snps_r4_v2.vcf.recode.vcf	 | 4,674,973 | 21.38


Use vcftools to get some summary stats.

###### Calculate average depth
script: `GATK_excl_filterSNPs_r4_v2_depth.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module load picard/2.1.1
ml vcftools
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SNPs=filtered_excl_snps_r4_v2.vcf.recode.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "running stats"
echo "using VCFtools"
#
#FILTER SNPs
vcftools --vcf $SNPs --depth --out mgen_final_depth
#
echo "complete"
```
Job report: `slurm-2402588.err-kp292`

Output: `mgen_final_depth.idepth`

| INDV | N_SITES | MEAN_DEPTH |
| --- | --- | --- |
| 2612_06 | 4670507 | 4.55998
| 2622_6  | 4674930 | 8.34166
| C12_4   | 4672436 | 6.01104
| C14_4   | 4669059 | 6.25304
| C1_9    | 4669948 | 4.83978
| C20_4   | 4672398 | 5.76416
| C33_6   | 4669469 | 5.3139
| C36_6   | 4674913 | 8.43256
| C43_2   | 4668971 | 5.68056
| C49_9   | 4673029 | 5.64995
| C73_7   | 4670478 | 6.26582
| C8_1    | 4668863 | 5.18427
| GT118_2 | 4672553 | 4.85596
| GT119_4 | 4669609 | 4.4925
| GT253_4 | 4669887 | 4.55024
| GT27_8  | 4671615 | 5.66747
| GT51_2  | 4670051 | 5.924
| GT83_6  | 4674948 | 8.65683 |

###### Calculate missing SNPs per individuals
Script: `GATK_excl_filterSNPs_r4_v2_missing_indv.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module load picard/2.1.1
ml vcftools
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SNPs=filtered_excl_snps_r4_v2.vcf.recode.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "running stats"
echo "using VCFtools"
#
#FILTER SNPs
vcftools --vcf $SNPs --missing-indv --out mgen_final_missing_indv
#
echo "complete"
```
job report: `slurm-2402614.err-kp292`
Output: `mgen_final_missing_indv.imiss`

| INDV | N_DATA | N_GENOTYPES_FILTERED | N_MISS | F_MISS |
| --- | --- | :---: | :---: | :---: |
| 2612_06 | 4674973 | 0 | 399260 | 0.0854037
| 2622_6  | 4674973 | 0 | 40596  | 0.00868369
| C12_4   | 4674973 | 0 | 121514 | 0.0259924
| C14_4   | 4674973 | 0 | 394521 | 0.08439
| C1_9    | 4674973 | 0 | 462805 | 0.0989963
| C20_4   | 4674973 | 0 | 136698 | 0.0292404
| C33_6   | 4674973 | 0 | 428759 | 0.0917137
| C36_6   | 4674973 | 0 | 38049  | 0.00813887
| C43_2   | 4674973 | 0 | 525255 | 0.112355
| C49_9   | 4674973 | 0 | 130252 | 0.0278616
| C73_7   | 4674973 | 0 | 221437 | 0.0473665
| C8_1    | 4674973 | 0 | 597737 | 0.127859
| GT118_2 | 4674973 | 0 | 209982 | 0.0449162
| GT119_4 | 4674973 | 0 | 677579 | 0.144938
| GT253_4 | 4674973 | 0 | 517368 | 0.110668
| GT27_8  | 4674973 | 0 | 190722 | 0.0407964
| GT51_2  | 4674973 | 0 | 289936 | 0.0620188
| GT83_6  | 4674973 | 0 | 36154  | 0.00773352

###### Calculate missing information on a per-SNP basis
script: `GATK_excl_filterSNPs_r4_v2_missing_site.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module load picard/2.1.1
ml vcftools
#
#SET VARS
SCRATCH=/scratch/general/lustre
GENOME=/uufs/chpc.utah.edu/common/home/kapheim-group1/Mgen/genome/Mgen.soapdenovo.fa
PICARDPATH=/uufs/chpc.utah.edu/sys/installdir/picard/2.1.1
GATKPATH=/uufs/chpc.utah.edu/sys/pkg/GATK
SNPs=filtered_excl_snps_r4_v2.vcf.recode.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "running stats"
echo "using VCFtools"
#
#FILTER SNPs
vcftools --vcf $SNPs --missing-site --out mgen_final_missing_site
#
echo "complete"
```
job report: `slurm-2402637.err-kp292`
Output: `mgen_final_missing_site.lmiss`
*Very long file - only showing the top 10 rows*

| CHR |  POS |  N_DATA | N_GENOTYPE_FILTERED | N_MISS | F_MISS |
| --- | --- | --- | --- | --- | --- |
| scaffold453  |   155  |   36   |   0   |    8    |   0.222222
| scaffold453  |   200  |   36   |   0   |    6    |   0.166667
| scaffold453  |   408  |   36   |   0   |    4    |   0.111111
| scaffold453  |   412  |   36   |   0   |    8    |   0.222222
| scaffold453  |   432  |   36   |   0   |    4    |   0.111111
| scaffold1081 |   281  |   36   |   0   |    8    |   0.222222
| scaffold1081 |   443  |   36   |   0   |    4    |   0.111111
| scaffold1081 |   448  |   36   |   0   |    4    |   0.111111
| scaffold1081 |   451  |   36   |   0   |    2    |   0.0555556

These are bundled in a tarball: `mgen_snp_stats.tar.gz`

## Impute missing genotypes
Used BEAGLE for the imputation,bu t ended up not using these, because they were causinng bias in the tests of gamma.

script: `beagle_mgen.slurm`
```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module load beagle
module load jdk/1.8.0_25
#
#SET VARS
BEAGLEPATH=/uufs/chpc.utah.edu/sys/installdir/beagle/4.1/
SNPs=filtered_excl_snps_r4_v2.vcf.recode.vcf
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
OUTDIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen/beagle
#
#
cd $DATADIR
#
echo "impute missing genotypes with BEAGLE"
#
#
java -Xmx3000m -jar /uufs/chpc.utah.edu/sys/installdir/beagle/4.1/beagle.22Apr16.1cf.jar gtgl=$SNPs out=${OUTDIR}/mgen_beagle_out nthreads=16
#java –jar $BEAGLEPATH/beagle.22Apr16.1cf.jar gtgl=$SNPs out=${OUTDIR}/mgen_beagle_out nthreads=16 phase-its=10 impute-its=10
#
echo "complete"
```

job report: `slurm-1453288.out-kp292`

## Fst
Use vcftools to calculate Fst between solitary and social females over 15 kb windows.

Lamichhaney et al 2015 used VCFtools to calculated Fst in non-overlapping 15kb windows between two phenotypes. Further divide significant regions into 5kb windows. Use PLINK to calc pairwise genetic distance between indivs and generate neighbor joining trees with Phylip. Then generated phased haplotypes for the region of interest using BEAGLE and use these to generate Juckes-Cantor corrected nt distances among the alleles for the diff phenotypes with PHYLIP.

script: `VCF_Fst_mgen.slurm`

```bash
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
module load vcftools
#
#SET VARS
SNPs=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen/beagle/mgen_beagle_out.vcf
WIN=15000
DATADIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/halictid_reseq/mgen
#
#
cd $DATADIR
#
echo "calculate fst between soc and sol mgen fems in 15kb windows"
echo "using VCFtools"
#
#FILTER SNPs
vcftools --vcf $SNPs --weir-fst-pop social.txt --weir-fst-pop solitary.txt --fst-window-size ${WIN}  --out mgen_socsol_fst
#
echo "complete"
```
job report: `slurm-1453433.err-kp292`
