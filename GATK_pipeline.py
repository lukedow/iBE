# MAP WITH STAR 
STAR \
--runThreadN 6 \
--outFilterMismatchNmax 5 \
--genomeDir /path/to/folder/GRCm39/STAR103 \
--readFilesCommand zcat \
--readFilesIn [comma separate FASTQs] \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ./[sample_name]
	   
#MARK DUPLICATES	
for i in [list of samples]
do
java -jar ~/.conda/envs/kal_env/share/picard-2.25.2-0/picard.jar MarkDuplicates -I $i\Aligned.sortedByCoord.out.bam -O $i\_marked_duplicates.bam -M $i\_marked_dup_metrics.txt
done	

# INDEX BAMs
for i in [list of samples]
do
java -jar ~/.conda/envs/kal_env/share/picard-2.25.2-0/picard.jar BuildBamIndex -I $i\_marked_duplicates.bam
done	

# SPLIT AND CIGAR READS
for i in [list of samples]
do	
java -jar ~/.conda/envs/kal_env/share/gatk4-4.2.0.0-1/gatk-package-4.2.0.0-local.jar SplitNCigarReads \
	-R /path/to/folder/mm10.dna.primary_assembly.fa \
	-I $i\_marked_duplicates.bam \
	-O $i\_split.bam
done

# ADD READ GROUP (required for Base Recalibration)
for i in [list of samples]
do
java -jar ~/.conda/envs/kal_env/share/picard-2.25.2-0/picard.jar AddOrReplaceReadGroups \
       -I $i\_split.bam \
       -O $i\_splitout.bam \
       -RGID 4 \
       -RGLB lib1 \
       -RGPL ILLUMINA \
       -RGPU unit1 \
       -RGSM 20
done

# Downloaded variant file for known SNPS (from NCBI dbsnp: ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/)
# INDEX .vcf file
java -jar ~/.conda/envs/kal_env/share/gatk4-4.2.0.0-1/gatk-package-4.2.0.0-local.jar IndexFeatureFile -F mm10_All.vcf

# BASE RECALIBRATION
for i in [list of samples]
do	
java -jar ~/.conda/envs/kal_env/share/gatk4-4.2.0.0-1/gatk-package-4.2.0.0-local.jar BaseRecalibrator \
   -I $i\_splitout.bam \
   -R /path/to/folder/mm10.dna.primary_assembly.fa \
   --known-sites /path/to/folder/dbSNP/mm10_All.vcf \
   -O $i\_recal_data.table
done

# APPLY RECALIBRATION
for i in [list of samples]
do	
java -jar ~/.conda/envs/kal_env/share/gatk4-4.2.0.0-1/gatk-package-4.2.0.0-local.jar ApplyBQSR \
   -R /path/to/folder/mm10.dna.primary_assembly.fa \
   -I $i\_splitout.bam \
   --bqsr-recal-file $i\_recal_data.table \
   -O $i\_recal.bam
done

gatk AnalyzeCovariates \
     -bqsr recal1.table \
     -before recal2.table \
     -after recal3.table \
     -plots AnalyzeCovariates.pdf

# CALL VARIANTS
for i in [list of samples] 
do		
java -jar ~/.conda/envs/kal_env/share/gatk4-4.2.0.0-1/gatk-package-4.2.0.0-local.jar HaplotypeCaller  \
   -R /path/to/folder/mm10.dna.primary_assembly.fa \
   -I $i\_recal.bam \
   -O $i\_output.vcf.gz \
   -bamout $i\_bamout.bam
done


# FILTER VARIANTS (add filter tag to variants with less than 10 reads covering)
for i in [list of samples]
do
java -jar ~/.conda/envs/kal_env/share/gatk4-4.2.0.0-1/gatk-package-4.2.0.0-local.jar VariantFiltration \
   -R /path/to/folder/mm10.dna.primary_assembly.fa \
   -V $i\_output.vcf.gz \
   -O $i\_filtered.vcf.gz \
   --filter-name "depth" \
   --filter-expression "DP < 10"
done

# OUTPUT VCF FILES TO TABLEs
for i in [list of samples]
do
java -jar ~/.conda/envs/kal_env/share/gatk4-4.2.0.0-1/gatk-package-4.2.0.0-local.jar VariantsToTable \
   -V $i\_filtered.vcf.gz \
   -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F DP -F QD -GF AD \
   -SMA \
   -O $i\_filtered_variants.tab
done

