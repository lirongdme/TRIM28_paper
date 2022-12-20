#!/bin/bash
#	Run on strix (64X)
gunzip *.fastq.gz;
#	Step 3:	QC with trim_galore
/ddn/gs1/home/lir8/Software/ATACseq/trim_galore_v0.4.4/trim_galore --paired --stringency 2 --dont_gzip --length 20 --quality 20 --trim1 --output_dir . *.1.fastq *.2.fastq;
#	Step 4:	Map the P.E. reads
nohup bowtie -p 5 -q --phred33-quals -m 1 -v 2 --best --strata --maxins 10000 --sam /ddn/gs1/home/lir8/Software/reference/hg38 -1 *1.fq -2 *2.fq Uniquely_mapped.sam;
#	Step 5:	Count reads mapped to chrM (> 50%)
grep chrM Uniquely_mapped.sam > chrM.txt;
wc -l chrM.txt > chrM-size.txt;
#	Step 6:	Remove duplicated reads
nohup java -Xms2g -jar /ddn/gs1/home/lir8/Software/Picard/picard-tools-1.96/SortSam.jar INPUT= Uniquely_mapped.sam OUTPUT= Uniquely_mapped_Sorted.sam SORT_ORDER=coordinate;
nohup java -Xms2g -jar /ddn/gs1/home/lir8/Software/Picard/picard-tools-1.96/MarkDuplicates.jar INPUT=Uniquely_mapped_Sorted.sam OUTPUT=Uniquely_mapped_Sorted_DeDUPLICATE.bam REMOVE_DUPLICATES=TRUE METRICS_FILE=DeDUPLICATE.txt;
#	Step 7:	Remove unmapped reads
nohup samtools view -b -F 4 Uniquely_mapped_Sorted_DeDUPLICATE.bam > Uniquely_mapped_Sorted_DeDUPLICATE-Clean.bam;
nohup samtools view -h Uniquely_mapped_Sorted_DeDUPLICATE-Clean.bam > Uniquely_mapped_Sorted_DeDUPLICATE-Clean.sam;	
# 	Step 8:	Remove reads mapped to chrM:	
grep chrM Uniquely_mapped_Sorted_DeDUPLICATE-Clean.sam > chrM-clean.txt;
wc -l chrM-clean.txt > chrM-clean-size.txt;
grep -v chrM Uniquely_mapped_Sorted_DeDUPLICATE-Clean.sam > Uniquely_mapped_Sorted_DeDUPLICATE-Clean_noChrM.sam;
# 	Step 9:	Prepare SAM file based on trimmed 9bp reads.
nohup perl /ddn/gs1/shared/dirib/scripts_to_share/sam2extend.pl Uniquely_mapped_Sorted_DeDUPLICATE-Clean_noChrM.sam 9;
# 	Step 10:	Random sampling to same depth 
samtools view -Sb Uniquely_mapped_Sorted_DeDUPLICATE-Clean_noChrM.ext9.sam > Uniquely_mapped_Sorted_DeDUPLICATE-Clean_noChrM.ext9.bam;
bamToBed -i Uniquely_mapped_Sorted_DeDUPLICATE-Clean_noChrM.ext9.bam > Uniquely_mapped_Sorted_DeDUPLICATE-Clean_noChrM.ext9.bed;
wc -l Uniquely_mapped_Sorted_DeDUPLICATE-Clean_noChrM.ext9.bed > Final-usable-size.txt;

# downsize all the ATAC-seq files to the same depth
shuf Uniquely_mapped_Sorted_DeDUPLICATE-Clean_noChrM.ext9.bed | head -63112882 | sort -k1,1 -k2,2n -k3,3n > Final_read_sameDepth.bed;	
macs2 callpeak --treatment Final_read_sameDepth.bed -f BED --name Peak --nomodel --extsize 9 --qvalue 0.0001 --gsize hs --keep-dup all;	
genomeCoverageBed -i Final_read_sameDepth.bed -bg -g hg38.chrom.sizes > Uniquely_mapped_Sorted_DeDUPLICATE-Clean_noChrM.ext9.bam.bedgraph;
bedGraphToBigWig Uniquely_mapped_Sorted_DeDUPLICATE-Clean_noChrM.ext9.bam.bedgraph hg38.chrom.sizes Uniquely_mapped_Sorted_DeDUPLICATE-Clean_noChrM.ext9.bam.bw;



