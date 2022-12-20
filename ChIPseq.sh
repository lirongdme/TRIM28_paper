#ChIPseq analysis of human data
gunzip *.gz;
nohup perl /ddn/gs1/home/lir8/Software/ChIPseq/trim_and_filter_SE.pl -i *.fastq -a 1 -b 50 -m 20 -q sanger -o Trim;
nohup bowtie -p 10 -q --phred33-quals -m 1 -v 2 -S /ddn/gs1/home/lir8/Software/reference/hg38 Trim.trim_1_50.minQS_20.fastq Uniquely_mapped.sam;
nohup java -Xmx2g -jar /ddn/gs1/home/lir8/Software/Picard/picard-tools-1.96/SortSam.jar INPUT= Uniquely_mapped.sam OUTPUT= Uniquely_mapped_Sorted.sam SORT_ORDER=coordinate;
nohup java -Xmx2g -jar /ddn/gs1/home/lir8/Software/Picard/picard-tools-1.96/MarkDuplicates.jar INPUT=Uniquely_mapped_Sorted.sam OUTPUT=Uniquely_mapped_Sorted_DeDUPLICATE.bam REMOVE_DUPLICATES=TRUE METRICS_FILE=DeDUPLICATE.txt;
samtools view -h Uniquely_mapped_Sorted_DeDUPLICATE.bam > Uniquely_mapped_Sorted_DeDUPLICATE.sam;
perl /ddn/gs1/home/lir8/Software/ChIPseq/sam2extend.pl Uniquely_mapped_Sorted_DeDUPLICATE.sam 300;
nohup java -Xmx2g -jar /ddn/gs1/home/lir8/Software/Picard/picard-tools-1.96/SortSam.jar INPUT= Uniquely_mapped_Sorted_DeDUPLICATE.ext300.sam OUTPUT= Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam SORT_ORDER= coordinate;
samtools view -b Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam | genomeCoverageBed -ibam stdin -bg -g /ddn/gs1/home/lir8/Software/reference/hg38.chrom.sizes > Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam.bedgraph;
sort -k1,1 -k2,2n Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam.bedgraph > Uniquely_mapped_Sorted_DeDUPLICATE_sorted.ext300.bam.bedgraph;
nohup bedGraphToBigWig Uniquely_mapped_Sorted_DeDUPLICATE_sorted.ext300.bam.bedgraph /ddn/gs1/home/lir8/Software/reference/hg38.chrom.sizes Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam.bw;

################ Peak calling with input  #############################
nohup macs2 callpeak --treatment Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam --control /ddn/gs1/home/lir8/Chipseq/input/Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam --name Peak --nomodel --qvalue 0.0001 --gsize hs;	
bedGraphToBigWig Uniquely_mapped_Sorted_DeDUPLICATE_sorted.ext300.bam.bedgraph /ddn/gs1/home/lir8/Software/reference/hg38.chrom.sizes Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam.bw;
sort -k1,1 -k2,2n Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam.bedgraph > Uniquely_mapped_Sorted_DeDUPLICATE_sorted.ext300.bam.bedgraph
nohup macs2 callpeak --treatment CutRun_bowtie2.Sorted_DeDUPLICATE_sorted.bam --control /ddn/gs1/home/lir8/Chipseq/CUTRUN/Pilot/Igg/CutRun_bowtie2.Sorted_DeDUPLICATE_sorted.bam --name Peak --nomodel --qvalue 0.0001 --gsize hs;	


#ChIPseq analysis of mouse data
nohup perl /ddn/gs1/shared/fargod/scripts_to_share/trim_and_filter_SE.pl -i *.fastq -a 1 -b 50 -m 20 -q sanger -o Trim;
nohup bowtie -p 10 -q --phred33-quals -m 1 -v 2 -S /ddn/gs1/shared/dirib/reference_genomes/mm10/indexBowtie_mm10ordered/mm10 Trim.trim_1_50.minQS_20.fastq Uniquely_mapped.sam;
nohup java -Xmx2g -jar /ddn/gs1/home/lir8/Software/Picard/picard-tools-1.96/SortSam.jar INPUT= Uniquely_mapped.sam OUTPUT= Uniquely_mapped_Sorted.sam SORT_ORDER=coordinate;
nohup java -Xmx2g -jar /ddn/gs1/home/lir8/Software/Picard/picard-tools-1.96/MarkDuplicates.jar INPUT=Uniquely_mapped_Sorted.sam OUTPUT=Uniquely_mapped_Sorted_DeDUPLICATE.bam REMOVE_DUPLICATES=TRUE METRICS_FILE=DeDUPLICATE.txt;
samtools view -h Uniquely_mapped_Sorted_DeDUPLICATE.bam > Uniquely_mapped_Sorted_DeDUPLICATE.sam;
perl /ddn/gs1/home/lir8/Software/ChIPseq/sam2extend.pl Uniquely_mapped_Sorted_DeDUPLICATE.sam 300;
nohup java -Xmx2g -jar /ddn/gs1/home/lir8/Software/Picard/picard-tools-1.96/SortSam.jar INPUT= Uniquely_mapped_Sorted_DeDUPLICATE.ext300.sam OUTPUT= Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam SORT_ORDER= coordinate;
samtools view -b Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam | genomeCoverageBed -ibam stdin -bg -g /ddn/gs1/home/lir8/Software/reference/mm10.chr.size > Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam.bedgraph;
sort -k1,1 -k2,2n Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam.bedgraph > Uniquely_mapped_Sorted_DeDUPLICATE_sorted.ext300.bam.bedgraph;
nohup bedGraphToBigWig Uniquely_mapped_Sorted_DeDUPLICATE_sorted.ext300.bam.bedgraph /ddn/gs1/home/lir8/Software/reference/mm10.chr.size Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam.bw;

################ Peak calling with input  #############################
nohup macs2 callpeak --treatment Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam --control /ddn/gs1/home/lir8/Chipseq/mTRIM28/input/Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam --name Peak --nomodel --qvalue 0.001 --gsize mm;	
bedGraphToBigWig Uniquely_mapped_Sorted_DeDUPLICATE_sorted.ext300.bam.bedgraph /ddn/gs1/home/lir8/Software/reference/hg38.chrom.sizes Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam.bw;
sort -k1,1 -k2,2n Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam.bedgraph > Uniquely_mapped_Sorted_DeDUPLICATE_sorted.ext300.bam.bedgraph
nohup macs2 callpeak --treatment Uniquely_mapped_Sorted_DeDUPLICATE.ext300.bam --name Peak --nomodel --qvalue 0.001 --gsize mm;	


####Rstudio for peak overlap
library('class');
library(GenomicRanges);
#	Read in files
List1 <- read.delim("~/peak/List1.txt", stringsAsFactors=FALSE,header=FALSE)
List2 <- read.delim("~/peak/List2.txt", stringsAsFactors=FALSE,header=FALSE)
List3 <- read.delim("~/peak/List3.txt", stringsAsFactors=FALSE,header=FALSE)
#	Format file
List1GR <- GRanges(seqnames=List1[,1], ranges=IRanges(List1[,2],List1[,3]));
List2GR <- GRanges(seqnames=List2[,1], ranges=IRanges(List2[,2],List2[,3]));
List3GR <- GRanges(seqnames=List3[,1], ranges=IRanges(List3[,2],List3[,3]));
# 	Find Overlap peaks
Overlap <- countOverlaps(List1GR, List2GR, ignore.strand=TRUE);
write.table(data.frame(Peak=unname(Overlap)), sep = "\t", file = "~/peak/Overlap-List1+List2.txt");
Overlap <- countOverlaps(List1GR, List3GR, ignore.strand=TRUE);
write.table(data.frame(Peak=unname(Overlap)), sep = "\t", file = "~/peak/Overlap-List1+List3.txt");
Overlap <- countOverlaps(List2GR, List3GR, ignore.strand=TRUE);
write.table(data.frame(Peak=unname(Overlap)), sep = "\t", file = "~/peak/Overlap-List2+List3.txt");


####Homer motif anlaysis
findMotifsGenome.pl TRIM28.bed hg38 ./TRIM28/ -size 200 -p 5