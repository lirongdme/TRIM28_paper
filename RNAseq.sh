#!/bin/bash
######RNAseq for human data
gunzip *.gz;
#	Step 1:	Remove low quality reads (QC < 20)
nohup perl /ddn/gs1/home/lir8/Software/RNAseq/trim_and_filter_PE.pl -1 *.1.fastq -2 *.2.fastq -a 1 -b 75 -c 1 -d 75 -m 20 -q sanger -o RNAseq;
#	Step 2:	Map reads to reference genome (tophat):	uniquely mapped to mm10
tophat --b2-sensitive --library-type fr-unstranded -g 1 -p 10 --segment-length 25 --mate-inner-dist -40 --mate-std-dev 50 -o TopHat-Out  /ddn/gs1/shared/dirib/reference_genomes/hg38 RNAseq.1.trim_1_75.minQS_20.fastq RNAseq.2.trim_1_75.minQS_20.fastq  > runningLog.txt;
#	Step 3:	Remove duplicated reads (Picard)
cd ./TopHat-Out;
nohup java -Xmx2g -jar /ddn/gs1/home/lir8/Software/Picard/picard-tools-1.96/SortSam.jar INPUT= accepted_hits.bam OUTPUT= Mapped_Sorted.bam SORT_ORDER=coordinate;
nohup java -Xmx2g -jar /ddn/gs1/home/lir8/Software/Picard/picard-tools-1.96/MarkDuplicates.jar INPUT= Mapped_Sorted.bam OUTPUT= Mapped_Sorted_DeDUPLICATE.bam REMOVE_DUPLICATES=TRUE METRICS_FILE=DeDUPLICATE.txt;
#	Generate a QC summary table
#	Step 4:	Display the mapped read coverage on UCSC Genome Browser
nohup genomeCoverageBed -ibam Mapped_Sorted_DeDUPLICATE.bam -g /ddn/gs1/home/lir8/Software/RNAseq/mm10.chr.size -bg -split > RNAseq_alignment.bedgraph;
bedGraphToBigWig RNAseq_alignment.bedgraph /ddn/gs1/home/lir8/Software/RNAseq/mm10.chr.size RNAseq_alignment.bw;

####Pari-wise comparison for DEG using Rstudio
library(edgeR);
x <- read.delim("HTSeq-Max50.txt", check.names=FALSE, stringsAsFactors=FALSE);
# siNTveh1	siNTeh2	siNTveh3	siTRIM28veh1 siTRIM28veh2 siTRIM28veh3
y <- DGEList(counts=x[,2:7],genes=x[,1]);
y <- DGEList(counts=x[,2:7],genes=x[,1]);
Cell <- factor(c("v1","v2","v3","v1","v2","v3"));
Treatment <- factor(c("NT","NT","NT","TR","TR","TR"));
data.frame(Sample=colnames(y),Cell,Treatment);
design <- model.matrix(~Cell+Treatment);
rownames(design) <- colnames(y);
design;
y <- estimateDisp(y, design, robust=TRUE);
y$common.dispersion;
plotBCV(y)
dev.off();
fit <- glmFit(y, design);
lrt <- glmLRT(fit);
topTags(lrt);
FDR <- p.adjust(lrt$table$PValue, method="BH");
sum(FDR < 0.05);
Output <- topTags(lrt,n=50000);
write.table(data.frame(Output), sep = "\t", file = "DEGs-NT-TR.txt");


######RNAseq for mouse data
gunzip *.gz;
#	Step 1:	Remove low quality reads (QC < 20)
nohup perl /ddn/gs1/home/lir8/Software/RNAseq/trim_and_filter_PE.pl -1 *.1.fastq -2 *.2.fastq -a 1 -b 75 -c 1 -d 75 -m 20 -q sanger -o RNAseq;
#	Step 2:	Map reads to reference genome (tophat):	uniquely mapped to mm10
tophat --b2-sensitive --library-type fr-unstranded -g 1 -p 10 --segment-length 25 --mate-inner-dist -40 --mate-std-dev 50 -o TopHat-Out  /ddn/gs1/shared/dirib/reference_genomes/mm10/indexBowtie2_mm10assembly/mm10 RNAseq.1.trim_1_75.minQS_20.fastq RNAseq.2.trim_1_75.minQS_20.fastq  > runningLog.txt;
#	Step 3:	Remove duplicated reads (Picard)
cd ./TopHat-Out;
nohup java -Xmx2g -jar /ddn/gs1/home/lir8/Software/Picard/picard-tools-1.96/SortSam.jar INPUT= accepted_hits.bam OUTPUT= Mapped_Sorted.bam SORT_ORDER=coordinate;
nohup java -Xmx2g -jar /ddn/gs1/home/lir8/Software/Picard/picard-tools-1.96/MarkDuplicates.jar INPUT= Mapped_Sorted.bam OUTPUT= Mapped_Sorted_DeDUPLICATE.bam REMOVE_DUPLICATES=TRUE METRICS_FILE=DeDUPLICATE.txt;
#	Generate a QC summary table
#	Step 4:	Display the mapped read coverage on UCSC Genome Browser
nohup genomeCoverageBed -ibam Mapped_Sorted_DeDUPLICATE.bam -g /ddn/gs1/home/lir8/Software/RNAseq/mm10.chr.size -bg -split > RNAseq_alignment.bedgraph;
bedGraphToBigWig RNAseq_alignment.bedgraph /ddn/gs1/home/lir8/Software/RNAseq/mm10.chr.size RNAseq_alignment.bw;

#Generate DEG
cd /ddn/gs1/home/lir8/RNAseq/AC/DEG;
nohup /ddn/gs1/home/lir8/Software/RNAseq/cufflinks-2.2.1.Linux_x86_64/cuffdiff -L Wnk1ff,Wnk1dd, -p 20 --max-bundle-frags 2000000 --dispersion-method per-condition --no-update-check /ddn/gs1/home/lir8/Software/RNAseq/knownGene-mm10+Gene.gtf /ddn/gs1/home/lir8/RNAseq/PRTRIM28/F1/TopHat-Out/accepted_hits.bam,/ddn/gs1/home/lir8/RNAseq/PRTRIM28/F2/TopHat-Out/accepted_hits.bam,/ddn/gs1/home/lir8/RNAseq/PRTRIM28/F3/TopHat-Out/accepted_hits.bam,/ddn/gs1/home/lir8/RNAseq/PRTRIM28/F4/TopHat-Out/accepted_hits.bam /ddn/gs1/home/lir8/RNAseq/PRTRIM28/D1/TopHat-Out/accepted_hits.bam,/ddn/gs1/home/lir8/RNAseq/PRTRIM28/D2/TopHat-Out/accepted_hits.bam,/ddn/gs1/home/lir8/RNAseq/PRTRIM28/D3/TopHat-Out/accepted_hits.bam,/ddn/gs1/home/lir8/RNAseq/PRTRIM28/D4/TopHat-Out/accepted_hits.bam;
perl extractFPKM+4x2.pl;	