#####install cell ranger to process the fastq files
Download cellranger,
tar xzvf cellranger-7.1.0.tar.xz

export PATH="/storage/hpc/group/ircf/users/rlgch/Software/scRNA/cellranger-7.0.1:$PATH"
source ~/.bash_profile

cellranger testrun --id=tiny

/storage/hpc/group/ircf/users/rlgch/Software/scRNA/cellranger-7.0.1/cellranger count --id=cellranger-DD \
                   --transcriptome=/storage/hpc/group/ircf/users/rlgch/Software/reference/10X/refdata-gex-mm10-2020-A \
                   --fastqs=/storage/hpc/group/ircf/users/rlgch/scRNA/TRIM28/DD/ \
                   --sample=DD \

/storage/hpc/group/ircf/users/rlgch/Software/scRNA/cellranger-7.0.1/cellranger count --id=cellranger-Flox \
                   --transcriptome=/storage/hpc/group/ircf/users/rlgch/Software/reference/10X/refdata-gex-mm10-2020-A \
                   --fastqs=/storage/hpc/group/ircf/users/rlgch/scRNA/TRIM28/Flox/ \
                   --sample=Flox \


####install Velocyto in conda environment "RNAvelocity" and generate loom file
srun --partition Interactive --qos interactive --pty /bin/bash
module load miniconda3
conda env create --file TutorialEnvironment.yml
conda activate RNAvelocity
module load samtools/1.14
velocyto run10x /storage/hpc/group/ircf/users/rlgch/scRNA/TRIM28/cellranger-DD /storage/hpc/group/ircf/users/rlgch/Software/reference/10X/refdata-gex-mm10-2020-A/genes/genes.gtf
velocyto run10x /storage/hpc/group/ircf/users/rlgch/scRNA/TRIM28/cellranger-Flox /storage/hpc/group/ircf/users/rlgch/Software/reference/10X/refdata-gex-mm10-2020-A/genes/genes.gtf
