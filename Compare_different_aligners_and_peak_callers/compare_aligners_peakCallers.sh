#!/bin/bash

# # diff ./bowtie_reference_files/bowtie_hs_ref_GRCh38.p2.fa ./new_chr/new_hs_ref_GRCh38.p2.fa ## no output means these two files are exactly the same. 
# # use new_hs_ref_GRCh38.p2.fa as the reference genome across the analysis.
# hg38_genome_fa=/allen/programs/celltypes/workgroups/rnaseqanalysis/RNAseq/indexes/GRCh38/new_chr/new_hs_ref_GRCh38.p2.fa
# mm10_genome_fa=/allen/programs/celltypes/workgroups/rnaseqanalysis/RNAseq/indexes/mm10_bowtie2/mm10_bowtie2.fa

# hg38_indices_folder=/allen/programs/celltypes/workgroups/0285/Personal_folders/YD/hKeratinocyte_mKidney/STAR_hg38_indices
# mm10_indices_folder=/allen/programs/celltypes/workgroups/0285/Personal_folders/YD/hKeratinocyte_mKidney/STAR_mm10_indices

# fastq_folder=/allen/programs/celltypes/workgroups/0285/Personal_folders/YD/hKeratinocyte_mKidney
# # Run STAR to generate genome indices specifying correct path to the genome FASTA files (GTF files are no required but recommended). It takes long time and large RAM to generate indices.
# STAR --runThreadN 30 --runMode genomeGenerate --genomeDir "$hg38_indices_folder" --genomeFastaFiles "$hg38_genome_fa"
# #STAR --runThreadN 30 --runMode genomeGenerate --genomeDir "$mm10_indices_folder" --genomeFastaFiles "$mm10_genome_fa"
# # make a "run directory" and switch to it
# STAR --runThreadN 30 --genomeDir "$hg38_indices_folder" --readFilesIn "$fastq_folder"/humanKeratinocyte_R1.fastq "$fastq_folder"/humanKeratinocyte_R2.fastq --outFileNamePrefix humanKeratinocyte_STAR_

# STAR --runThreadN 30 --genomeDir "$mm10_indices_folder" --readFilesIn "$fastq_folder"/mouseKidney_R1.fastq "$fastq_folder"/mouseKidney_R2.fastq --outFileNamePrefix mouseKidney_STAR_



# #BWA Aligner (works on ibs-lurctd-ux2)
# bwa index -p hg38_genome_BWA_index -a bwtsw "$hg38_genome_fa"
# for i in $( ls | grep ^h | grep R1.fastq$ ); 
# #first make .sai files
# do bwa aln -t 8 hg38_genome_BWA_index "$i" > "$i".sai;
# bwa aln -t 8 hg38_genome_BWA_index "${i/1.fastq/}"2.fastq  > "${i/1.fastq/}"2.fastq.sai;
# #then align
# bwa sampe -a 2000 hg38_genome_BWA_index "$i".sai "${i/1.fastq/}"2.fastq.sai "$i" "${i/1.fastq/}"2.fastq > "${i/.fastq/}"_bwaSampe_aln.sam; done

# # BWA aligner for mm10 genome.
# bwa index -p mm10_genome_BWA_index -a bwtsw "$mm10_genome_fa"
# for i in $( ls | grep ^m | grep R1.fastq$ ); 
# #first make .sai files
# do bwa aln -t 8 mm10_genome_BWA_index "$i" > "$i".sai;
# bwa aln -t 8 mm10_genome_BWA_index "${i/1.fastq/}"2.fastq  > "${i/1.fastq/}"2.fastq.sai;
# #then align
# bwa sampe -a 2000 mm10_genome_BWA_index "$i".sai "${i/1.fastq/}"2.fastq.sai "$i" "${i/1.fastq/}"2.fastq > "${i/.fastq/}"_bwaSampe_aln.sam; done



# # convert SAM to BAM and drop unmapped reads (-F 4) and secondary reads (-F 256) and low quality alignments (-q 10)
# # -F 4 -F 256. This is documented in http://www.samformat.info/sam-format-flag (SAM Format Flag).
# for i in $( ls | grep _aln.sam$ ); do samtools view -bS -q 10 -F 4 -F 256 "$i" > "${i/_aln.sam/}"_aln.bam; done


# notice that the outputs from STAR are named differently (humanKeratinocyte_STAR_Aligned.out.sam). So change samtools commands accordingly.
for i in $(ls | grep .out.sam$); do samtools view -bS -q 10 -F 4 -F 256 "$i" > "${i/.out.sam/}"_aln.bam; done

#sort BAMs by chromosomal location
for i in $( ls | grep _aln.bam$ ); do samtools sort  "$i"  -o "${i/_aln.bam/}"_aln_srt.bam; done
#remove duplicates
for i in $( ls | grep _aln_srt.bam$ ); do samtools rmdup "$i" "${i/_aln_srt.bam/}"_aln_srt_rmdup.bam; done
rm *_aln.sam
rm *_aln.bam
rm *_aln_srt.bam




# The following codes are run on ibs-lurctd-ux2 because macs2 and epic2 commands are installed.
#MACS2 Peak Caller: use ibs-lurctd-ux2.
for i in $(ls | grep ^h| grep .bam$); do macs2 callpeak -t "$i" -f BAM -g hs -n "${i/.bam}"_macs -q 0.01; done
for i in $(ls | grep ^m| grep .bam$); do macs2 callpeak -t "$i" -f BAM -g mm -n "${i/.bam}"_macs -q 0.01; done


#HOMER Peak Caller 
#add path
PATH=$PATH:/allen/programs/celltypes/workgroups/mct-t200/T502/Lab_Notebook/tools/HOMER/bin
#call peaks with Homer -region on all .bam files in a directory
folder=peaks_homer_region
mkdir sam 
mkdir peaks_homer_region 
for f in $( ls | grep .bam$ ); do
             echo 'now processing:' $f
                echo 'creating sam files'
                samtools view -h $f > ./sam/$f.sam
                echo 'creating tag directories'
                makeTagDirectory ./$folder/$f ./sam/$f.sam  -format sam         
                echo 'calling peaks'
                findPeaks ./$folder/${f%.} -o auto -region
done




# SICER Peak Calling: use epic2 tool (a SICER wrapper for easier usage)
#convert to bed format. epics works on sam, single-end bam, and bed format. Here I prefer to use bed format.
for i in $( ls | grep .bam$ );
	do bamToBed -i "$i" > "${i/.bam/}".bed;		# this is using bedtools. bamToBed [options] -i <BAM>
done

#change chromosome names of human samples. This step is very important for correct outputs from SICER peak calling programs.
for i in $( ls | grep ^h | grep -v macs | grep .bed$); do python changeToBed_bamToBed.py "$i" "${i/.bed/}".rename.bed "1"; done

# Call SICER using epic2 tool.
for i in $( ls | grep ^h | grep -v macs | grep .rename.bed$); do epic2 --treatment "$i" --genome hg38 > "${i/.rename.bed/}"_SICER_peaks.txt; done

for i in $( ls | grep ^m | grep -v macs | grep .bed$); do epic2 --treatment "$i" --genome mm10 > "${i/.bed/}"_SICER_peaks.txt; done



