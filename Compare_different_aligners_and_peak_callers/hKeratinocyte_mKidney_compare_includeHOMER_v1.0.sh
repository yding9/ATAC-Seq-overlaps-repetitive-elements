#!/bin/bash
## align fastq files using John Mich's method. Revised from allInstructions.sh
#gunzip *.gz;

TRIM_GALORE=//allen/programs/celltypes/workgroups/mct-t200/T502/Lab_Notebook/tools/trim_galore
BOWTIE=//allen/programs/celltypes/workgroups/mct-t200/T502/Lab_Notebook/tools/bowtie-1.1.0/bowtie
BOWTIE1_INDEX_MM10=//allen/programs/celltypes/workgroups/rnaseqanalysis/RNAseq/indexes/mm10/genome
BOWTIE1_INDEX_HG38=//allen/programs/celltypes/workgroups/rnaseqanalysis/RNAseq/indexes/GRCh38/bowtie_reference_files/bowtie_hs_ref_GRCh38.p2

for i in $( ls | grep ^h | grep _R1.fastq$ );
do bowtie2 -p 8 --no-mixed --no-discordant -X 2000 -x /allen/programs/celltypes/workgroups/rnaseqanalysis/RNAseq/indexes/GRCh38/bowtie2_reference_files/bowtie2_hs_ref_GRCh38.p2 -1 "$i" -2 "${i/1.fastq/}"2.fastq -S "${i/_R1.fastq/}"_aln.sam;
done

for i in $( ls | grep ^m | grep _R1.fastq$ );
do bowtie2 -p 8 --no-mixed --no-discordant -X 2000 -x /allen/programs/celltypes/workgroups/rnaseqanalysis/RNAseq/indexes/mm10_bowtie2/mm10_bowtie2 -1 "$i" -2 "${i/1.fastq/}"2.fastq -S "${i/_R1.fastq/}"_aln.sam;
done
# convert SAM to BAM and drop unmapped reads (-F 4) and secondary reads (-F 256) and low quality alignments (-q 10)
# -F 4 -F 256. This is documented in http://www.samformat.info/sam-format-flag (SAM Format Flag).
for i in $( ls | grep _aln.sam$ ); do samtools view -bS -q 10 -F 4 -F 256 "$i" > "${i/_aln.sam/}"_aln.bam; done
#sort BAMs by chromosomal location
for i in $( ls | grep _aln.bam$ ); do samtools sort  "$i"  -o "${i/_aln.bam/}"_aln_srt.bam; done
#remove duplicates
for i in $( ls | grep _aln_srt.bam$ ); do samtools rmdup "$i" "${i/_aln_srt.bam/}"_aln_srt_rmdup.bam; done
rm *_aln.sam
rm *_aln.bam
rm *_aln_srt.bam


##################################################################################################################################################################################################
# align fastq files using Lucas's method. Revised from Lucas' 00_config.sh and 01_alignment.sh.

$BOWTIE -p 8 -m 1 -S -X 2000 --chunkmbs 512 --un mKidney_bowtie1_unaligned.fastq $BOWTIE1_INDEX_MM10 -1 mouseKidney_R1.fastq -2 mouseKidney_R2.fastq mKidney_bowtie1_aligned.sam
$BOWTIE -p 8 -m 1 -S -X 2000 --chunkmbs 512 --un hKeratinocyte_bowtie1_unaligned.fastq $BOWTIE1_INDEX_HG38 -1 humanKeratinocyte_R1.fastq -2 humanKeratinocyte_R2.fastq hKeratinocyte_bowtie1_aligned.sam
# convert SAM results to BAM
samtools view -bS mKidney_bowtie1_aligned.sam > mKidney_bowtie1_aligned.bam
samtools view -bS hKeratinocyte_bowtie1_aligned.sam > hKeratinocyte_bowtie1_aligned.bam
#trim Nextera primer sequences from unaligned reads
$TRIM_GALORE --nextera --paired --three_prime_clip_R1 1 --three_prime_clip_R2 1 mKidney_bowtie1_unaligned_1.fastq mKidney_bowtie1_unaligned_2.fastq
$TRIM_GALORE --nextera --paired --three_prime_clip_R1 1 --three_prime_clip_R2 1 hKeratinocyte_bowtie1_unaligned_1.fastq hKeratinocyte_bowtie1_unaligned_2.fastq
# align trimmed reads.
$BOWTIE -p 8 -m 1 -S -X 2000 --chunkmbs 512 $BOWTIE1_INDEX_MM10 -1 mKidney_bowtie1_unaligned_1_val_1.fq -2 mKidney_bowtie1_unaligned_2_val_2.fq mKidney_bowtie1_ump_trim_aligned.sam
$BOWTIE -p 8 -m 1 -S -X 2000 --chunkmbs 512 $BOWTIE1_INDEX_HG38 -1 hKeratinocyte_bowtie1_unaligned_1_val_1.fq -2 hKeratinocyte_bowtie1_unaligned_2_val_2.fq hKeratinocyte_bowtie1_ump_trim_aligned.sam
# convert trimeed read results to BAM.
samtools view -bS mKidney_bowtie1_ump_trim_aligned.sam > mKidney_bowtie1_ump_trim_aligned.bam
samtools view -bS hKeratinocyte_bowtie1_ump_trim_aligned.sam > hKeratinocyte_bowtie1_ump_trim_aligned.bam
# concatenate initial and trimmed results
samtools cat -o mKidney_bowtie1_merged.bam mKidney_bowtie1_aligned.bam mKidney_bowtie1_ump_trim_aligned.bam
samtools cat -o hKeratinocyte_bowtie1_merged.bam hKeratinocyte_bowtie1_aligned.bam hKeratinocyte_bowtie1_ump_trim_aligned.bam
# retain only aligned pairs based on BAM flags
samtools view -b -F 4 mKidney_bowtie1_merged.bam > mKidney_bowtie1_mapped.bam
samtools view -b -F 4 hKeratinocyte_bowtie1_merged.bam > hKeratinocyte_bowtie1_mapped.bam
# sort mapped reads
samtools sort mKidney_bowtie1_mapped.bam -o mKidney_bowtie1_mapped.srt.bam
samtools sort hKeratinocyte_bowtie1_mapped.bam -o hKeratinocyte_bowtie1_mapped.srt.bam
# remove duplicates
samtools rmdup mKidney_bowtie1_mapped.srt.bam mKidney_bowtie1_mapped.rmd.srt.bam
samtools rmdup hKeratinocyte_bowtie1_mapped.srt.bam hKeratinocyte_bowtie1_mapped.rmd.srt.bam

rm *_aligned.sam
rm *_aligned.bam
rm *_mapped.bam
rm *_merged.bam
rm *_mapped.srt.bam
rm $(ls | grep unaligned)

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



### -c: Instead of printing the alignments, only count them and print the total number. All filter options, such as -f, -F, and -q, are taken into account.
for i in $(ls | grep .bam$); do samtools view -c $i; echo $i; done

### calculate the overlap with repeat elements.
reference=/allen/programs/celltypes/workgroups/0285/Personal_folders/YD/Peaks_from_JM/reference

for i in $( ls | grep ^h | grep .bed5$ );
do intersectBed -wa -u -a "$i" -b $reference/hg38.repeatMasker.bed > "${i/.bed5/_olpRM.bed5}";
intersectBed -wa -u -a "$i" -b $reference/hg38.LINE.bed > "${i/.bed5/_olpLINE.bed5}";
intersectBed -wa -u -a "$i" -b $reference/hg38.LTR.bed > "${i/.bed5/_olpLTR.bed5}";
intersectBed -wa -u -a "$i" -b $reference/hg38.Satellite.bed > "${i/.bed5/_olpSatellite.bed5}";
intersectBed -wa -u -a "$i" -b $reference/hg38.Simple_repeat.bed > "${i/.bed5/_olpSimple_repeat.bed5}";
intersectBed -wa -u -a "$i" -b $reference/hg38.SINE.bed > "${i/.bed5/_olpSINE.bed5}";
intersectBed -wa -u -a "$i" -b $reference/hg38.DNA_Transposon.bed > "${i/.bed5/_olpDNA_Transposon.bed5}";
done

for i in $( ls | grep ^m | grep .bed5$ );
do intersectBed -wa -u -a "$i" -b $reference/mm10.repeatMasker.bed > "${i/.bed5/_olpRM.bed5}";
intersectBed -wa -u -a "$i" -b $reference/mm10.LINE.bed > "${i/.bed5/_olpLINE.bed5}";
intersectBed -wa -u -a "$i" -b $reference/mm10.LTR.bed > "${i/.bed5/_olpLTR.bed5}";
intersectBed -wa -u -a "$i" -b $reference/mm10.Satellite.bed > "${i/.bed5/_olpSatellite.bed5}";
intersectBed -wa -u -a "$i" -b $reference/mm10.Simple_repeat.bed > "${i/.bed5/_olpSimple_repeat.bed5}";
intersectBed -wa -u -a "$i" -b $reference/mm10.SINE.bed > "${i/.bed5/_olpSINE.bed5}";
intersectBed -wa -u -a "$i" -b $reference/mm10.DNA_Transposon.bed > "${i/.bed5/_olpDNA_Transposon.bed5}";
done

# Count the percentage of overlap peaks to repetitive elements.
for i in $(ls | grep _olpRM.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> overlap_summary.txt; printf " " >> overlap_summary.txt; wc -l "${i/_olpRM.bed5/.bed5}" >> overlap_summary.txt; done

for i in $(ls | grep _olpLINE.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> overlap_summary.txt; printf " " >> overlap_summary.txt; wc -l "${i/_olpLINE.bed5/.bed5}" >> overlap_summary.txt; done

for i in $(ls | grep _olpLTR.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> overlap_summary.txt; printf " " >> overlap_summary.txt; wc -l "${i/_olpLTR.bed5/.bed5}" >> overlap_summary.txt; done

for i in $(ls | grep _olpSatellite.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> overlap_summary.txt; printf " " >> overlap_summary.txt; wc -l "${i/_olpSatellite.bed5/.bed5}" >> overlap_summary.txt; done

for i in $(ls | grep _olpSimple_repeat.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> overlap_summary.txt; printf " " >> overlap_summary.txt; wc -l "${i/_olpSimple_repeat.bed5/.bed5}" >> overlap_summary.txt; done

for i in $(ls | grep _olpSINE.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> overlap_summary.txt; printf " " >> overlap_summary.txt; wc -l "${i/_olpSINE.bed5/.bed5}" >> overlap_summary.txt; done

for i in $(ls | grep _olpDNA_Transposon.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> overlap_summary.txt; printf " " >> overlap_summary.txt; wc -l "${i/_olpDNA_Transposon.bed5/.bed5}" >> overlap_summary.txt; done

