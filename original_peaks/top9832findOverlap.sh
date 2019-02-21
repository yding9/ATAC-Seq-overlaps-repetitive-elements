#!/bin/bash
# This script is used to revise findOverlapRepE.sh and extract top 9832 peaks per condition --> do the intersections and calculate % again. 
# try to rename files with meaningful names for easy manipulation, and put all files in the same directory for easy coding.
cd $YIDINGPATH
cd Peaks_from_JM/
mkdir top9832PeaksFolder/
find . -name *mm10.bed5 2>/dev/null > temp
while read i; do cp $i top9832PeaksFolder/; done < temp
find . -name *hg38.bed5 2>/dev/null > temp
while read i; do cp $i top9832PeaksFolder/; done < temp
find . -name *hg38_hgCon.bed5 2>/dev/null > temp
while read i; do cp $i top9832PeaksFolder/; done < temp
find . -name *hg38_hgDiv.bed5 2>/dev/null > temp
while read i; do cp $i top9832PeaksFolder/; done < temp

find . -name *mm10_mmCon.bed5 2>/dev/null > temp
while read i; do cp $i top9832PeaksFolder/; done < temp
find . -name *mm10_mmDiv.bed5 2>/dev/null > temp
while read i; do cp $i top9832PeaksFolder/; done < temp



cd top9832PeaksFolder/
# shuffle chromosome positions and then sort based on 4th column (Hommer peak calling scores) descending orders, and then output the top 9832 peaks (lines) to a new file.
for i in $(ls); do shuf $i | sort -T . -k4,4rn | head -9832 > "${i/.bed5/_top9832.bed5}"; done

for i in $( ls | grep ^h | grep top9832.bed5$ ); 
do intersectBed -wa -u -a "$i" -b ../reference/hg38.repeatMasker.bed > "${i/.bed5/_RepetitiveElements.bed5}"; 
intersectBed -wa -u -a "$i" -b ../reference/hg38.LINE.bed > "${i/.bed5/_LINE.bed5}";
intersectBed -wa -u -a "$i" -b ../reference/hg38.LTR.bed > "${i/.bed5/_LTR.bed5}";
intersectBed -wa -u -a "$i" -b ../reference/hg38.Satellite.bed > "${i/.bed5/_Satellite.bed5}";
intersectBed -wa -u -a "$i" -b ../reference/hg38.Simple_repeat.bed > "${i/.bed5/_SimpleRepeat.bed5}";
intersectBed -wa -u -a "$i" -b ../reference/hg38.SINE.bed > "${i/.bed5/_SINE.bed5}";
intersectBed -wa -u -a "$i" -b ../reference/hg38.DNA_Transposon.bed > "${i/.bed5/_DNATransposon.bed5}";
done

for i in $( ls | grep ^m | grep top9832.bed5$ );
do intersectBed -wa -u -a "$i" -b ../reference/mm10.repeatMasker.bed > "${i/.bed5/_RepetitiveElements.bed5}"; 
intersectBed -wa -u -a "$i" -b ../reference/mm10.LINE.bed > "${i/.bed5/_LINE.bed5}";
intersectBed -wa -u -a "$i" -b ../reference/mm10.LTR.bed > "${i/.bed5/_LTR.bed5}";
intersectBed -wa -u -a "$i" -b ../reference/mm10.Satellite.bed > "${i/.bed5/_Satellite.bed5}";
intersectBed -wa -u -a "$i" -b ../reference/mm10.Simple_repeat.bed > "${i/.bed5/_SimpleRepeat.bed5}";
intersectBed -wa -u -a "$i" -b ../reference/mm10.SINE.bed > "${i/.bed5/_SINE.bed5}";
intersectBed -wa -u -a "$i" -b ../reference/mm10.DNA_Transposon.bed > "${i/.bed5/_DNATransposon.bed5}"; done

##########################################################################################################################################
# Count the number of overlap peaks. Note this is an important step to generate the right count matrix in R!
for i in $( ls | grep top9832 | grep RepetitiveElements | sort); do wc -l "$i" | tr -d '\n' >> ../overlap_top9832_summary.txt; printf " " >> ../overlap_top9832_summary.txt; wc -l "${i/_RepetitiveElements.bed5/.bed5}" >> ../overlap_top9832_summary.txt; done
for i in $( ls | grep top9832 | grep LINE | sort); do wc -l "$i" | tr -d '\n' >> ../overlap_top9832_summary.txt; printf " " >> ../overlap_top9832_summary.txt; wc -l "${i/_LINE.bed5/.bed5}" >> ../overlap_top9832_summary.txt; done
for i in $( ls | grep top9832 | grep LTR | sort); do wc -l "$i" | tr -d '\n' >> ../overlap_top9832_summary.txt; printf " " >> ../overlap_top9832_summary.txt; wc -l "${i/_LTR.bed5/.bed5}" >> ../overlap_top9832_summary.txt; done
for i in $( ls | grep top9832 | grep Satellite | sort); do wc -l "$i" | tr -d '\n' >> ../overlap_top9832_summary.txt; printf " " >> ../overlap_top9832_summary.txt; wc -l "${i/_Satellite.bed5/.bed5}" >> ../overlap_top9832_summary.txt; done
for i in $( ls | grep top9832 | grep SimpleRepeat | sort); do wc -l "$i" | tr -d '\n' >> ../overlap_top9832_summary.txt; printf " " >> ../overlap_top9832_summary.txt; wc -l "${i/_SimpleRepeat.bed5/.bed5}" >> ../overlap_top9832_summary.txt; done
for i in $( ls | grep top9832 | grep SINE | sort); do wc -l "$i" | tr -d '\n' >> ../overlap_top9832_summary.txt; printf " " >> ../overlap_top9832_summary.txt; wc -l "${i/_SINE.bed5/.bed5}" >> ../overlap_top9832_summary.txt; done
for i in $( ls | grep top9832 | grep DNATransposon | sort); do wc -l "$i" | tr -d '\n' >> ../overlap_top9832_summary.txt; printf " " >> ../overlap_top9832_summary.txt; wc -l "${i/_DNATransposon.bed5/.bed5}" >> ../overlap_top9832_summary.txt; done
