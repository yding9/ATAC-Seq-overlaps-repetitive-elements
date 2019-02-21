#!/bin/bash
##########################################################################################################################################################
head -50 hg38.fa.out > test1.out
# convert repeat masked out file (with header information) to bed file format      (https://www.biostars.org/p/128068/)
#awk 'BEGIN{OFS="\t"}{if(NR>3) {if($9=="C"){strand="-"}else{strand="+"};print $5,$6-1,$7,$10,".",strand}}' test1.out > test1.bed
#awk 'BEGIN{OFS="\t"}{if(NR>3) {if($9=="C"){strand="-"}else{strand="+"};print $5,$6-1,$7,$10,".",strand}}' hg38.fa.out > hg38.repeatMasker.bed
awk 'BEGIN{OFS="\t"}{if(NR>3) {if($9=="C"){strand="-"}else{strand="+"};print $5,$6-1,$7,$10,$11,strand}}' hg38.fa.out > hg38.repeatMasker.bed
# To convert repeated masked out file (without header information) to bed file format
#awk '{print $5,$6-1, $7,$10, $11}' hg38.repeatMasker.bed > test2_LINE.bed


# To see a frequency count for column 11 (repeat class/family)
#awk '{print $11}' hg38.fa.out | sort | uniq -c > hg38.classtype
awk '{print $5}' hg38.repeatMasker.bed | sort | uniq -c > hg38.classtype.txt
#To print the whole fields ($0) when the 11st field contains "LINE..".
awk '$5 ~ "LINE" {print $0}' hg38.repeatMasker.bed > hg38.LINE.bed
awk '$5 ~ "SINE" {print $0}' hg38.repeatMasker.bed > hg38.SINE.bed
awk '$5 ~ "LTR" {print $0}' hg38.repeatMasker.bed > hg38.LTR.bed
awk '$5 ~ "Simple_repeat" {print $0}' hg38.repeatMasker.bed > hg38.Simple_repeat.bed
awk '$5 ~ "Satellite" {print $0}' hg38.repeatMasker.bed > hg38.Satellite.bed
awk '$5 ~ "DNA" {print $0}' hg38.repeatMasker.bed > hg38.DNA_Transposon.bed

# convert repeat masked out file (with header information) to bed file format      (https://www.biostars.org/p/128068/)
awk 'BEGIN{OFS="\t"}{if(NR>3) {if($9=="C"){strand="-"}else{strand="+"};print $5,$6-1,$7,$10,$11,strand}}' mm10.fa.out > mm10.repeatMasker.bed
# To see a frequency count for column 11 (repeat class/family)
awk '{print $5}' mm10.repeatMasker.bed | sort | uniq -c > mm10.classtype.txt
#To print the whole fields ($0) when the 11st field contains "LINE..".
awk '$5 ~ "LINE" {print $0}' mm10.repeatMasker.bed > mm10.LINE.bed
awk '$5 ~ "SINE" {print $0}' mm10.repeatMasker.bed > mm10.SINE.bed
awk '$5 ~ "LTR" {print $0}' mm10.repeatMasker.bed > mm10.LTR.bed
awk '$5 ~ "Simple_repeat" {print $0}' mm10.repeatMasker.bed > mm10.Simple_repeat.bed
awk '$5 ~ "Satellite" {print $0}' mm10.repeatMasker.bed > mm10.Satellite.bed
awk '$5 ~ "DNA" {print $0}' mm10.repeatMasker.bed > mm10.DNA_Transposon.bed
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################


####################################################################################################################################################################################################
# The code below is trying to find conserved and divergent peaks in human and mouse respectively.
####################################################################################################################################################################################################

#liftOver mm10 to hg38
for i in $( ls | grep ^m | grep .bed5$ ); do liftOver -minMatch=0.6 -bedPlus=3 "$i" mm10ToHg38.over.chain "$i".hg38 unmapped; done   #grep ^m: search files(directories) starting with "m". | grep .bed5$: search .bed5 in files. $ indicate special characters like \t, \n.

#liftOver hg38 to mm10
for i in $( ls | grep ^h | grep .bed5$ ); do liftOver -minMatch=0.6 -bedPlus=3 "$i" hg38ToMm10.over.chain "$i".mm10 unmapped; done

# intersect human 11 subclass peaks to corresponding mouse 11 subclass peaks to find human conserved peaks for these 11 subclasses.  
for i in $( ls | grep ^h | grep hg38.bed5$ ); 
do intersectBed -wa -u -a "$i" -b m"$(echo $i | cut -d '_' -f 1 | cut -c 2- )"_peaks_mm10.bed5.hg38 > "${i/.bed5/_hgCon.bed5}"; done
    
# find divergent human peaks for 11 subclasses using the -v option (Reporting the absence of any overlapping features)
for i in $( ls | grep ^h | grep hg38.bed5$ ); 
do intersectBed -wa -v -a "$i" -b m"$(echo $i | cut -d '_' -f 1 | cut -c 2- )"_peaks_mm10.bed5.hg38 > "${i/.bed5/_hgDiv.bed5}"; done
    
# intersect mouse 11 subclass peaks to corresponding human 11 subclass peaks to find mouse conserved peaks for these 11 subclasses.
for i in $( ls | grep ^m | grep mm10.bed5$ ); 
do intersectBed -wa -u -a "$i" -b h"$(echo $i | cut -d '_' -f 1 | cut -c 2- )"_peaks_hg38.bed5.mm10 > "${i/.bed5/_mmCon.bed5}"; done

# find divergent mouse peaks using the -v option.
for i in $( ls | grep ^m | grep mm10.bed5$ ); 
do intersectBed -wa -v -a "$i" -b h"$(echo $i | cut -d '_' -f 1 | cut -c 2- )"_peaks_hg38.bed5.mm10 > "${i/.bed5/_mmDiv.bed5}"; done

##########################################################################################################################################################
##########################################################################################################################################################
  
# ONLY REPORT original A entry once if any overlaps found in B. In other words, just report the fact at least one overlap was found in B.
#sort and output only UNIQUE entities in olpRM.bed5 files. This is the byproduct by bedtools intersectBed (if there are multiple hits for A, it will output A multiple times.)
# actually this problem can be solved by using -u option with intersectBed function. So use that option instead. 
#for i in $( ls | grep _olpRM.bed5$ );
	#do sort -k1,1 -k2,2n -u "$i" > "${i/_olpRM.bed5/_olpRM_UNIQ.bed5}";
#done

# calculate intersections between peak bed files and hg38 subliraries from hg38 repeatMasker
for i in $( ls | grep ^h | grep .bed5$ );
do intersectBed -wa -u -a "$i" -b ./reference/hg38.repeatMasker.bed > "${i/.bed5/_olpRM.bed5}";
intersectBed -wa -u -a "$i" -b ./reference/hg38.LINE.bed > "${i/.bed5/_olpLINE.bed5}";
intersectBed -wa -u -a "$i" -b ./reference/hg38.LTR.bed > "${i/.bed5/_olpLTR.bed5}";
intersectBed -wa -u -a "$i" -b ./reference/hg38.Satellite.bed > "${i/.bed5/_olpSatellite.bed5}";
intersectBed -wa -u -a "$i" -b ./reference/hg38.Simple_repeat.bed > "${i/.bed5/_olpSimple_repeat.bed5}";
intersectBed -wa -u -a "$i" -b ./reference/hg38.SINE.bed > "${i/.bed5/_olpSINE.bed5}";
intersectBed -wa -u -a "$i" -b ./reference/hg38.DNA_Transposon.bed > "${i/.bed5/_olpDNA_Transposon.bed5}";
done

# calculate intersections between peak bed files and mm10 subliraries from mm10 repeatMasker
for i in $( ls | grep ^m | grep .bed5$ );
do intersectBed -wa -u -a "$i" -b ./reference/mm10.repeatMasker.bed > "${i/.bed5/_olpRM.bed5}";
intersectBed -wa -u -a "$i" -b ./reference/mm10.LINE.bed > "${i/.bed5/_olpLINE.bed5}";
intersectBed -wa -u -a "$i" -b ./reference/mm10.LTR.bed > "${i/.bed5/_olpLTR.bed5}";
intersectBed -wa -u -a "$i" -b ./reference/mm10.Satellite.bed > "${i/.bed5/_olpSatellite.bed5}";
intersectBed -wa -u -a "$i" -b ./reference/mm10.Simple_repeat.bed > "${i/.bed5/_olpSimple_repeat.bed5}";
intersectBed -wa -u -a "$i" -b ./reference/mm10.SINE.bed > "${i/.bed5/_olpSINE.bed5}";
intersectBed -wa -u -a "$i" -b ./reference/mm10.DNA_Transposon.bed > "${i/.bed5/_olpDNA_Transposon.bed5}";
done



# Note: uniq isn’t able to detect the duplicate lines unless they are adjacent. The content in the file must be therefore sorted before using uniq or you can simply use sort -u instead f uniq.
#presort your data by chromosome and then by start position (e.g., sort -k1,1 -k2,2n in.bed > in.sorted.bed for BED files)


##########################################################################################################################################
# Count the number of overlap peaks. Note this is an important step to generate the right count matrix in R!
# Use the following count_lines.sh

#!/bin/bash

for i in $(ls | grep _olpRM.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> ../overlap_summary.txt; printf " " >> ../overlap_summary.txt; wc -l "${i/_olpRM.bed5/.bed5}" >> ../overlap_summary.txt; done

for i in $(ls | grep _olpLINE.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> ../overlap_summary.txt; printf " " >> ../overlap_summary.txt; wc -l "${i/_olpLINE.bed5/.bed5}" >> ../overlap_summary.txt; done

for i in $(ls | grep _olpLTR.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> ../overlap_summary.txt; printf " " >> ../overlap_summary.txt; wc -l "${i/_olpLTR.bed5/.bed5}" >> ../overlap_summary.txt; done

for i in $(ls | grep _olpSatellite.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> ../overlap_summary.txt; printf " " >> ../overlap_summary.txt; wc -l "${i/_olpSatellite.bed5/.bed5}" >> ../overlap_summary.txt; done

for i in $(ls | grep _olpSimple_repeat.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> ../overlap_summary.txt; printf " " >> ../overlap_summary.txt; wc -l "${i/_olpSimple_repeat.bed5/.bed5}" >> ../overlap_summary.txt; done

for i in $(ls | grep _olpSINE.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> ../overlap_summary.txt; printf " " >> ../overlap_summary.txt; wc -l "${i/_olpSINE.bed5/.bed5}" >> ../overlap_summary.txt; done

for i in $(ls | grep _olpDNA_Transposon.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> ../overlap_summary.txt; printf " " >> ../overlap_summary.txt; wc -l "${i/_olpDNA_Transposon.bed5/.bed5}" >> ../overlap_summary.txt; done


