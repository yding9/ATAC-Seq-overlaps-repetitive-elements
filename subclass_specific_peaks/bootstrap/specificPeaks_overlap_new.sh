#!/bin/bash
# to make specific peak files in bash, assuming you have a directory with all your human and mouse .bed5 files
# when generating this script on Windows, it always has carriafe return at the end of lines. So use the following command line in the terminal to run the following script.
# cat specificPeak_OverlapRepE.sh | tr -d "\r" | bash
### If this script is generated on Mac or Unix system, directly use ./specificPeaks_overlap_new.sh to run in terminal.
cd ./OriginalPeaks_fromJohn
for letter in $( echo h m ); do
for i in $(ls | grep ^"$letter" | grep .bed5$ );
do ls | grep ^"$letter" | grep .bed5$ | grep -v "$i" > temp;
while read line; 
do cat "$line" >> temp2;
done < temp
intersectBed -v -wa -a "$i" -b temp2 > "${i/.bed5/_specific.bed5}";
rm temp2;
done
done

####################################################################################################################################################################################################
## The code below is trying to find conserved and divergent peaks in human and mouse respectively.
####################################################################################################################################################################################################

#liftOver mm10 to hg38
for i in $( ls | grep ^m | grep _specific.bed5$ );
do liftOver -minMatch=0.6 -bedPlus=3 "$i" mm10ToHg38.over.chain "$i".hg38 unmapped;
done   
#grep ^m: search files(directories) starting with "m". | grep .bed5$: search .bed5 in files. $ indicate special characters like \t, \n.

#liftOver hg38 to mm10
for i in $( ls | grep ^h | grep _specific.bed5$ );
do liftOver -minMatch=0.6 -bedPlus=3 "$i" hg38ToMm10.over.chain "$i".mm10 unmapped;
done

# intersect human peaks to corresponding mouse peaks to find human conserved peaks.  
for i in $( ls | grep ^h | grep hg38_specific.bed5$ ); 
do intersectBed -wa -u -a "$i" -b m"$(echo $i | cut -d '_' -f 1 | cut -c 2- )"_peaks_mm10_specific.bed5.hg38 > "${i/.bed5/_hgCon.bed5}"; done
    
# find divergent human peaks for 11 subclasses using the -v option (Reporting the absence of any overlapping features)
for i in $( ls | grep ^h | grep hg38_specific.bed5$ ); 
do intersectBed -wa -v -a "$i" -b m"$(echo $i | cut -d '_' -f 1 | cut -c 2- )"_peaks_mm10_specific.bed5.hg38 > "${i/.bed5/_hgDiv.bed5}"; done
    
# intersect mouse peaks to corresponding human peaks to find mouse conserved peaks.
for i in $( ls | grep ^m | grep mm10_specific.bed5$ ); 
do intersectBed -wa -u -a "$i" -b h"$(echo $i | cut -d '_' -f 1 | cut -c 2- )"_peaks_hg38_specific.bed5.mm10 > "${i/.bed5/_mmCon.bed5}"; done

# find divergent mouse peaks using the -v option.
for i in $( ls | grep ^m | grep mm10_specific.bed5$ ); 
do intersectBed -wa -v -a "$i" -b h"$(echo $i | cut -d '_' -f 1 | cut -c 2- )"_peaks_hg38_specific.bed5.mm10 > "${i/.bed5/_mmDiv.bed5}"; done

#########################################################################################################################################################
# use bedtools to map ATAC-seq peaks to repetitive elements in the genome.
#########################################################################################################################################################
 
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


##########################################################################################################################################
# Count the percentage of overlap peaks to repetitive elements.
for i in $(ls | grep _olpRM.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> specific_overlap_summary.txt; printf " " >> specific_overlap_summary.txt; wc -l "${i/_olpRM.bed5/.bed5}" >> specific_overlap_summary.txt; done

for i in $(ls | grep _olpLINE.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> specific_overlap_summary.txt; printf " " >> specific_overlap_summary.txt; wc -l "${i/_olpLINE.bed5/.bed5}" >> specific_overlap_summary.txt; done

for i in $(ls | grep _olpLTR.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> specific_overlap_summary.txt; printf " " >> specific_overlap_summary.txt; wc -l "${i/_olpLTR.bed5/.bed5}" >> specific_overlap_summary.txt; done

for i in $(ls | grep _olpSatellite.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> specific_overlap_summary.txt; printf " " >> specific_overlap_summary.txt; wc -l "${i/_olpSatellite.bed5/.bed5}" >> specific_overlap_summary.txt; done

for i in $(ls | grep _olpSimple_repeat.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> specific_overlap_summary.txt; printf " " >> specific_overlap_summary.txt; wc -l "${i/_olpSimple_repeat.bed5/.bed5}" >> specific_overlap_summary.txt; done

for i in $(ls | grep _olpSINE.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> specific_overlap_summary.txt; printf " " >> specific_overlap_summary.txt; wc -l "${i/_olpSINE.bed5/.bed5}" >> specific_overlap_summary.txt; done

for i in $(ls | grep _olpDNA_Transposon.bed5$ | sort); do wc -l "$i" | tr -d '\n' >> specific_overlap_summary.txt; printf " " >> specific_overlap_summary.txt; wc -l "${i/_olpDNA_Transposon.bed5/.bed5}" >> specific_overlap_summary.txt; done



##########################################################################################################################################
#First do 100X bootstrapping and calculate intersections between subsets of peaks and reference repetitive elements files. (Real data)
for i in $(ls | grep ^h | grep .bed5$);
do number_lines=$((` wc -l "$i" | cut -f1 -d " "`))
number_to_subsample=$(($number_lines * 8 / 10))
for j in $( seq 1 100);
do shuf -n $number_to_subsample "$i" > temp;
intersectBed -wa -u -a temp -b ./reference/hg38.repeatMasker.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_repeatMaskerintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/hg38.LINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_LINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/hg38.LTR.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_LTRintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/hg38.Satellite.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_Satelliteintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/hg38.Simple_repeat.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_SimpleRepeatintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/hg38.SINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_SINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/hg38.DNA_Transposon.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_DNATransposonintsxns100xCounts.txt;
done
done

# Second do 100X shuffling and calculate intersections between randomized peaksets and reference repetitive elements files (Randomized data).
for i in $(ls | grep ^h | grep .bed5$);
do for j in $( seq 1 100);
do shuffleBed -noOverlapping -i "$i" -g hg38ChromSizes.txt > temp;
number_lines_shuffle=$((` wc -l temp | cut -f1 -d " "`))
intersectBed -wa -u -a temp -b ./reference/hg38.repeatMasker.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_repeatMaskerintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/hg38.LINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_LINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/hg38.LTR.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_LTRintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/hg38.Satellite.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_Satelliteintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/hg38.Simple_repeat.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_SimpleRepeatintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/hg38.SINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_SINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/hg38.DNA_Transposon.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_DNATransposonintsxns100xCounts.txt;
done
done

################ The following is to calculate mouse peaks overlapping with mouse repetitive elements libraries.
#First do 100X bootstrapping and calculate intersections between subsets of peaks and reference repetitive elements files. (Real data)
for i in $(ls | grep ^m | grep .bed5$);
do number_lines=$((` wc -l "$i" | cut -f1 -d " "`))
number_to_subsample=$(($number_lines * 8 / 10))
for j in $( seq 1 100);
do shuf -n $number_to_subsample "$i" > temp;
intersectBed -wa -u -a temp -b ./reference/mm10.repeatMasker.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_repeatMaskerintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/mm10.LINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_LINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/mm10.LTR.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_LTRintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/mm10.Satellite.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_Satelliteintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/mm10.Simple_repeat.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_SimpleRepeatintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/mm10.SINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_SINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/mm10.DNA_Transposon.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_DNATransposonintsxns100xCounts.txt;
done
done


# Second do 100X shuffling and calculate intersections between randomized peaksets and reference repetitive elements files (Randomized data).
for i in $(ls | grep ^m | grep .bed5$);
do for j in $( seq 1 100);
do shuffleBed -noOverlapping -i "$i" -g mm10ChromSizes.txt > temp;
number_lines_shuffle=$((` wc -l temp | cut -f1 -d " "`))
intersectBed -wa -u -a temp -b ./reference/mm10.repeatMasker.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_repeatMaskerintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/mm10.LINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_LINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/mm10.LTR.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_LTRintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/mm10.Satellite.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_Satelliteintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/mm10.Simple_repeat.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_SimpleRepeatintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/mm10.SINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_SINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ./reference/mm10.DNA_Transposon.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_DNATransposonintsxns100xCounts.txt;
done
done
rm temp;
rm temp2






