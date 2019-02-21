# This is an advanced version of bootstrap_overlap_repetitiveElements.sh. Directly output the percentage of overlapping in txt files.
# ls | grep -v top9832
# cp $(ls | grep -v top9832) ../bootstrap
#First do 100X bootstrapping and calculate intersections between subsets of peaks and reference repetitive elements files. (Real data)
for i in $(ls | grep ^h | grep .bed5$);
do number_lines=$((` wc -l "$i" | cut -f1 -d " "`))
number_to_subsample=$(($number_lines * 8 / 10))
for j in $( seq 1 100);
do shuf -n $number_to_subsample "$i" > temp;
intersectBed -wa -u -a temp -b ../reference/hg38.repeatMasker.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_repeatMaskerintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/hg38.LINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_LINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/hg38.LTR.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_LTRintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/hg38.Satellite.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_Satelliteintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/hg38.Simple_repeat.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_SimpleRepeatintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/hg38.SINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_SINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/hg38.DNA_Transposon.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_DNATransposonintsxns100xCounts.txt;
done
done

# Second do 100X shuffling and calculate intersections between randomized peaksets and reference repetitive elements files (Randomized data).
for i in $(ls | grep ^h | grep .bed5$);
do for j in $( seq 1 100);
do shuffleBed -noOverlapping -i "$i" -g hg38ChromSizes.txt > temp;
number_lines_shuffle=$((` wc -l temp | cut -f1 -d " "`))
intersectBed -wa -u -a temp -b ../reference/hg38.repeatMasker.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_repeatMaskerintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/hg38.LINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_LINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/hg38.LTR.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_LTRintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/hg38.Satellite.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_Satelliteintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/hg38.Simple_repeat.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_SimpleRepeatintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/hg38.SINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_SINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/hg38.DNA_Transposon.bed > temp2;
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
intersectBed -wa -u -a temp -b ../reference/mm10.repeatMasker.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_repeatMaskerintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/mm10.LINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_LINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/mm10.LTR.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_LTRintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/mm10.Satellite.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_Satelliteintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/mm10.Simple_repeat.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_SimpleRepeatintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/mm10.SINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_SINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/mm10.DNA_Transposon.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_to_subsample. | python >> "${i/.bed5/}"_DNATransposonintsxns100xCounts.txt;
done
done


# Second do 100X shuffling and calculate intersections between randomized peaksets and reference repetitive elements files (Randomized data).
for i in $(ls | grep ^m | grep .bed5$);
do for j in $( seq 1 100);
do shuffleBed -noOverlapping -i "$i" -g mm10ChromSizes.txt > temp;
number_lines_shuffle=$((` wc -l temp | cut -f1 -d " "`))
intersectBed -wa -u -a temp -b ../reference/mm10.repeatMasker.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_repeatMaskerintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/mm10.LINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_LINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/mm10.LTR.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_LTRintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/mm10.Satellite.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_Satelliteintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/mm10.Simple_repeat.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_SimpleRepeatintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/mm10.SINE.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_SINEintsxns100xCounts.txt;
intersectBed -wa -u -a temp -b ../reference/mm10.DNA_Transposon.bed > temp2;
echo print $(wc -l temp2 | cut -f1 -d " ")/$number_lines_shuffle. | python >> "${i/.bed5/}"_shuffled_DNATransposonintsxns100xCounts.txt;
done
done
rm temp;
rm temp2



