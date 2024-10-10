summ()
{
#Run dorado summary
$summary $1 > ${1}.Summary.tsv
#Get read length
awk '{print $10}' ${1}.Summary.tsv > ${1}.length.distribution
#Remove column name
sed -i '1d' ${1}.length.distribution
#extract lowest and highest values
lowest=$(head -n 1 ${1}.length.distribution)
highest=$(tail -n 1 ${1}.length.distribution)

#get mean
sum=$(awk '{s+=$1} END {print s}' ${1}.length.distribution)
count=$(wc -l < ${1}.length.distribution)
mean=$(echo "scale=2; $sum / $count" | bc)

#get median
if (( count % 2 == 0 )); then
    middle1=$(awk "NR==$count/2" ${1}.length.distribution)
    middle2=$(awk "NR==($count/2)+1" ${1}.length.distribution)
    median=$(echo "scale=2; ($middle1 + $middle2) / 2" | bc)
else
    median=$(awk "NR==($count+1)/2" ${1}.length.distribution)
fi
LengDist ${1}.length.distribution
rm ${1}.length.distribution
#write all out
echo -e "SAMPLE_ID\tMIN\tMAX\tMEAN\tMEDIAN" >> ${1}.length.stats.tsv
echo -e "${1}\t${lowest}\t${highest}\t${mean}\t${median}" >> ${1}.length.stats.tsv
}

filter()
{
awk '{print $11}' ${1}.Summary.tsv > $1.quality.distribution
#Remove title
sed -i '1d' $1.quality.distribution
#Extract fail and pass with 10 cut off
LC_NUMERIC=C awk '{if ($1>=10) print $1}' $1.quality.distribution > $1.PASS
LC_NUMERIC=C awk '{if ($1<10) print $1}' $1.quality.distribution > $1.FAIL
#QC
PASS=$(wc -l < $1.PASS)
FAIL=$(wc -l < $1.FAIL) 
TOTAL=$(wc -l < $1.quality.distribution)
#MEAN
sum=$(awk '{s+=$1} END {print s}' ${1}.quality.distribution)
count=$(wc -l < ${1}.quality.distribution)
mean=$(echo "scale=2; $sum / $count" | bc)
BasecallQC=$(echo "scale=3; $PASS/ $TOTAL" | bc)
echo -e "SAMPLE_ID \t PASS_READS \t FAIL_READS \t TOTAL_READS \t AVERAGE_Q \t ACCURACY" >> $1_QC_Stats.txt
echo -e  "${1}\t${PASS}\t${FAIL}\t${TOTAL}\t${mean}\t${BasecallQC}" >> $1_QC_Stats.txt
ReadQual $1.quality.distribution
rm $1.PASS $1.FAIL $1.quality.distribution

}

LengDist()
{
Rscript -e "
library(ggplot2)
ReadLength <- scan('$1')
ReadLengthPlot <- ggplot(data.frame(value = ReadLength), aes(x = value)) +
  geom_histogram(binwidth = 400, fill = 'lightblue', color = 'black') +
  scale_x_continuous(breaks = seq(0, max(ReadLength), by = 5000)) +
  labs(title = 'Read length distribution', x = 'Length', y = 'Read count')+
  theme(plot.title = element_text(hjust = 0.5))

png('$1.png')
print(ReadLengthPlot)
dev.off()
"
}

ReadQual()
{
Rscript -e "
library(ggplot2)
ReadQuality <-scan('$1')
ReadQualityPlot <- ggplot(data.frame(value = ReadQuality), aes(x = value)) +
  geom_histogram(binwidth = 1, fill = 'lightblue', color = 'black') +
  scale_x_continuous(breaks = seq(0, max(ReadQuality), by = 1)) +
  labs(title = 'Read quality distribution', x = 'Quality', y = 'Read count')+
  theme(plot.title = element_text(hjust = 0.5))

png('$1.png')
print(ReadQualityPlot)
dev.off()
"
}

QCV5BAM()
{
$call "$model1" "$input" > ${input}.fast_v5.0.0.bam
$call --min-qscore 10 "$model1" "$input" > ${input}.fast_v5.0.0.min_Q_score_10.bam
RawProc ${input}.fast_v5.0.0.bam
summ ${input}.fast_v5.0.0.min_Q_score_10.bam
}

QCV43BAM()
{
$call "$model2" "$input" > ${input}.fast_v4.3.0.bam
$call --min-qscore 10 "$model2" "$input" > ${input}.fast_v4.3.0.min_Q_score_10.bam
RawProc ${input}.fast_v4.3.0.bam
summ ${input}.fast_v4.3.0.min_Q_score_10.bam
}

QCHACBAM()
{
$call "$model3" "$input" > ${input}.hac_v5.0.0.bam
$call --min-qscore 10 "$model3" "$input" > ${input}.hac_v5.0.0.min_Q_score_10.bam
RawProc ${input}.hac_v5.0.0.bam
summ ${input}.hac_v5.0.0.min_Q_score_10.bam
}

RawProc()
{
summ $1
filter $1
}


fqCallV5()
{
$call --emit-fastq "$model1" "$input" > ${input}.fast_v5.0.0.fastq
$call --emit-fastq --min-qscore 10 "$model1" "$input" > ${input}.fast_v5.0.0.min_Q_score_10.fastq
}


fqCallV43()
{
$call --emit-fastq "$model2" "$input" > ${input}.fast_v4.3.0.fastq
$call --emit-fastq --min-qscore 10 "$model2" "$input" > ${input}.fast_v4.3.0.min_Q_score_10.fastq
}

fqCallHAC()
{
$call --emit-fastq "$model3" "$input" > ${input}.hac_v5.0.0.fastq
$call --emit-fastq --min-qscore 10 "$model3" "$input" > ${input}.hac_v5.0.0.min_Q_score_10.fastq
}

filterFQ()
{
#Extract headers
grep -E "^@[A-Za-z0-9_-]{5}" $1*.0.fastq > $1.nofilt.headers
awk '$1 ~ /^@[A-Za-z0-9-]+$/ && $1 !~ /[?.]/ { print $0} ' $1.nofilt.headers | sort > $1.nofilter.headers
grep -E "^@[A-Za-z0-9_-]{5}" $1*10.fastq > $1.min10.headers
awk '$1 ~ /^@[A-Za-z0-9-]+$/ && $1  !~ /[?.]/ { print $0 } ' $1.min10.headers | sort > $1.filtered.headers
#comm
comm -23 $1.nofilter.headers $1.filtered.headers > $1.uniqueID
#Extract failures
#QC
filtered_count=$(wc -l < $1.filtered.headers)
nofilter_count=$(wc -l < $1.nofilter.headers)
fail_count=$(wc -l < $1.uniqueID)
BasecallQC=$(echo "scale=3; $filtered_count / $nofilter_count" | bc)
echo -e  "${1}.fast.v5.0.0\t${filtered_count}\t${fail_count}\t${nofilter_count}\t${BasecallQC}" >> $1_Stats.txt
rm $1.filtered.headers $1.nofilter.headers $1.uniqueID $1.nofilt.headers $1.min10.headers
}



align()
{
~/Task/Software/minimap2-2.28_x64-linux/minimap2  -x map-ont -a ../Genome/GCA_001708105.1_ASM170810v1_genomic.fna $1
}



input="$1"
# Show model choices
echo "Which model would you like to use for basecalling?"
echo "1. dna_r10.4.1_e8.2_400bps_fast@v5.0.0"
echo "2. dna_r10.4.1_e8.2_400bps_fast@v4.3.0"
echo "3. dna_r10.4.1_e8.2_400bps_hac@v5.0.0"
# Select model
read -p "Enter the number (1, 2 or 3): " model_choice
#Calling

call="../Software/dorado-0.8.0-linux-x64/bin/dorado basecaller"
summary="../Software/dorado-0.8.0-linux-x64/bin/dorado summary"
model1="../Software/dorado-0.8.0-linux-x64/dna_r10.4.1_e8.2_400bps_fast@v5.0.0"
model2="../Software/dorado-0.8.0-linux-x64/dna_r10.4.1_e8.2_400bps_fast@v4.3.0"
model3="../Software/dorado-0.8.0-linux-x64/dna_r10.4.1_e8.2_400bps_hac@v5.0.0"
# Analysis based on model choice
case $model_choice in
1)
        echo "Model = dna_r10.4.1_e8.2_400bps_fast@v5.0.0."
	#fqCallV5 $1
	#filterFQ $1
	QCV5BAM $1
	#align $1.fast_v5.0.0.min_Q_score_10.fastq > $1.fast_v5.0.0.minimap2.sam
	
;;
2)
        echo "Model = dna_r10.4.1_e8.2_400bps_fast@v4.3.0."
        #fqCallV43 $1
	#filterFQ $1
	QCV43BAM $1
	#align $1.fast_v4.3.0.min_Q_score_10.fastq > $1.fast_v4.3.0.minimap2.sam
;;
3)
	echo "Model = dna_r10.4.1_e8.2_400bps_hac@v5.0.0"
	fqCallHAC $1
	filterFQ $1
	QCHACBAM
	align $1.hac_v5.0.0.min_Q_score_10.fastq > $1.hac_v5.0.0.minimap2.sam
;;
*)
    echo "Invalid choice. Please select 1, 2 or 3."
    exit 1
    ;;
esac

mkdir Pipe
mkdir Pipe/Basecalls
mkdir Pipe/Alignment
mkdir Pipe/Reports
mkdir Pipe/Graphs
mv ${input}*fastq ${input}*bam Pipe/Basecalls
mv ${input}*sam Pipe/Alignment
mv ${input}*txt ${input}*tsv Pipe/Reports
mv ${input}*png Pipe/Graphs
