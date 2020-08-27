#!/usr/bin/env bash

# Usage info
show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-d DATASET] [-s EXTRACT_SNP_DIR] [-p PFB]
This script must be run after the SNP extraction commands with infer_snp_allele.pl has been run. It does convenient parsing of the outputs and logs preparing for the next step.
       -h                     Display this help and exit
       -d DATASET             String denoting dataset name. Must be consistent throughout the project
       -s EXTRACT_SNP_DIR     Full path to extracted SNP data folder (after infer_snp_allele.pl commands have been run)
       -p PFB                 Full path to PennCNV PFB file


EOF
}

  unset BATCH
  OPTIND=1 # Reset is necessary if getopts was used previously in the script.
  while getopts "hd:s:p:" opt; do
       case "$opt" in
           h)
               show_help
               exit 0
               ;;
           d)  DATASET=$OPTARG
               ;;
           s)  EXTRACT_SNP_DIR=$OPTARG
               ;;
           p)  PFB=$OPTARG
               ;;
           '?')
               show_help >&2
               exit 1
               ;;
       esac
   done
   shift "$((OPTIND-1))" # Shift off the options and optional --.

PFB=$1
PENN=$2
FILE=$3


# collapse state1 and 2, state5 and 6
# also append family relation to sample id for convenience
# this is needed for downstream analysis
# since infer_snp_allele only understands copynumber 1 and 3

PREFIX=${FILE%.*}
BASE=$(basename $PREFIX)
dir_name=$(dirname $FILE)


merge_dir=${dir_name}/merge_${BASE}
mkdir -p $merge_dir
# make sure we do not overrun anything
rm $merge_dir/*
cp $FILE $merge_dir/original.triocnv


# file for merging
f=$merge_dir/original.triocnv
# initial gap fraction 50%
PERCBP=50;
DIFF=100
i=0
while [ $DIFF -gt 0 ];
do
#
echo $DIFF
pfix=${f%.*};
base=$(basename $f)
sfix=${base#*.}
NROW_IN=$(wc -l $f | awk '{print $1}')

out=${pfix}_${i}.${sfix};
log=${pfix}_${i}.log
echo "Merging iteration ${i}";
echo "file in: ${pfix}.${sfix}"
echo "file out: ${out}"
# translate percentage into the fraction
# FRACBP="0"$(echo "scale=1;${PERCBP} / 100" | bc -l)
FRACBP=$(awk -v n=$PERCBP 'BEGIN{print(n/100)}')

if (( $i % 2 )); then
  # start with conservative mind merging on markers
  ${PENN}/clean_cnv.pl --signalfile ${PFB} --fraction ${FRACBP} --output ${out} combineseg ${pfix}.${sfix} 2>&1 | tee ${log}
else
  # do most merging on long distance
  # max 50%
  ${PENN}/clean_cnv.pl --signalfile ${PFB} --fraction ${FRACBP} --bp --output ${out} combineseg ${pfix}.${sfix} 2>&1 | tee ${log}
  # decrease the gap fraction in consequent merge calls
  # to avoid overmerging
  if [ $PERCBP > 20 ]; then
  PERCBP=$(($PERCBP-5))
  fi
fi
f=$out
NROW_OUT=$(wc -l $out | awk '{print $1}')
DIFF=$(($NROW_IN-$NROW_OUT))
i=$(($i+1))
pfix=${f%.*}
sfix=${f#*.}
done

cp $f  ${merge_dir}/merged.triocnv
rm ${merge_dir}/original*.triocnv
cp ${merge_dir}/merged.triocnv ${dir_name}/${BASE}_merged.triocnv
## when doing bed files intersect
## we must separate the types!
# masked
${PENN}/visualize_cnv.pl -format bed -out ${merge_dir}/${prefix}_maskedstatus.merged.triocnv.bed <(awk -v OFS="\t" '{k=$4;gsub(".*=","",k); $5=$5k}1' ${merge_dir}/${prefix}_maskedstatus.merged.triocnv)
# unmasked: make sure that the sample id is constructed in the same way!
${PENN}/visualize_cnv.pl -format bed \
-out ${prefix}_unmaskedstatus.triocnv.bed  <(awk -v OFS="\t" '{gsub(".*\\/[a-z]+([0-9])*(_|_.|.)", "", $5); gsub(".*=","",$NF); k=$4;gsub(".*=","",k); $5=$5k":"($NF)}1' ${prefix}_unmaskedstatus.triocnv)


# intersect the two files
# for each segment in masked file (the most merged)
# collect the segments and status from unmasked file
## REMEMBER bed files are 0-based!!!
# need to add 1 to start when switching back... (line 2, field $2 for masked and $11 for unmasked)

## we merge the status as follows:
#if there are more than one CNV, all with the same status - take that status
# otherwise indicate it as ambiguous
# we only intersect events of the same type (first two elements of b and a arrays)
bedtools intersect -loj -a ${prefix}_maskedstatus.merged.triocnv.bed -b ${prefix}_unmaskedstatus.triocnv.bed |\
  awk '{split($4,b,":"); split($13,a,":"); if (b[1]":"b[2]==a[1]":"a[2]) print $0}' > ${prefix}_maskedstatus.merging.pbed

  awk -v OFS="\t" '{split($13, a, ":"); $2=$2+1; $11=$11+1; print $1":"$2"-"$3, $5, $4, $10":"$11"-"$12, $14, a[3]}' ${prefix}_maskedstatus.merging.pbed |\
  awk '{k=$1"__"$2"__"$3; a[k]=a[k]"|"$4"__"$5"__"$6}END{ for (i in a) print i, a[i]}' |\
  awk '{gsub("^\\|", "", $2)}1' |\
  awk '{l=split($2,b, "|");
  if (l==1) {len=1} else {
  split("", arr, ":");
  for (i in b) {
    split(b[i],c,"__"); st=c[3]; arr[st]=arr[st]+1};
    len=length(arr);} print $0, len, l}' |\
  awk '{if ($3==1){
    k=$2; gsub("\\|.*","",k);split(k,a,"__");
    s=a[3]=="denovo-inherited"?"ambiguous":a[3]
  } else {s="ambiguous"}; print $0,s}' > ${prefix}_maskedstatus.merge.log


# make bed12 track! start positions must be sorted numerically! the starts are relative to main start, not absolute
sort -n -k1,2 -k11 ${prefix}_maskedstatus.merging.pbed |\
  awk -v OFS="\t" '{split($13, a, ":"); n=$1":"$2"-"$3; s=$12-$11; o=$11-$2;print $1,$2,$3, n, $4, $5, $6, $9, o, s, a[2]}' |\
  awk '{k=$4"__"$5"__"$6"__"$8; a[k]=a[k]","$9; b[k]=b[k]","$10; c[k]=c[k]+1}END{ for (i in a) print i, a[i],b[i],c[i]}' |\
  awk -v OFS="\t" '{split($1, a, "__"); split(a[1], b, ":"); split(b[2], c, "-"); gsub("^,","",$2); gsub("^,","",$3); print b[1], c[1], c[2], a[1]"__"a[2], a[3], ".", c[1], c[2], a[4], $4, $3, $2}' > ${prefix}_maskedstatus.merging_ucsc.bed

awk '$10>1' ${prefix}_maskedstatus.merging_ucsc.bed > ${prefix}_maskedstatus.merging_ucsc_mult.bed


### track the PEnnCNV trio status through merging
bedtools intersect -loj -a <(awk '$4~"offspring"' ${prefix}_maskedstatus.merged.triocnv.bed) -b ${dir_name}/${base}_pennstatus.bed |\
awk '{split($4,a,":"); split($13,b,","); f=a[1]==b[1]?1:0}f' |\
awk '{split($13, b, ","); id=$1":"($2+1)"-"$3"__"$4; a[id]=a[id]"|"b[2]; c[id]=c[id]+1}END{for (i in a) {s=a[i]; cc=c[i]; gsub("^\\|","",s); gsub("__", " ", i); print i, s,cc}}' > ${dir_name}/${base}_pennstatus_merge.log

awk 'FNR==NR{arr[$1]=$2;next}{l=$1"__"$5; $5=$5":"arr[l]; print $0}' <(awk '{split($1,a,"__"); k=a[1]"__"a[3]; gsub("[0-9]$", "", k); gsub("[0-9]$", "", $NF); print k, $NF}' ${prefix}_maskedstatus_merging.log) ${prefix}_maskedstatus_merged.triocnv > ${prefix}_masked_merged_annot.triocnv
# final merged file with statuses for all individuals in the trio
cp ${prefix}_masked_merged_annot.triocnv ${dir_name}/${base}_trio.triocnv
# only offspring, as downstream infer_snp_allele is applied to calls in offspring
awk '$5~"offspring"' ${dir_name}/${base}_trio.triocnv > ${dir_name}/${base}.triocnv
#awk '{split($5, b, ":"); a[b[1]]++}END{for (i in a) print i, a[i]}' ${prefix}_merged.triocnv > ${prefix}_merged_NumCNV.txt



#awk 'NR==FNR{b[$1]=$2;next}{split($5, a, ":"); if (b[a[1]]<46) gsub(".*\\/[a-z]+([0-9])*(_|_.|.)", "", $5); split($5,d,":"); $5=d[1]":"d[2]":triostate="d[3]; print $0}' ${prefix}_merged_NumCNV.txt  ${prefix}_merged.triocnv > ${prefix}_merged_filterN.triocnv

