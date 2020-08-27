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


OUT=${EXTRACT_SNP_DIR}/${DATASET}_candidateCnv.txt
OUTLOG=${EXTRACT_SNP_DIR}/${DATASET}_candidateCnvLOG.txt

touch $OUT
for f in $(find ${EXTRACT_SNP_DIR}/${DATASET}*.out | grep -v runcommands | grep -v log | grep -v candidateCnv); do
  b=$(basename ${f});
  echo $b;
  awk -v x=$b 'NR==1{next}{print $0, x}' $f >> $OUT;
done

touch $OUTLOG
for f in $(find ${EXTRACT_SNP_DIR}/${DATASET}_* | grep log | grep -v candidateCnv); do
  b=$(basename ${f});
  px=${b%.*};
  echo "-------";
  cat $f | grep "vidence" | awk -v x=$px -v OFS=',' '{print $0, x}'  >> $OUTLOG;
done

# rm ${EXTRACT_SNP_DIR}/*.out
# rm ${EXTRACT_SNP_DIR}/*.log

awk 'NR==FNR {id[$1]; next} $1 in id' <(cut -f1 ${OUT}) $PFB > ${EXTRACT_SNP_DIR}/${DATASET}_probecoord.txt

