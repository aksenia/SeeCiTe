#!/usr/bin/env python3
'''
This script is to extract the signal information around the CNV calls
For analysis and visualisation - the context of a call
This script must be run after the trio CNV calling is done
It uses predefined project structure and filenames
'''

import json
import sys
import os
import subprocess
import argparse
import re

from datastructutils import write_list
from collections import defaultdict


def getFlankingProbes(chr_array, start, end, flank):
  '''
  Given an array of probe coordinates for a chromosome
  A start and an end of a region, and
  Length of desired flanks, find the probes closest to the flanking region
  :param chr_array:
  :param start:
  :param end:
  :param flank:
  :return:
  '''
  # extract the actual coordinates
  probearray = [int(x) for x in chr_array.keys()]
  probearray.sort()

  leftflank = int(start) - flank
  rightflank = int(end) + flank

  larray = probearray[:probearray.index(int(start))]
  if (len(larray)==0):
    leftprobe=chr_array[start]
  else:
    idx = min(flank, len(larray))
    leftprobe = chr_array[str(larray[-idx])]

  rarray = probearray[(probearray.index(int(end)) + 1):]
  if (len(rarray)==0):
    rightprobe=chr_array[end]
  else:
    idx = min(flank, len(rarray))
    rightprobe = chr_array[str(rarray[idx-1])]
  result = [leftprobe, rightprobe]
  return (result)

def readPfb(pfb_file):
  '''
  Read the PFB file into an ordered dictionary
  :param pfb_file: Dataset-specific PFB file with probe coordinates
  :return:
  '''
  with open(pfb_file) as file:
    probe_lines = file.read().splitlines()
  # make dictionary for coordinates
  probe_dict = defaultdict(list)
  # store coordinates per chr
  for line in probe_lines:
      s=line.strip().split('\t')
      probe_dict[s[1]].append(s)
  return(probe_dict)

def readSampleList(sample_file):
  '''
  Make a dictionary of trios with offspring as id
  :param sample_file: List file of trio samples, one per line in order father mother offspring
  :return:
  '''
  with open(sample_file) as file:
    sample_lines = file.read().splitlines()

  samples = {}
  for line in sample_lines:
    l = line.rstrip('\n').split('\t')
    prefix=os.path.dirname(l[2])
    offspring=os.path.basename(l[2])
    samples[offspring] = line.replace('\t', ' ')
  return(samples)

def createExtractCommandsFlank(args):
  '''
  Given an array of probe coordinates for a chromosome
  A start and an end of a region, and
  Length of desired flanks, identify the probes at the boundary of the given flank
  Create the commands to run PennCNV infer_snp_allele script to extract the probe data
  :param args:
  :param flank:

  :return:
  '''

  # read CNV calls
  cnv_file = args['cnv']
  pfb_file = args['pfb']
  sample_file = args['list']
  dataset = args['dataset']
  penn = args['penn']
  flank=int(args['flank_size'])

  command_prefix = penn + '/infer_snp_allele.pl -pfb ' + pfb_file + ' -hmm ' + penn + '/lib/hhall.hmm -denovocn '

  with open(cnv_file) as file:
    dcnvs = file.read().splitlines()
  # output folder
  cdir = args['out']
  # read probe coordinates file to a dictionary
  d = readPfb(pfb_file = pfb_file)
  # read samples in trios
  samples = readSampleList(sample_file = sample_file)

  # parse each line for cnv and create corresponding command
  cnvs = []
  print("Parsing CNV lines")
  for line in dcnvs:
    l = line.rstrip('\n').split('\t')
    coordcnv = l[0]
    chr,coord = coordcnv.split(":")
    chr = re.sub('[a-z]', '', chr)
    chrarr = defaultdict(list)
    for p in d[chr]:
      chrarr[p[2]] = p[:-1]
    cnvstart, cnvend = coord.split("-")
    offspring = l[4].split(":")[0]
    copynumber = l[3].split("=")[1]
    startsnp = l[5].split("=")[1]
    endsnp = l[6].split("=")[1]
    print("Offspring " + offspring + " CNV " + coordcnv)
    r = getFlankingProbes(chr_array=chrarr, start=cnvstart, end=cnvend, flank=flank)
    outf = cdir + "/flanks/" + dataset + "_flanks_" + offspring + "_" + coordcnv
    cmnd = command_prefix + copynumber + ' '  + ' -start ' + r[0][0] + ' -end ' + r[1][0] + ' -out ' + outf + ".out" + ' -log ' + outf + ".log " + samples[offspring]
    cnvs.append(cmnd)
  return(cnvs)


def parseExtractedOutput(args):
  '''
  Given an array of probe coordinates for a chromosome
  A start and an end of a region, and
  Length of desired flanks, identify the probes at the boundary of the given flank
  Create the commands to run PennCNV infer_snp_allele script to extract the probe data
  :param args:
  :return:
  '''
  cmd_lines=[]
  outd=args['out'] + '/flanks/'
  out=outd + args['dataset'] + '.snp_flank.txt'
  outlog=outd + args['dataset'] + '.snp_flankLOG.txt'
  cmd_lines.append('touch ' + out)
  cmd_lines.append('for f in $(find ' + outd + args['dataset']+ '*.out | grep -v runcommands | grep -v log | grep -v snp_flank); do')
  cmd_lines.append('b=$(basename ${f});')
  cmd_lines.append('awk -v x=$b \'NR==1{next}{print $0, x}\' $f >> ' + out + ';done')
  cmd_lines.append('touch ' + outlog)
  cmd_lines.append('for f in $(find ' + outd + args['dataset']+ '*.out | grep log | grep -v snp_flank); do')
  cmd_lines.append('b=$(basename ${f});')
  cmd_lines.append('px=${b%.*};')
  cmd_lines.append('cat $f | grep "vidence" | awk -v x=$px -v OFS=\',\' \'{print $0, x}\'  >> ' + outlog + ';done')
  cmd_lines.append('rm ' + outd + '*.out')
  cmd_lines.append('rm ' + outd + '*.log')
  cmd_lines.append('awk \'NR==FNR{id[$1]; next} $1 in id\' <(cut -f1 ' + out + ') ' + args['pfb'] + ' > ' + outd + args['dataset'] + '_probecoord.txt')
  return(cmd_lines)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Create commands for SNP extraction for regions flanking CNV')
  parser.add_argument('-l','--list', nargs = '?', help='Full path to the list of family samples, tab-separated, ordered father mother offspring ', required=True)
  parser.add_argument('-c','--cnv', nargs = '?', help='Full path to the list of CNV calls in offspring only in custom format', required=True)
  parser.add_argument('-d','--dataset', nargs = '?', help='Dataset id', required=True)
  parser.add_argument('-p','--pfb', nargs = '?', help='Full path to dataset-specific PFB file ', required=True)
  parser.add_argument('-s','--penn', nargs = '?', help='Full path to PennCNV installation ', required=True)
  parser.add_argument('-o','--out', nargs = '?', help='Full path to output folder ', required=True)
  parser.add_argument('-f','--flank_size', nargs = '?', help='Number of probes to include in each flank, default 60 ', default=60, required=False)

  args = vars(parser.parse_args())
  print("Arguments supplied: ")
  print(args)

  if not os.path.exists(args['out'] + "/flanks"):
      os.makedirs(args['out'] + "/flanks")

  command_file = args['out'] + "/flanks/" + args['dataset'] +"_runcommands_FLANKS.txt"
  cmd_lines = parseExtractedOutput(args=args)
  cnvs = createExtractCommandsFlank(args=args)
  print("Writing file " + command_file)
  write_list(filename = command_file , listvar = cnvs + cmd_lines)

