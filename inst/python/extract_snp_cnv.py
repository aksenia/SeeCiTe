#!/usr/bin/env python3
'''
Prepare the commands for extracting SNP LRR and BAF values with Penn CNV script infer_snp_allele.pl
This script must be run after the trio CNV calling is done
It uses predefined project structure and filenames
'''

import json
import sys
import os
import subprocess
import argparse

from datastructutils import write_list

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

def createExtractCommands(args):
  '''
  Create the commands to run PennCNV infer_snp_allele script to extract the probe data
  for the list of CNVs
  :param args:
  :return:
  '''

  # read CNV calls
  cnv_file = args['cnv']
  pfb_file = args['pfb']
  sample_file = args['list']
  dataset = args['dataset']
  penn = args['penn']

  command_prefix = penn + '/infer_snp_allele.pl -pfb ' + pfb_file + ' -hmm ' + penn + '/lib/hhall.hmm -denovocn '

  with open(cnv_file) as file:
    dcnvs = file.read().splitlines()
  # output folder
  cdir = args['out']
  # read samples in trios
  samples = readSampleList(sample_file = sample_file)
  cnvs = []
  for line in dcnvs:
    l = line.rstrip('\n').split('\t')
    coordcnv = l[0]
    offspring = l[4].split(":")[0]
    copynumber = l[3].split("=")[1]
    startsnp = l[5].split("=")[1]
    endsnp = l[6].split("=")[1]
    outf = cdir + "/" + dataset + "_" + offspring + "_" + coordcnv
    cmnd = command_prefix + copynumber + ' '  + ' -start ' + startsnp + ' -end ' + endsnp + ' -out ' + outf + '.out' +  ' -log ' + outf + '.log ' + samples[offspring]
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
  outd=args['out'] + "/"
  out=outd + args['dataset'] + '.snp_cnv.txt'
  outlog=outd + args['dataset'] + '.snp_cnvLOG.txt'
  cmd_lines.append('##############################')
  cmd_lines.append('if [  -f  ' + out + ' ]; then')
  cmd_lines.append('echo "File ' + out + ' exists. Removing the file!" ')
  cmd_lines.append('rm ' + out + ';')
  cmd_lines.append('fi')
  cmd_lines.append('touch ' + out)
  cmd_lines.append('for f in $(find ' + outd + args['dataset']+ '*.out | grep -v runcommands | grep -v log | grep -v snp_cnv); do')
  cmd_lines.append('b=$(basename ${f});')
  cmd_lines.append('pf=${b%.*};')
  cmd_lines.append('awk -v x=$pf \'NR==1{next}{print $0, x}\' $f >> ' + out + ';done')
  cmd_lines.append('##############################')
  cmd_lines.append('if [  -f  ' + outlog + ' ]; then')
  cmd_lines.append('echo "File ' + outlog + ' exists. Removing the file!" ')
  cmd_lines.append('rm ' + outlog + ';')
  cmd_lines.append('fi')
  cmd_lines.append('touch ' + outlog)
  cmd_lines.append('for f in $(find ' + outd + args['dataset']+ '*.log | grep log | grep -v snp_cnv); do')
  cmd_lines.append('b=$(basename ${f});')
  cmd_lines.append('px=${b%.*};')
  cmd_lines.append('cat $f | grep "vidence" | awk -v x=$px -v OFS=\',\' \'{print $0, x}\'  >> ' + outlog + ';done')
  cmd_lines.append('##############################')
  cmd_lines.append('# uncomment for clean up')
  cmd_lines.append('#rm ' + outd + '*.out')
  cmd_lines.append('#rm ' + outd + '*.log')
  return(cmd_lines)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Create commands for SNP extraction for CNV region')
  parser.add_argument('-l','--list', nargs = '?', help='Full path to the list of trio samples, tab-separated, ordered father mother offspring, one trio per line ', required=True)
  parser.add_argument('-c','--cnv', nargs = '?', help='Full path to the list of CNV calls in offspring only in custom format ', required=True)
  parser.add_argument('-d','--dataset', nargs = '?', help='Dataset id ', required=True)
  parser.add_argument('-p','--pfb', nargs = '?', help='Full path to dataset-specific PFB file ', required=True)
  parser.add_argument('-s','--penn', nargs = '?', help='Full path to PennCNV installation ', required=True)
  parser.add_argument('-o','--out', nargs = '?', help='Full path to output folder ', required=True)

  args = vars(parser.parse_args())
  print("Arguments supplied: ")
  print(args)

  if not os.path.exists(args['out']):
      os.makedirs(args['out'])

  command_file = args['out'] + "/" + args['dataset'] +"_runcommands.txt"
  cmd_lines = parseExtractedOutput(args=args)
  cnvs = createExtractCommands(args=args)
  print("Writing file " + command_file)
  write_list(filename = command_file , listvar = cnvs + cmd_lines)
