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
import argparse
import re

import pandas as pd

from collections import defaultdict

def write_list(filename, listvar):
    '''
    Write list to a file, one word per row
    Args:
        filename: name of the file out
        listvar: list of strings

    Returns: NA, writes out a file

    '''
    out = open(filename, "w")
    for s in listvar:
        out.write("%s\n" % s)
    out.close()

def arrayToDF(chr_array, idx_array):
  '''
  Given an array of probe coordinates for a chromosom,
  A selected index array holding all SNPs to extract,
  Extract the LRR and BAF values from the chromosome array and
  format it as a pandas data frame
  :param chr_array:
  :param idx_array:
  :return:
  '''
  df = pd.DataFrame([])
  for i in idx_array:
      j = [chr_array[str(i)]]
      dfj = pd.DataFrame(j, columns = ['Name', 'Chr','Position', 'LRR', 'BAF'])
      df = df.append(dfj)
  return (df)


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
  # define an index array holding CNV snps
  cnvarray = probearray[probearray.index(int(start)):(probearray.index(int(end)) + 1)]
  cnvdf = arrayToDF(chr_array=chr_array, idx_array=cnvarray)
  cnvdf['locus'] = 'CNV'
  # snps upstream the CNV    
  larray = probearray[:probearray.index(int(start))]
  # if there is less SNPs left than the chosen flank size
  if (len(larray)==0):
    leftprobe=chr_array[start]
  else:
    idx = min(flank, len(larray))
    leftprobe = chr_array[str(larray[-idx])]
  # extract the index snps for left flank
  leftflank = probearray[probearray.index(int(leftprobe[2])):probearray.index(int(start))]
  # format it to pandas data frame
  leftdf = arrayToDF(chr_array=chr_array, idx_array=leftflank)
  # tag the locus
  leftdf['locus'] = 'flank'
  # snp downstrem the CNV
  rarray = probearray[(probearray.index(int(end)) + 1):]
  if (len(rarray)==0):
    rightprobe=chr_array[end]
  else:
    idx = min(flank, len(rarray))
    rightprobe = chr_array[str(rarray[idx-1])]
  rightflank = probearray[(probearray.index(int(end)) + 1):(probearray.index(int(rightprobe[2]))+1)]
  rightdf = arrayToDF(chr_array=chr_array, idx_array=rightflank)
  rightdf['locus'] = 'flank'
  
  frames = [leftdf, cnvdf, rightdf]
  result = pd.concat(frames)
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
  idx = 0
  for line in probe_lines:
      s=line.strip().split('\t')
      if idx > 0:
          probe_dict[s[1]].append(s)
      idx+=1
  return(probe_dict)

def readPenn(penn_file):
  '''
  Read the cnv file from PennCNV after individual calling or merging
  format as dictionary of dictionary by sample id and then cnv coord
  :param penn_file: Path to a file
  :return:
  '''
  with open(penn_file) as file:
    dcnvs = file.read().splitlines()

  indivs = defaultdict(list)
  for line in dcnvs:
      s=line.strip().split('\t')
      s[:] = [item for item in s if item != '']
      indivs[s[4]].append(s)
  cnvs = defaultdict(list)
  for key in indivs.keys():
      cnv = defaultdict(list)
      for item in indivs[key]:
          cnv[item[0]].append(item)
      cnvs[key].append(cnv)
  return(cnvs)
      
def extractSNPdata(cnv_dict, signal_file, flank):
  '''
  Given a subset of Penn CNV file for one individual, formatted as dict,
  a corresponding sample signal file and
  a number of snps in flanking regions,
  extract the corresponding raw signal data from the signal file
  for each CNV in the Penn file
  :param cnv_file:
  :param signal_file:
  :param flank:
  :return:
  '''
  
  # read signal file to a dictionary
  print("Reading in signal file " + signal_file)
  pfb = readPfb(pfb_file = signal_file)

  penni = cnv_dict
  # parse each line for cnv and extract the coresponding SNP data
  cnvs = []
  print("Parsing CNV lines")
  for k in penni.keys():
    l = penni[k][0]
    chrom = k.split(':')[0]
    cnvstart, cnvend = k.split(':')[1].split('-')
    chrom = re.sub('[a-z]', '', chrom)
    chrarr = defaultdict(list)
    for p in pfb[chrom]:
      chrarr[p[2]] = p
    individual = l[4]
    copynum = re.sub('.*=', '', l[3])
    numsnp = re.sub('.*=', '', l[1])
    print("Individual " + individual + " CNV " + k)
    r = getFlankingProbes(chr_array=chrarr, start=cnvstart, end=cnvend, flank=flank)
    r['coordcnv'] = k
    r['sample'] = individual
    r['copynumber'] = copynum
    r['numsnp'] = numsnp
    cnvs.append(r)
  result = pd.concat(cnvs)
  return(result)


def runExtractSNPdata(args):
  '''
  Given an array of script arguments: cnv for PennCNV file; map for comma-separated mapping of
  samples and the corresponding raw data file (relative of full path);
  flank_size - number of probes to extract from each side of CNV
  Run the extract SNP function for each individual and pool the result in one data frame
  :param args:
  :return:
  '''
  # file with CNV calls 
  cnv_file = args['cnv']
  # read map file
  map_file = args['map']
  with open(map_file) as mp:
    map_lines = mp.read().splitlines()
  smap = {}
  for l in map_lines:
      s = l.split(',')
      smap[s[0]] = s[1]
  # read number of probes in flanks to extract
  flank=int(args['flank_size'])
  
 # read CNV calls
  pennd = readPenn(cnv_file)
  outdata = []
  for key in smap.keys():
      print("Processing sample " + key)
      cnv_indiv = pennd[key][0]
      pfb_file = smap[key]
      sig = extractSNPdata(cnv_dict = cnv_indiv, signal_file = pfb_file, flank = flank)
      outdata.append(sig)
  total_out = pd.concat(outdata)
  return(total_out)

def writeOutput(args, outdata):
  '''
  Write the resulted pandas data frame with data into a flat file
  :param args:
  :return:
  '''
  outfile = args['out'] + "/" + args['dataset'] + ".signal_flanks_" + args['flank_size'] + ".txt"
  print("Done! Writing the results into file " + outfile)
  outdata.to_csv(outfile, index = False, header=True, sep='\t')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Extract data for SNPs in regions flanking CNV')
  parser.add_argument('-c','--cnv', nargs = '?', help='Full path to the list of CNV calls in PennCNV format', required=True)
  parser.add_argument('-m','--map', nargs = '?', help='Full path to map file in format sampleid,path/to/signalfile; one per line', required=True)
  parser.add_argument('-o','--out', nargs = '?', help='Full path to output folder ', required=True)
  parser.add_argument('-d','--dataset', nargs = '?', help='A string for dataset name used for outputs', required=True)
  parser.add_argument('-f','--flank_size', nargs = '?', help='Number of probes to include in each flank, default 50 ', default=50, required=False)

  args = vars(parser.parse_args())
  print("Arguments supplied: ")
  print(args)

  if not os.path.exists(args['out']):
      print("Creating output directory")
      os.makedirs(args['out'])
  print("Extracting the SNP data...")
  result = runExtractSNPdata(args = args)
  writeOutput(args = args, outdata = result)
                      
                    

