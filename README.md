This vignette describes all steps necessary to run SeeCiTe analysis,
using the public HapMap trio data, with intermidiate outputs supplied
with the package.

Installation.
-------------

Make sure the dependencies are installed fist:

``` r
generic_packages <- c("magrittr", "dplyr", "tidyr", "tools", "purrr", "utils", "rlang", "bedr")
plotting_packages <- c("ggplot2", "scales", "gridExtra", "cowplot", "rogme", "ggpubr")
stat_packages <- c("statip", "outliers", "effsize", "lawstat", "ks")
packages <- c(generic_packages, plotting_packages, stat_packages)

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
```

Then use devtools package to install directly from GitHub:

``` r
library(devtools)
devtools::install_github("aksenia/SeeCiTe")
```

Step I. Preparing the input files.
----------------------------------

The preparation step takes in 1) an original PennCNV-trio output
(produced by running PennCNV’s detect\_cnv.pl with the -trio flag) and
2) merged and/or filtered by frequency and size file in a standard
PennCNV format (PennCNV’s clean\_cnv.pl will do the segment merging
automatically and output a file in such format). The merged file defines
the CNVs to analyse in terms of boundaries and loci covered, e.g. all
CNVs that are in the first file but do not overlap with the CNVs in the
second file will be ignored.

The function *runExtractInheritance* will take these two input files and
produce additional intermediate files, necessary for consequent steps,
shown below, with the prefix of the merged file:

    # Loading SeeCiTe

``` r
library(SeeCiTe)
# PennCNV-trio output
file_original <- system.file("extdata", "affy6ceu.original.triocnv", package = "SeeCiTe")
# PennCNV merge output
file_merged <- system.file("extdata", "affy6ceu.merged.filtered.triocnv", package = "SeeCiTe")
# Input files for SeeCiTe
input_files <- runExtractInheritance(filename_orig = file_original, filename_merged = file_merged)
```

The input files now contain CNVs to analyse for each offspring and
inheritance, decoded from PennCNV-trio HMM state. The merging log is
also created to keep track whether a CNV was merged and if so, how many
segments were merged.

``` r
print(input_files)
# $triocnv_file
# [1] "/Users/alpaca/Documents/uib/dev/SeeCiTe/inst/extdata/affy6ceu.merged.filtered_annot_offspring.triocnv"
# 
# $merge_trace
# [1] "/Users/alpaca/Documents/uib/dev/SeeCiTe/inst/extdata/affy6ceu.merged.filtered_merge.log"
dir <- dirname(file_original)
# Intermidiate files for the reference of inheritance mapping.
list.files(dir, pattern = tools::file_path_sans_ext(basename(file_merged)), full.names = F)
# [1] "affy6ceu.merged.filtered_annot_offspring.triocnv"
# [2] "affy6ceu.merged.filtered_annot.triocnv"          
# [3] "affy6ceu.merged.filtered_merge.log"              
# [4] "affy6ceu.merged.filtered_unmaskedstatus.triocnv" 
# [5] "affy6ceu.merged.filtered.triocnv"
```

Step II. Prepare and extract SNP data.
--------------------------------------

For the extraction of the SNP-level data for each individual in a trio
the following inputs are needed: 1) PFB file (same file with probe
coordinates used when running PennCNV); 2) full path to PennCNV
installation; 3) File with paths to LRR and BAF signal files in a
tab-separated format in the order father, mother, offspring (same as in
PennCNV-trio); 4) A parameter setting how many probes in flanks to
extract; 5) dataset name - it must be consistent and will be used for
file naming throughout the project; 6) full path to the directory in
which the extracted SNP data will be stored, for each CNV (must be
created in advance).

``` r
pfb_file <- file.path("~/Documents/uib/dev/toydata/affygw6.hg19.sorted.pfb")
penn_path <- "~/local/PennCNV1.0.4"
penn_trio_list <- file.path("~/Documents/uib/dev/toydata/affy6hm_trio.tab")
n_flanking_snp <- 5
run_dir <- "~/Documents/uib/dev/toydata/dev"

commands <- makePythonCommands(penn_path=penn_path, 
                               pfb_file=pfb_file, 
                               penn_trio_list=penn_trio_list, 
                               triocnv_file=input_files[["triocnv_file"]],
                               n_flanking_snp=5, 
                               dataset="affy6ceu", 
                               run_dir=run_dir)
print(commands)
# $cnv
# [1] "python3 /Users/alpaca/Documents/uib/dev/SeeCiTe/inst/python/extract_snp_cnv.py -l ~/Documents/uib/dev/toydata/affy6hm_trio.tab -c /Users/alpaca/Documents/uib/dev/SeeCiTe/inst/extdata/affy6ceu.merged.filtered_annot_offspring.triocnv -d affy6ceu -p ~/Documents/uib/dev/toydata/affygw6.hg19.sorted.pfb -s ~/local/PennCNV1.0.4 -o ~/Documents/uib/dev/toydata/dev"
# 
# $flank
# [1] "python3 /Users/alpaca/Documents/uib/dev/SeeCiTe/inst/python/extract_snp_flanks.py -l ~/Documents/uib/dev/toydata/affy6hm_trio.tab -c /Users/alpaca/Documents/uib/dev/SeeCiTe/inst/extdata/affy6ceu.merged.filtered_annot_offspring.triocnv -d affy6ceu -p ~/Documents/uib/dev/toydata/affygw6.hg19.sorted.pfb -s ~/local/PennCNV1.0.4 -o ~/Documents/uib/dev/toydata/dev -f 5"
```

The result will be two script files with one line per CNV with a command
for PennCNV infer\_snp\_allele.pl that will do the extraction,
terminating by lines that collect the data into one table for the whole
cohort. This must be run by the user. For large cohorts one can split
the commands into batches or submit to a cluster.

Step III. Gather and read in all input.
---------------------------------------

The previous step extracts SNP data into files in the provided
*run\_dir*: files with prefix *probecoord.txt*, *snp\_flank.txt* and
*snp\_cnv.log*. The main CNV file is *triocnv\_file* in *input\_data*
object, while *merge\_trace* is the merging log in the same object.
Finally, *cnv\_qcsum\_file* is the QC summary output of PennCNV. The
*cache\_id* tells where R should store cache for core calculations.

``` r
args <- list(triocnv_file=input_files[["triocnv_file"]],
             probecoord_file=system.file("extdata", "affy6ceu.probecoord.txt", package = "SeeCiTe"),
             snp_flank_file=system.file("extdata", "affy6ceu.snp_flank.txt", package = "SeeCiTe"),
             snp_cnv_log_file=system.file("extdata", "affy6ceu.snp_cnv.log", package = "SeeCiTe"),
             cnv_qcsum_file=system.file("extdata", "affy6ceu.qcsum", package = "SeeCiTe"),
             dataset="affy6ceu",
             cache_id="~/Documents/uib/dev/toydata",
             merge_trace=input_files[["merge_trace"]])
```

Now all inputs are in order and can be read and formatted:

``` r
main_dt <- readInputs(args = args)
# [1] "Reading merged formatted PennCNV trio file for offspring /Users/alpaca/Documents/uib/dev/SeeCiTe/inst/extdata/affy6ceu.merged.filtered_annot_offspring.triocnv"
# [1] "Reading probes file /Users/alpaca/Documents/uib/dev/SeeCiTe/inst/extdata/affy6ceu.probecoord.txt"
# [1] "Collapsing duplicate probes"
# [1] "Annotating probes in CNV loci"
# [1] "Reading extracted intensity for flanking regions /Users/alpaca/Documents/uib/dev/SeeCiTe/inst/extdata/affy6ceu.snp_flank.txt"
# 
#    denovo inherited 
#         6        43 
# [1] "Reading denovo test results for CNV regions /Users/alpaca/Documents/uib/dev/SeeCiTe/inst/extdata/affy6ceu.snp_cnv.log"
# [1] "Reading quality summary file /Users/alpaca/Documents/uib/dev/SeeCiTe/inst/extdata/affy6ceu.qcsum"
# [1] "Reading merging log file /Users/alpaca/Documents/uib/dev/SeeCiTe/inst/extdata/affy6ceu.merged.filtered_merge.log"
candidateCnvs <- main_dt[["data"]]
```

Step IV. Run SeeCiTe classification.
------------------------------------

First, a summary statistic collection step, for each CNV in offspring:

``` r
clu_baf <- runAnalyzeSignal(input_data = candidateCnvs, args = args, use_cache = T)
# [1] "Found cached file ~/Documents/uib/dev/toydata/affy6ceu_clu_baf.rds, reading in"
head(clu_baf[,c(1:4)], n=3)
# # A tibble: 3 x 4
#   cnvTypeBAF    relation nprobes33_66 nprobes_cn4
#   <chr>         <chr>           <int>       <int>
# 1 F:possiblyLOH F                   0           0
# 2 M:possiblyLOH M                   0           0
# 3 O:possiblyLOH O                   0           0
```

The clasification is the final step in the analysis which annotates each
CNV with suggested inheritance and SeeCiTe quality class.

``` r
cnv_class <- classifyTrios(clu_baf)
with(cnv_class, table(seecite, inheritanceTest))
#             inheritanceTest
# seecite      denovo inherited mosaic unclear
#   borderline      3         0      3      19
#   probable        0        18      0       0
#   unlikely        3         3      0       0
```

Step V. Visualize and write summary files.
------------------------------------------

The results can be visualized either for each single CNV region or for a
whole cohort:

``` r
Sample <- "affy6.scale.NA12878"
Cnv <- "chr19:20596206-20716389"
plotRawTrio(input_data = candidateCnvs %>% dplyr::filter(sample==Sample, coordcnv==Cnv), 
            sifted_data = clu_baf %>% dplyr::filter(sample==Sample, coordcnv==Cnv), 
            penn_qcsum = main_dt[["qcsum"]] %>% dplyr::filter(sample==Sample),
            merge_trace = main_dt[["merge"]] %>% dplyr::filter(sample==Sample, coordcnv==Cnv))
```

![](man/figures/step-v-indiv-plot-1.png)

This will write a pdf file with such plots per SeeCiTe category:

``` r
plotCohort(main_data=main_dt,
           sifted_data=clu_baf,
           classified_data=cnv_class,
           output_dir = "~/Documents/uib/dev/toydata/affy6ceu_viz",
           dataset="affy6ceu",
           subset_nprobes=20,
           subset_length=150000)
```

Finally, the summary statistics and SeeCiTe classifications can be
written out as plain text files, together with bed (UCSC 6-column style)
and plink formatted CNV regions:

``` r
writeSeecite(classified_data=cnv_class,
          output_dir = "~/Documents/uib/dev/toydata/affy6ceu_viz",
          dataset="affy6ceu")
```
