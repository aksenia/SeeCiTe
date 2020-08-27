#' Create commands to run helper python scripts that automate SNP extraction via PennCNV infer_snp_allele.pl for a cohort.
#'
#' @param penn_path Full path to PennCNV installation folder
#' @param pfb_file Full path to PFB file used in the analysis
#' @param penn_trio_list Full path to tab-separated file, listing one trio per line in order father,mother,offspring with full paths to raw data files for each, as used in PennCNV trio
#' @param triocnv_file Full path to annotated trio file for offspring (after inheritance map)
#' @param n_flanking_snp Integer. A number of probes to extract in upstream and downstream flanks of CNV. Default is 60.
#' @param dataset Dataset string used throughout the analysis. Must be consistent. Default is TEST.
#' @param run_dir Full path to where extract the snp data. Will be incorporated into extract commands.
#'
#' @return A data table with merged sample and CNV ids as well as original tristates and counts of segments merged
#' @examples
makePythonCommands <- function(penn_path, pfb_file, penn_trio_list, triocnv_file, n_flanking_snp=60, dataset="TEST", run_dir){
  helper_script_cnv <- system.file("python", "extract_snp_cnv.py", package = "SeeCiTe")
  helper_script_flank <- system.file("python", "extract_snp_flanks.py", package = "SeeCiTe")
  command_flank <-paste0("python3 ", helper_script_flank,
                         " -l ", penn_trio_list,
                         " -c ", triocnv_file,
                         " -d ", dataset,
                         " -p ", pfb_file ,
                         " -s ", penn_path,
                         " -o ", run_dir,
                         " -f ", n_flanking_snp)
  command_cnv <- paste0("python3 ", helper_script_cnv,
                        " -l ", penn_trio_list,
                        " -c ", triocnv_file,
                        " -d ", dataset,
                        " -p ", pfb_file ,
                        " -s ", penn_path,
                        " -o ", run_dir)
  return(list(cnv=command_cnv, flank=command_flank))
}
