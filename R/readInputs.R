#' Load and combine multiple input files
#'
#' @param args A list holding all necessary data files. eg  "triocnv_file" (input formatted penn triocnv file for offspring),
#' "probecoord_file" (file created by _runcommands script in flanks ),
#' "snp_flank_file" (_candidateCnv.txt file in flanks),
#' "snp_cnv_log_file" (_candidateCnvLOG.txt file in cnv locus),
#' "cnv_qcsum_file" (*.qcsum file from PennCNV), "dataset" (dataset string used throughout the analysis),
#' "cache_id" (full path to where to store rds cache files for computation havier steps),
#' "merge_trace" (table in a custom format with original PennCNV trio status and merging statistics, use runExtractInheritance to create this table)
#'
#' @return A list with "data": a data table with SNP level values of LRR and BAF for each individual in a trio;
#' if "cnv_qcsum_file" in args is supplied, a "qcsum" table with the genome-wide LRR and BAF summaries for each individual;
#' if "merge_trace" in args is supplied, a "merge" table with original PennCNV trio states and merging segment counts
#' @export
#'
#'
readInputs <- function(args){
  penn_cn_state <- decodeHMMstates()

    print(paste0("Reading merged formatted PennCNV trio file for offspring ", args[["triocnv_file"]]))
    penncnv <- readPennTrioFormat(args[["triocnv_file"]]) %>%
      dplyr::left_join(penn_cn_state, by = c("state")) %>%
      dplyr::select(-dplyr::matches("snp"), -relation)  %>%
      tidyr::unite(cnvID, sample, coordcnv, sep="_", remove=F) %>%
      dplyr::mutate(CNVID=paste(args[["dataset"]], cnvID, sep="_"))

    print(paste0("Reading probes file ", args[["probecoord_file"]]))
    probescoord <- read.table(args[["probecoord_file"]],
                            stringsAsFactors = F,
                            col.names = c("Name",	"Chr", "Position","PFB"))

    ## some arrays, e.g. Affy6.0 have duplicate probes
    print("Collapsing duplicate probes")
    probescollapse <- probescoord %>%
      dplyr::group_by(Chr, Position, PFB) %>%
      dplyr::summarise(Name=Name[[1]])

    #### find out which probes fall into predicted CNV segments
    penncnv_dict <- penncnv %>%
      tidyr::separate(coordcnv, c("Chr", "coord"), sep=":", remove=F) %>%
      tidyr::separate(coord, c("start", "end"), sep="-", remove = T) %>%
      dplyr::select(CNVID, sample, coordcnv, start, end, Chr) %>%
      tidyr::gather(key="pos", value="Position", -sample, -coordcnv, -Chr, -CNVID) %>%
      dplyr::mutate(Chr=as.integer(gsub("[a-z]*","",Chr)),
                  Position=as.integer(Position))

      lpdv <- split(penncnv_dict, penncnv_dict$CNVID)

    print("Annotating probes in CNV loci")
    suppressMessages(cnvProbes <- do.call("rbind", lapply(names(lpdv), function(n) {
      s <- lpdv[[n]]
      st <- as.numeric(gsub(".*:|-.*","", n))
      sn <- as.numeric(gsub(".*-","", n))
      schr <- as.numeric(gsub("[a-z]*","",gsub(".*_|:.*","", n)))
      s %>%
        dplyr::right_join(probescoord) %>%
        dplyr::mutate(locus="CNV",
                      CNVID=n) %>%
        dplyr::filter(Position >=st, Position<=sn, Chr==schr) %>%
        dplyr::select(CNVID, Name, locus) %>%
        as.data.frame()})))

    print(paste0("Reading extracted intensity for flanking regions ", args[["snp_flank_file"]]))
    suppressMessages(candCnvFlanks <- read.table(args[["snp_flank_file"]],
                              stringsAsFactors = F,
                              col.names = c("Name", paste("LRR", c("F", "M", "O"), sep = "_"),
                                            paste("BAF", c("F", "M", "O"), sep = "_"),
                                            paste("GENO", c("F", "M", "O"), sep = "_"), "Origin", "CNVID")) %>%
      dplyr::mutate(CNVID = gsub("_flanks_", "_", CNVID),
                    CNVID = gsub("\\.out$", "", CNVID)) %>%
      dplyr::left_join(cnvProbes, by = c("CNVID", "Name")) %>%
      dplyr::mutate(locus = ifelse(is.na(locus), "flank", locus)) %>%
      dplyr::filter(CNVID %in% penncnv$CNVID))

    print(table(penncnv$status))
    # We need assessment for the CNV locus only!
    print(paste0("Reading denovo test results for CNV regions ", args[["snp_cnv_log_file"]]))
    candCnvLog <- read.table(args[["snp_cnv_log_file"]],
                           stringsAsFactors = F,
                           sep = ",",
                           col.names = c("N", "CNVID")) %>%
      tidyr::separate(N, c("Note", "Message"), sep = "\\(de novo ") %>%
      dplyr::mutate(Note = gsub("NOTICE: ", "", Note),
                  Message = gsub(" \\)", "", Message),
                  Support = ifelse(grepl(":", Message), gsub("in trio.*:", "", Message),
                                   gsub("in trio.*$", "", Message))) %>%
      dplyr::select(-Message) %>%
      dplyr::filter(CNVID %in% penncnv$CNVID)

    suppressMessages(candidateCnvs <- candCnvFlanks %>%
      dplyr::left_join(candCnvLog, by = c("CNVID")) %>%
      dplyr::left_join(penncnv %>% dplyr::select(CNVID, status,type)) %>%
      tidyr::separate(CNVID, c("dataset", "sample", "coordcnv"), sep = "_") %>%
      #  separate(CNVID, c("dataset", "sample1", "sample2", "coordcnv"), sep = "_") %>%
      #  unite(sample, sample1, sample2, sep = "_") %>%
      dplyr::inner_join(probescollapse, by = c("Name")))


    suppressWarnings(candidateCnvsl <- candidateCnvs %>%
      dplyr::select(-Note, -PFB) %>%
      tidyr::gather(key, value, -Name, -Chr, -Position, -dataset,
                  -sample, -coordcnv, -locus, -Origin, -Support, -status, -type) %>%
      tidyr::separate(key, c("parameter", "relation")) %>%
      dplyr::mutate(value = as.numeric(value)) %>%
      dplyr::mutate(Support = gsub("= ", "=", Support),
                  Support = gsub("  ", " ", Support)) %>%
      tidyr::separate(Support, c("copynumber", "numsnp", "numsnp_F", "numsnp_M", "Pval"), sep = " ", remove = F) %>%
      dplyr::mutate(copynumber=as.integer(gsub(".*=", "", copynumber)),
                  numsnp=as.integer(gsub(".*=", "", numsnp)),
                  numsnp_F=as.integer(gsub(".*=", "", numsnp_F)),
                  numsnp_M=as.integer(gsub(".*=", "", numsnp_M)),
                  Pval=as.numeric(gsub(".*=", "", Pval))) %>%
      dplyr::group_by(dataset, sample, coordcnv, locus) %>%
      dplyr::mutate(numsnp = dplyr::n_distinct(Name)) %>%
      dplyr::ungroup())

  out <- list(data=candidateCnvsl)
  if (sum(grepl("cnv_qcsum_file", names(args)))>0) {
    print(paste0("Reading quality summary file ", args[["cnv_qcsum_file"]]))
    qcsum <- read.table(args[["cnv_qcsum_file"]],
                            stringsAsFactors = F,
                            header = T) %>%
      dplyr::mutate(sample=basename(File)) %>%
      dplyr::select(-File, -dplyr::matches("median"))
    out[["qcsum"]] <- qcsum
  }

  if (sum(grepl("merge_trace", names(args)))>0) {
    print(paste0("Reading merging log file ", args[["merge_trace"]]))
    merge <- read.table(args[["merge_trace"]],
                            stringsAsFactors = F,
                            col.names = c("coordcnv", "sample", "tstate", "mlog")) %>%
      dplyr::mutate(sample=gsub(":.*","",sample))
    out[["merge"]] <- merge
  }
  out
}


