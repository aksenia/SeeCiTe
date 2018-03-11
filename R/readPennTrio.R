readPennTrio <- function(filename) {
  res <- read.table(filename,
                    stringsAsFactors = F,
                    col.names = c("coordcnv", "numsnp", "length", "state",
                                  "samplefile", "startsnp", "endsnp", "relation", "triostate")) %>%
    separate(state, c("HMMstate", "copynumber"), sep = ",") %>%
    mutate(sample = gsub(".*\\.", "", samplefile),
           numsnp = as.integer(gsub(".*=", "", numsnp)),
           length = as.integer(gsub(".*=|,", "", length)),
           copynumber = as.integer(gsub(".*=", "", copynumber)),
           startsnp = gsub(".*=", "", startsnp),
           endsnp = gsub(".*=", "", endsnp),
           triostate = gsub(".*=", "", triostate)) %>%
    separate(coordcnv, c("chr", "coord"), sep = ":", remove = F) %>%
    separate(coord, c("start", "end"), sep = "-", remove = T) %>%
    mutate(start = as.integer(start),
           end = as.integer(end),
           length = as.integer(length))
  res
}
