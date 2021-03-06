% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runExtractInheritance.R
\name{runExtractInheritance}
\alias{runExtractInheritance}
\title{Annotate original PennCNV trio data table with inheritance status using trio state decomposition}
\usage{
runExtractInheritance(filename_orig, filename_merged)
}
\arguments{
\item{filename_orig}{A path to output file of PennCNV trio}

\item{filename_merged}{A path to a result of merging with PennCNV of the filename_orig}
}
\value{
Returns paths to and writes files with status of inheritance, one of "denovo"/"inherited"/"ambiguous" for offspring and "call_in_parent" for parents (with and without parental data),
as well as file with no inheritance information for further merging. In case if no merging is planned, a log file with original PennCNV trio states is also written for further in the pipeline.
}
\description{
Annotate original PennCNV trio data table with inheritance status using trio state decomposition
}
