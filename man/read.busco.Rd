% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.busco.R
\name{read.busco}
\alias{read.busco}
\title{Read the *_full_table.tsv output from BUSCO}
\usage{
read.busco(path)
}
\arguments{
\item{path}{path to BUSCO *_full_table.tsv output (or whatever the user may have renamed it)}
}
\description{
Read the *_full_table.tsv output from BUSCO
}
\details{
A simple wrapper for read.table() that will read in the *_full_table.tsv file output by BUSCO
and give sensible column names
}
