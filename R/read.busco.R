#' Read the *_full_table.tsv output from BUSCO
#'
#' @param path path to BUSCO *_full_table.tsv output (or whatever the user may have renamed it)
#'
#' @details
#' A simple wrapper for read.table() that will read in the *_full_table.tsv file output by BUSCO
#' and give sensible column names
#'

#' @export
read.busco<-function(path) {
  dat<-read.table(file = path, sep = '\t', quote = '', header = F, fill = T)
  names(dat)<-c('busco_id', 'status', 'sequence', 'gene_start',
                'gene_end', 'strand', 'score', 'length', 'orthodb_url', 'description')
  return(dat)
}
