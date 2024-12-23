# CSBfilter
`CSBfilter` is a set of functions written in R to generate a list of sequences that should be kept based on gene content completeness and sequence lengths. Gene content completeness is assessed using [BUSCO](https://busco.ezlab.org/) or the recently released [compleasm](https://github.com/huangnengCSU/compleasm). Artifactual duplicates in genome assemblies are somewhat of an inevitability, and many assemblies pipelines include a duplicate contig purging step. Even still, purging is often imperfect, and even the highest quality genome assemblies still have some number of duplicates. The functions I present here, are intended to add a "finishing touch", to limit the number of BUSCO duplicates (D) in a genome assembly.

## When should CSBfilter be used?
The functions in `CSBfilter` would probably be best if used just prior to scaffolding and just after. The function `busco.contig.filter()` is written with the intention that the user will run it after purging their assembly of haplotigs.  The function `scaffold.filter()` is written specifically for assemblies that have been curated using the [Juicebox ](https://github.com/aidenlab/Juicebox) HiC analysis/visualization suite.

### `busco.contig.filter`
This function takes two main arguments and an optional third. 
1. `contig.info` which is a two column dataframe of sequence names and their respective lengths
2. `busco.out` or `compleasm.out` which are the tsv outputs from either BUSCO or compleasm. Please supply just one or the other.
3. `out.prefix` (optional; defaults to NULL) an output prefix for the `.txt` format list of the sequence names that should be *KEPT*. Not supplying one will print the list to console without saving a file.

Note: if the user wants to output the inverse list (that is a list of sequences that should be *REMOVED*), then the user should do something like the following:

```
out.list<-busco.contig.filter(contig.info = contig.info, busco.out = busco.out)
remove.list<-out.list[!(out.list %in% contig.info[,1])]
```
Where `contig.info` contains sequences names in the 1st column.

### `scaffold.filter`
This function takes two main arguments and a third optional one.
1. `jbat.list` is a list object created from the helper function `read.jbat.review`
2. `busco.dat` the full_tsv output from BUSCO. Unfortunately, I wrote this prior to learning about compleasm, so this only works with BUSCO full_table.tsv. compleasm does output a BUSCO formatted list, so you can just use that for now.
3. `out` (optional, defaults to NULL): an output prefix for the lists that will be saved to disk.

This function will write 3 files to disk: a list of sequences to keep, a list of sequences to remove, and a summary of the number of BUSCO genes that are on each sequence.

### How the *.filter functions work
Both functions follow a simple logic: if a sequence has a BUSCO single copy gene or a fragmented gene on it, keep it. If a sequence doesn't have ANY BUSCO genes on it, keep it. If a sequence contains duplicated BUSCO genes, keep the sequence that is longest. `scaffold.filter` has an additional layer to it, which is that it will remove sequences labelled 'debris' during the manual curation process of Juicebox. If a 'debris' sequence contains a single or fragmented BUSCO gene on it, it will be a put in the 'keep' list.

### Helper functions
I've also written wrapper functions that will help the user easily read in `*_full_table.tsv` from BUSCO and the `*.review.assembly` from Juicebox after manual curation. The `read.busco()` function is somewhat trivial, though please note that I've written both `*.filter` functions to use the columns I named in `read.busco`.

#### Additional notes
If you have a genome assembly that is not scaffolded using HiC data, then you must use `busco.contig.filter` since you won't have a `*.review.assembly` file.
