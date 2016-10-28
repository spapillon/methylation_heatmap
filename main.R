# meth_gr : methylation GRanges object
# it needs to have a numeric mcol nammed "value"
# that represents the methylation level at that position

# feature_gr : GRange object of the features of interest
# it needs to have a character (or factor) mcol nammed "feature"
# that represents the feature unique identifier



# Our WGBS is fucked up and the start and end are inverted in the bedGraph file
# Load methylation data
library(rtracklayer)
source("functions.R")

meth_gr = import("/Volumes/vdu-032-aa/papillon/K36M/wgbs/methylation/bedgraph/C1_MPC_BS.bedgraph")
reference_gr = import("~/Documents/reference/mm10/mm10_refseq.bed")


# Format a lightweight object that is useful for faster processing
# (gets rid of meth_gr entries that are not within reference_gr) 
gr_list = prune_dataset(meth_gr,sample(reference_gr,1000),upstream=10000,downstream=10000)

system.time(feature_methylation <- process_data(gr_list$meth_gr_list, gr_list$feature_gr_list, upstream = 10000, downstream = 10000, feature_perc = 0.01, bin_size = 100, mc.cores=6))


