library(parallel)
library(rtracklayer)
library(ggplots2)

#
# meth_gr : methylation GRanges object
# it needs to have a numeric mcol nammed "value"
# that represents the methylation level at that position

# feature_gr : GRange object of the features of interest
# it needs to have a character (or factor) mcol nammed "feature"
# that represents the feature unique identifier



# Our WGBS is fucked up and the start and end are inverted in the bedGraph file
# Load methylation data
meth_gr = import("/Volumes/abacus/spapillo/K27M_peaks/wgbs/human/BT245-C19/BT245-C19_cells_BS_1.profile.bedGraph")
meth_gr$value = as.numeric(meth_gr$score)
start(meth_gr) = end(meth_gr)
end(meth_gr) = start(meth_gr) + 1

# Below I'm manipulating the NGSPLOT database to use it directly
# for re-ordering later. But you can use any properly formatted GRanges object
# Load annotation
load("~/bin/ngsplot/database/hg19/hg19.ensembl.genebody.protein_coding.RData")
load("~/Desktop/BT245.RData")

# This is formatting specific to my stuff
genome.coord = genome.coord[paste(genome.coord$gname, genome.coord$tid, sep=":") %in% gene_names, ]
ensemble_gr = GRanges(seqnames=genome.coord$chr, strand=genome.coord$strand,
        ranges = IRanges(start = genome.coord$start, end = genome.coord$end))
ensemble_gr$gname = genome.coord$gname
ensemble_gr$tid = genome.coord$tid
ensemble_gr$feature = paste(ensemble_gr$gname, ensemble_gr$tid, sep=":")

# Remove unecessary things..
gr_list = prune_dataset(meth_gr,ensemble_gr,upstream=10000,downstream=10000)

# Remove the large objects that aren't used anymore
rm(meth_gr,ensemble_gr,genome.coord)


gr_list$feature_gr_list = lapply(gr_list$feature_gr_list, function(i) {
  i$feature = paste(i$gname, i$tid, sep=":")
  i })
# Run the thing.
tmp_gr_list = gr_list
tmp_gr_list[[1]] = gr_list[[1]][c('chr21','chr22')]
tmp_gr_list[[2]] = gr_list[[2]][c('chr21','chr22')]


feature_methylation = process_data(tmp_gr_list$meth_gr_list, tmp_gr_list$feature_gr_list, upstream = 10000, downstream = 10000, feature_perc = 0.01, bin_size = 100, mc.cores=4)

width(tmp) = 1
  data_list = test(chr1_meth,head(ensemble_gr,n=100))
  df = Reduce(rbind.data.frame, data_list[-2])


