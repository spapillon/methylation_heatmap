
###############################
#
###############################


# Run this first to save a lot of time and memory
prune_dataset = function(meth_gr, feature_gr, upstream = 1000, downstream = 1000) {
  # Drop datapoints outside of regions of interest
  upstream_gr = flank(feature_gr, width=upstream, start=T, ignore.strand=F)
  downstream_gr = flank(feature_gr, width=downstream, start=F, ignore.strand=F)
  overlaps = findOverlaps(meth_gr, reduce(reduce(c(upstream_gr, feature_gr, downstream_gr))))
  meth_gr = meth_gr[queryHits(overlaps)]
  # Split by chromosomes for faster processing
  # and parallelization
  meth_gr_list = split(meth_gr, seqnames(meth_gr))
  feature_gr_list = split(feature_gr, seqnames(feature_gr))
  meth_gr_list = meth_gr_list[intersect(names(meth_gr_list), names(feature_gr_list))]
  #feature_gr_list = feature_gr_list[intersect(names(meth_gr_list), names(feature_gr_list))]
  return(list(meth_gr_list = meth_gr_list, feature_gr_list = feature_gr_list))
}

# Removes unecessary entries in meth_gr and feature_gr 
# splits GRanges objects in chromosomes and launches
# main method on chromosomes individually
process_data = function(meth_gr_list, feature_gr_list, upstream = 1000, downstream = 1000,
                feature_perc = 0.01, bin_size = 100, mc.cores = 4 ) {

  # Fetch methylation values and adjust the feature column
  # to fit whatever there is in the feature mcol of feature_gr
  # (useful for post-processing)
  list_of_data_per_chr <- mclapply(sample(names(feature_gr_list)), function(chr) {
    message(paste("Processing",chr))
    data_list = process(meth_gr_list[[chr]], feature_gr_list[[chr]], upstream, downstream,
      feature_perc, bin_size)
    chr_data_frame = Reduce(rbind.data.frame, data_list)
    chr_data_frame$feature = feature_gr_list[[chr]][chr_data_frame$feature]$feature
    chr_data_frame
  }, mc.cores = mc.cores)
  # Crunch everything back together
  Reduce(rbind.data.frame, list_of_data_per_chr)
}


# Fetch methylation values as a function of distance for upstream & downstream
# and as a function of % for features. Methylation values are binned.
methylation_per_feature = function(meth_gr, feature_gr, upstream = 1000, downstream = 1000,
		feature_perc = 0.01, bin_size = 100) {
  upstream_gr = flank(feature_gr, width=upstream, start=T, ignore.strand=F)
  downstream_gr = flank(feature_gr, width=downstream, start=F, ignore.strand=F)

  # TSS plots are generated simply by setting width(feature_gr) = 1.
  # Thus, we can skip the within-feature fetching of methylation.
  if(any(width(feature_gr) != 1)) {
  feature_overlaps = findOverlaps(feature_gr, meth_gr)
  feature_meth = meth_by_percent(meth_gr, feature_gr, feature_overlaps)
  feature_meth_scaled = bin_methylation(feature_meth, 1, feature_perc)
  } else {
    feature_meth_scaled = NA
  }
  upstream_overlaps = findOverlaps(upstream_gr, meth_gr)
  upstream_meth = meth_by_distance(meth_gr, feature_gr, upstream_overlaps)
  upstream_meth_scaled = bin_methylation(upstream_meth, upstream, bin_size) 
  # Since distance is always positive, this will order everything correctly
  upstream_meth_scaled$bin = -upstream_meth_scaled$bin

  downstream_overlaps = findOverlaps(downstream_gr, meth_gr)
  downstream_meth = meth_by_distance(meth_gr, feature_gr, downstream_overlaps) 
  downstream_meth_scaled = bin_methylation(downstream_meth, upstream, bin_size) 

  list(upstream = upstream_meth_scaled, feature = feature_meth_scaled, downstream = downstream_meth_scaled)
}

# Averages methylation levels in bins
bin_methylation = function(feature_list, max, by) {
  Reduce(rbind.data.frame,lapply(1:length(feature_list), function(i) {
      bins = cut(feature_list[[i]][,'distance'], breaks=seq(0, max, by),include.lowest=T, labels=seq((0+by/2), (max-by/2), by))
      mean_meth_per_bin = tapply(feature_list[[i]][,'value'], bins, mean)
      data.frame(bin=as.numeric(names(mean_meth_per_bin)), meth= mean_meth_per_bin, feature=i)
  }))
}

# Loop through features and obtain every meth_gr entry that overlaps
# and compute the distance.
# This function represents majority of compute time.
meth_by_distance = function(meth_gr, feature_gr, overlaps) {
  sapply(1:length(feature_gr), function(i) {
    meth_i = subjectHits(overlaps[queryHits(overlaps) == i])
    distance =  distance(meth_gr[meth_i], feature_gr[i])
    cbind(distance=distance, value=meth_gr[meth_i]$value)
  })
}

# Uses meth_by_distance function with a temporary
# feature_gr object with width == 1. Doing so makes feature_gr
# end == start + 1. This allows to compute distance to feature 
# start and the divide by the feature witdh, thus obtaining a %.
# Since reverse-strand features start and actually the end(feature_gr)
# those entries should be 1 - x.
meth_by_percent = function(meth_gr, feature_gr, overlaps) {

  tmp_feature_gr = feature_gr
  width(tmp_feature_gr) = 1
  feature_meth = meth_by_distance(meth_gr, tmp_feature_gr, overlaps) 
  sapply(1:length(feature_gr), function(i) {
    percentage = feature_meth[[i]][,'distance'] / width(feature_gr[i])
    # 1 - x for reverse strand
    if(as.character(strand(feature_gr[i])) == "-") {
      percentage = 1 - percentage
    }
    feature_meth[[i]][,'distance'] = percentage
    feature_meth[[i]]
   })
}

plotit = function(df) {

  p = ggplot(data=df,aes(x=ordered(bin),y=ordered(feature), fill=meth))
  p = p + geom_tile()
  p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p = p + scale_fill_gradientn(colours=c("blue","white","red"), na.value = "lightgray")
  p

}

