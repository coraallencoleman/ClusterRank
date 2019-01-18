
#Clustered Rankings
An R package for creating complete rankings with mixture model clustering & visualizations

##ClusterRank(dataframe, datatype, weighted, scale, n.iter, row_names)
creates a full ranking from raw normal, binomial, or Poisson data.
cassigns clusters to items within the ranking using a mixture model.

##PlotClusterRank(result of cluster rank, xlab, maintitle)
creates a visualization showing the ranks with clusters and confidence intervals of ranks