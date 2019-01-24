
# Clustered Rankings  
An R package for creating complete rankings for a list of items with mixture model clustering & visualizations   
 
## Functions  
**ClusterRank(y,n=NULL,se=NULL,ti=rep(1,length(y)),k=NULL,datatype, scale=identity,weighted=TRUE,n.iter=1000,n.samp=10000,row_names=NULL)**  
creates a full ranking from raw normal, binomial, or Poisson data  
cassigns clusters to items within the ranking using a mixture model   

**PlotClusterRank(ClusteredRanking, xlab=NULL, maintitle=NULL)**     
creates a visualization using the result of ClusterRank.
Shows ranks with clusters and confidence intervals of ranks. 

## Data Accepted 
See data folder for examples.  
### Binomial Data
requires: y count data, n trials
optional: row names for items

first column should be item names; second should be successes; third column should be number of trials.  

### Poisson Data
requires: y counts, t time
optional: row names for items

### Normal Data
requires: means, standard deviations
optional: row names for items

## Hidden Functions TODO
createClusters
ssignRanksClusters
cleanResults
getmode
plotClusterRanks


