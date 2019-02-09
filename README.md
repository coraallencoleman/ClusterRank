
# Clustered Rankings  
An R package for creating complete rankings for a list of items with mixture model clustering & visualizations   
 
## Functions  
**ClusterRank(y,n=NULL,se=NULL,ti=rep(1,length(y)),k=NULL,datatype, scale=identity,weighted=TRUE,n.iter=1000,n.samp=10000,row_names=NULL)**  
creates a full ranking from raw normal, binomial, or Poisson data  
cassigns clusters to items within the ranking using a mixture model   

**PlotClusterRank(ClusterRank, xlab=NULL, maintitle=NULL)**     
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

# TODO break into small

## Distribution-Specific:  
1. filling EM algorithm
2. generating posterior samples
3. generating CIs
4. final dataframe
5. posterior dist of theta for each county. comes from Ez (from EM)

## Common:  
sorting posterior samples (get order stat)
assigning clusters to order stats
assigning ranks

# TODO optional output
posterior distribution of theta for each order statistic. What is the distribution of theta for the min item, max item. We can use this to illustrate compromise how clusters get assigned to rank positions. (might be useful to save/return optionally) take existing sorted sample thing and use tapply to turn into empirrical proportion count.
df
order stats along rows 1:N
posterior samples along columns
table this ^ and divide by nsamp
conver to 1:N (rows) col: 1:theta

we want a table for each of the order stats:
ntheta_1/nsamp

this will be useful for presenting and paper. It will be a more clear picture of which cluster you belong to, but assignmting counties to order stats is harder.

# Deadline  
 - prototype package for March  
 - 15 slide talk  for March 25  
