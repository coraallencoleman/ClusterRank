
# Clustered Rankings  
An R package for creating complete rankings for a list of items with mixture model clustering & visualizations   
 
## Functions  
**ClusterRankBin**
  Assigns ranks then clusters to each item in a list based on Binomial data. Calls npmleBin()

  Args:
    y: number of binomial distributed events
    n: number of attempts
    k: number of starting clusters desired. Defaults to length(y)
    scale: scale on which to do ranking
    weighted: boolean indicating if inverse variance weighted is used
    n.iter: iterations used in EM algorithm
    n.samp: number of samples from posterior distribution
    row_names: optional row names argument

  Returns:
      list including ranked_table, posterior, theta, pr_theta


**ClusterRankPois**
  Assigns ranks then clusters to each item in a list based on Poisson data. Calls npmlePois()

  Args:
    y: number of poisson distributed events
    ti: time vector. Length of time.
    k: number of starting clusters desired. Defaults to length(y)
    scale: scale for ranking
    weighted: boolean indicating if inverse variance weighted is used
    n.iter: iterations used in EM algorithm
    n.samp: number of samples from posterior distribution
    row_names: optional row names argument

  Returns:
      list including ranked_table, posterior, theta, pr_theta

**ClusterRankNorm**
  Assigns ranks then clusters to each item in a list based on Normal data. Calls npmleNorm()

  Args:
    y: means for each item to be ranked
    se: standard errors for each mean y
    k: number of starting clusters desired. Defaults to length(y)
    scale: scale for ranking
    weighted: boolean indicating if inverse variance weighted is used
    n.iter: iterations used in EM algorithm
    n.samp: number of samples from posterior distribution
    row_names: optional row names argument

  Returns:
      list including ranked_table, posterior, theta, pr_theta

**PlotClusterRank(ClusterRank, xlab=NULL, maintitle=NULL)**
creates a visualization using the result of ClusterRankBin, ClusterRankPois or ClusterRankNorm. Shows ranks with clusters and confidence intervals of ranks. 

## Data Accepted 
See data folder for examples.  

**Binomial Data**
requires: y count data, n trials
optional: row names for items

**Poisson Data**
requires: y counts, t time
optional: row names for items

**Normal Data**
requires: means, standard deviations
optional: row names for items

## Testing
To test using test_that package, run: 
devtools::test()

