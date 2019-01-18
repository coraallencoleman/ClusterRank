
# Clustered Rankings  
An R package for creating complete rankings for a list of items with mixture model clustering & visualizations   
 
## Functions  
**ClusterRank(dataframe, datatype, weighted, scale, n.iter, row_names)**  
creates a full ranking from raw normal, binomial, or Poisson data  
cassigns clusters to items within the ranking using a mixture model   

**PlotClusterRank(result of cluster rank, xlab, maintitle)**     
creates a visualization showing the ranks with clusters and confidence intervals of ranks  

## Data Format  
See data folder for examples.  
### Binomial Data
first column should be item names; second should be successes; third column should be number of trials.  

### Poisson Data
first column should be item names; second should be successes; third column should be number of trials.  

### Normal Data
first column should be item names; second should be means; third column should be number of standard deviations.  

