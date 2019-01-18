# Rank and Cluster Pisson CT Example Data 

library(tidyverse)
library(coda)
library(reshape2)
library(clue)
library(Hmisc)
library(RColorBrewer)

setwd("/Users/cora/git_repos/RankingMethods")

normData <- read_csv("data/normal_ct.csv")

normData <- normData %>% filter(!is.na(mean)) #removes missing

mean_rc <- rank_cluster.norm(normData$mean,normData$sd,n.iter=1000,row_names=normData$county)
mean_rc2 <- rank_cluster.norm(normData$mean,normData$sd,row_names=normData$county,weighted=FALSE)
mean_x10_rc <- rank_cluster.norm(normData$mean*10,normData$sd*10,row_names=normData$county)
mean_x100_rc <- rank_cluster.norm(normData$mean*100,normData$sd*100,row_names=normData$county)

mean_rc_rnk <- rank_cluster.norm(normData$mean,normData$sd,row_names=normData$county,scale=rank)
mean_rc_rnk2 <- rank_cluster.norm(normData$mean,normData$sd,row_names=normData$county,scale=rank,weighted=FALSE)
mean_x10_rc <- rank_cluster.norm(normData$mean*10,normData$sd*10,row_names=normData$county)

plot_rt(rc = mean_rc_rnk)

#TODO
# teen_wi <- read_csv("mean_teen_wi.csv")
# teen_wi <- teen_wi %>% filter(!is.na(teen_sd))
# 
# teen_rc <- rank_cluster.norm(teen_wi$teen_sd,teen_wi$teen_pop,row_names=teen_wi$county)
# teen_x10_rc <- rank_cluster.norm(teen_wi$teen_sd*10,teen_wi$teen_pop*10,row_names=teen_wi$county)
