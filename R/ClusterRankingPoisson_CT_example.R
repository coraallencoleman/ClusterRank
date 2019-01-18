# Rank and Cluster Pisson CT Example Data 

library(tidyverse)
library(coda)
library(reshape2)
library(clue)
library(Hmisc)
library(RColorBrewer)

setwd("/Users/cora/git_repos/RankingMethods")

poisData <- read_csv("data/lbw_ct.csv")

poisData <- poisData %>% mutate(nbw=births-lbw) %>% filter(!is.na(lbw))

lbw_rc <- rank_cluster.pois(poisData$lbw,poisData$births,row_names=poisData$county)
lbw_rc2 <- rank_cluster.pois(poisData$lbw,poisData$births,row_names=poisData$county,weighted=FALSE)
lbw_x10_rc <- rank_cluster.pois(poisData$lbw*10,poisData$births*10,row_names=poisData$county)
lbw_x100_rc <- rank_cluster.pois(poisData$lbw*100,poisData$births*100,row_names=poisData$county)

lbw_rc_rnk <- rank_cluster.pois(poisData$lbw,poisData$births,row_names=poisData$county,scale=rank)
lbw_rc_rnk2 <- rank_cluster.pois(poisData$lbw,poisData$births,row_names=poisData$county,scale=rank,weighted=FALSE)
lbw_x10_rc <- rank_cluster.pois(poisData$lbw*10,poisData$births*10,row_names=poisData$county)

plot_rt(rc = lbw_rc_rnk)

#TODO
# teen_wi <- read_csv("lbw_teen_wi.csv")
# teen_wi <- teen_wi %>% filter(!is.na(teen_births))
# 
# teen_rc <- rank_cluster.pois(teen_wi$teen_births,teen_wi$teen_pop,row_names=teen_wi$county)
# teen_x10_rc <- rank_cluster.pois(teen_wi$teen_births*10,teen_wi$teen_pop*10,row_names=teen_wi$county)
