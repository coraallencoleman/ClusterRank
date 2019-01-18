# Rank and Cluster Binomial Data CT Example

library(tidyverse)
library(coda)
library(reshape2)
library(clue)
library(Hmisc)
library(RColorBrewer)

setwd("/Users/cora/git_repos/RankingMethods")

setwd("/Users/cora/git_repos/RankingMethods")

binData <- read_csv("data/lbw_ct.csv")

binData <- binData %>% mutate(nbw=births-lbw) %>% filter(!is.na(lbw))

lbw_rc <- rank_cluster.bin(binData$lbw,binData$births,row_names=binData$county)
lbw_rc2 <- rank_cluster.bin(binData$lbw,binData$births,row_names=binData$county,weighted=FALSE)
lbw_x10_rc <- rank_cluster.bin(binData$lbw*10,binData$births*10,row_names=binData$county)
lbw_x100_rc <- rank_cluster.bin(binData$lbw*100,binData$births*100,row_names=binData$county)

lbw_rc_rnk <- rank_cluster.bin(binData$lbw,binData$births,row_names=binData$county,scale=rank)
lbw_rc_rnk2 <- rank_cluster.bin(binData$lbw,binData$births,row_names=binData$county,scale=rank,weighted=FALSE)
lbw_x10_rc <- rank_cluster.bin(binData$lbw*10,binData$births*10,row_names=binData$county)

plot_rt(rc = lbw_rc_rnk)

#TODO
# teen_wi <- read_csv("lbw_teen_wi.csv")
# teen_wi <- teen_wi %>% filter(!is.na(teen_births))
# 
# teen_rc <- rank_cluster.bin(teen_wi$teen_births,teen_wi$teen_pop,row_names=teen_wi$county)
# teen_x10_rc <- rank_cluster.bin(teen_wi$teen_births*10,teen_wi$teen_pop*10,row_names=teen_wi$county)
