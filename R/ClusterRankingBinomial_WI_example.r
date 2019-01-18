# Rank and Cluster Binomial Data WI Example

library(tidyverse)
library(coda)
library(reshape2)
library(clue)
library(Hmisc)
library(RColorBrewer)

setwd("/Users/cora/git_repos/RankingMethods")

lbw_wi <- read_csv("lbw_wi.csv")
lbw_wi <- lbw_wi %>% mutate(nbw=births-lbw) %>% filter(!is.na(lbw))

lbw_rc <- rank_cluster.bin(lbw_wi$lbw,lbw_wi$births,row_names=lbw_wi$county)
lbw_rc2 <- rank_cluster.bin(lbw_wi$lbw,lbw_wi$births,row_names=lbw_wi$county,weighted=FALSE)
lbw_x10_rc <- rank_cluster.bin(lbw_wi$lbw*10,lbw_wi$births*10,row_names=lbw_wi$county)
lbw_x100_rc <- rank_cluster.bin(lbw_wi$lbw*100,lbw_wi$births*100,row_names=lbw_wi$county)

lbw_rc_rnk <- rank_cluster.bin(lbw_wi$lbw,lbw_wi$births,row_names=lbw_wi$county,scale=rank)
lbw_rc_rnk2 <- rank_cluster.bin(lbw_wi$lbw,lbw_wi$births,row_names=lbw_wi$county,scale=rank,weighted=FALSE)
lbw_x10_rc <- rank_cluster.bin(lbw_wi$lbw*10,lbw_wi$births*10,row_names=lbw_wi$county)

plot_rt(rc = lbw_rc_rnk)

#TODO
# teen_wi <- read_csv("lbw_teen_wi.csv")
# teen_wi <- teen_wi %>% filter(!is.na(teen_births))
# 
# teen_rc <- rank_cluster.bin(teen_wi$teen_births,teen_wi$teen_pop,row_names=teen_wi$county)
# teen_x10_rc <- rank_cluster.bin(teen_wi$teen_births*10,teen_wi$teen_pop*10,row_names=teen_wi$county)
