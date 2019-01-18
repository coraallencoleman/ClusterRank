#Cluster Ranking Poisson functions
 #note: replaced theta with lambda

library(tidyverse)
library(reshape2)
library(clue)
library(Hmisc)
library(RColorBrewer)

npmle.pois <- function(y,ti=rep(1,length(y)),k=NULL,n.iter=1000,row_names=NULL) {
  #y, persontime ti
  if (is.null(k)) {
    lambda<-sort(y/ti) #sorted probabilities.
    k<-length(lambda) #number of groups to start
  } else {
    lambda <- seq(min(y/ti),max(y/ti),length=k)
  }
  p_lambda <- rep(1/k,k) #evenly spaced probabilities between 0 and 1 for groups?

  E_z <- matrix(NA,length(y),k)
  for (j in 1:n.iter) {
    for (i in 1:k) {
      E_z[,i] <- log(p_lambda[i])+dpois(y, ti*lambda[i],log=TRUE) #E-step
    }
    E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x))))) #E_z will be
    p_lambda <- apply(E_z,2,mean) #M-step
    lambda <- y%*%E_z/ti%*%E_z
  }

  ord<-order(lambda)
  lambda<-c(lambda[ord])
  p_lambda<-p_lambda[ord]

  p_lambda <- tapply(p_lambda,cumsum(!duplicated(round(lambda,8))),sum)
  lambda <- lambda[!duplicated(round(lambda,8))]

  E_z <- matrix(NA,length(y),length(lambda))
  for (i in 1:length(lambda)) {
    E_z[,i] <- log(p_lambda[i])+dpois(y,ti*lambda[i],log=TRUE)
  }
  E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x)))))

  rownames(E_z)<-row_names
  colnames(E_z)<-signif(lambda,3)

  return(list(lambda=lambda, p_lambda=p_lambda, post_lambda=E_z))
}

rank_cluster.pois <- function(y,n,k=NULL,scale=identity,weighted=TRUE,n.iter=1000,n.samp=10000,row_names=NULL) {
  #if we always label parameters as theta, p_theta, posterior should be same. shouldnt depend on data type
  N <- length(y)

  npmle_res <- npmle.pois(y,n,k,n.iter,row_names)

  smp <- apply(npmle_res$post_lambda,1,
               function(x,lambda,n.samp)
                 sample(lambda,n.samp,replace=TRUE,prob=x),
               lambda=scale(npmle_res$lambda),n.samp=n.samp)
  smp <- t(smp)
  smp.ord <- apply(smp,2,sort)

  if (weighted) wgt <- 1/pmax(.Machine$double.eps,apply(smp,1,var)) else wgt <- rep(1,N)

  loss <- matrix(NA,N,N)
  for (i in 1:N) {
    for (j in 1:N) {
      loss[i,j] <- wgt[i] * mean((smp[i,]-smp.ord[j,])^2)
    }
  }

  rnk <- as.numeric(solve_LSAP(loss))
  grp <- match(apply(smp.ord,1,getmode),scale(npmle_res$lambda))[rnk]
  grp <- factor(grp)
  p_grp <- npmle_res$post_lambda[cbind(1:N,as.numeric(grp))]
  levels(grp) <- signif(npmle_res$lambda,3)

  ord <- order(rnk)

  #break this into a new function
  CI <- poisconf(y,n) #this depends on data type. TODO find new for pois

  ranked_table <- data_frame(name=row_names,rank=rnk,group=factor(grp),
                             y=y,n=n,p=y/n,
                             p_LCL=CI[,2],p_UCL=CI[,3],
                             pm=c(npmle_res$post_lambda%*%npmle_res$lambda),
                             p_grp=p_grp)
  ranked_table <- ranked_table[ord,]
  ranked_table$name <- factor(ranked_table$name,levels=ranked_table$name,ordered=TRUE)

  posterior <- npmle_res$post_lambda[ord,]

  return(list(ranked_table=ranked_table,posterior=posterior,lambda=npmle_res$lambda,pr_lambda=npmle_res$p_lambda))
}

#same for binomial
# getmode <- function(v) { #same
#   uniqv <- unique(v)
#   uniqv[which.max(tabulate(match(v, uniqv)))]
# }
#
# plot_rt <- function(rc,xlab="Proportion") { #same for each data type
#   post_df <- melt(rc$posterior)
#   post_df$group <- rc$ranked_table$group[match(post_df$Var1,rc$ranked_table$name)]
#   post_df$p_grp <- rc$ranked_table$p_grp[match(post_df$Var1,rc$ranked_table$name)]
#
#   return(ggplot(rc$ranked_table,aes(y=name,x=p,color=group,alpha=p_grp))+
#            geom_point(pch=3)+
#            geom_point(aes(x=pm),pch=4)+
#            geom_point(data=post_df,aes(y=Var1,x=as.numeric(Var2),color=group,size=value,alpha=value))+
#            geom_errorbarh(aes(xmin=p_LCL,xmax=p_UCL),height=0)+
#            scale_y_discrete("",limits=rev(levels(rc$ranked_table$name)))+
#            scale_x_continuous(xlab,breaks=rc$lambda[!duplicated(round(rc$lambda,2))],
#                               labels=round(rc$lambda[!duplicated(round(rc$lambda,2))],3),minor_breaks=rc$lambda)+
#            scale_color_manual(values=rep(brewer.pal(8,"Dark2"),1+floor(length(levels(rc$ranked_table$group))/8)))+
#            scale_size_area(max_size=5)+
#            scale_alpha(limits=c(0,1),range=c(0,1))+
#            theme_bw()+
#            guides(color=FALSE,size=FALSE,alpha=FALSE))
# }

