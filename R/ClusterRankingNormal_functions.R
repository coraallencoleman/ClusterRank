#Cluster Ranking functions for Normal Data
#note: replaced theta with mu

library(tidyverse)
library(coda)
library(reshape2)
library(clue)
library(Hmisc)
library(RColorBrewer)

#key point with normal cversion:
#comes in as y, se
#unkonwn mu, p_mu. We assume se is set here.
#mu i hat =
npmle.norm <- function(mean, se, n, k=NULL,n.iter=1000,row_names=NULL) {
  if (is.null(k)) {
    mu<-sort(y) #sorted means. TODO why do we sort here?
    sd <- sd
    #sd <-sd(mu) #TODO is this who we should calculate this? is this in the correct order?
    k<-length(mu) #k = number of units to rank
  } else {
    mu <- seq(min(y),max(y),length=k)
  }
  p_mu <- seq(min(y),max(y),length=k) #TODO how to do this for normal? is this to make evenly spaced means for groupings?

  E_z <- matrix(NA,length(y),k) #TODO use y for incoming data
  for (j in 1:n.iter) {
    for (i in 1:k) {
      #TODO check this
      E_z[,i] <- log(p_mu[i])+dnorm(y, mu[i], se, log=TRUE) #uses the se that comes in with data
    }
    E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x)))))
    p_mu <- apply(E_z,2,y)
    mu <- ((y/se^2)%*%E_z)/((1/se^2)%*%E_z)
  }

  ord<-order(mu)
  mu<-c(mu[ord])
  p_mu<-p_mu[ord]

  p_mu <- tapply(p_mu,cumsum(!duplicated(round(mu,8))),sum)
  mu <- mu[!duplicated(round(mu,8))]

  E_z <- matrix(NA,length(y),length(mu))
  for (i in 1:length(mu)) {
    E_z[,i] <- log(p_mu[i])+dnorm(y/n,mu[i],log=TRUE)
  }
  E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x)))))

  rownames(E_z)<-row_names
  colnames(E_z)<-signif(mu,3)

  return(list(mu=mu, p_mu=p_mu, post_mu=E_z))
}

rank_cluster.norm <- function(mean,sd,k=NULL,scale=identity,weighted=TRUE,n.iter=1000,n.samp=10000,row_names=NULL) {
  N <- length(mean)

  npmle_res <- npmle.norm(mean,sd,k,n.iter,row_names)

  smp <- apply(npmle_res$post_mean,1,
               function(x,mean,n.samp)
                 sample(mean,n.samp,replace=TRUE,prob=x), #TODO update this
               mu=scale(npmle_res$mean),n.samp=n.samp)
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
  grp <- match(apply(smp.ord,1,getmode),scale(npmle_res$mu))[rnk]
  grp <- factor(grp)
  p_grp <- npmle_res$post_mu[cbind(1:N,as.numeric(grp))]
  levels(grp) <- signif(npmle_res$mu,3)

  ord <- order(rnk)

  CI <- normconf(y,n)

  ranked_table <- data_frame(name=row_names,rank=rnk,group=factor(grp),
                             y=y,n=n,p=y/n,
                             p_LCL=CI[,2],p_UCL=CI[,3],
                             pm=c(npmle_res$post_mu%*%npmle_res$mu),
                             p_grp=p_grp)
  ranked_table <- ranked_table[ord,]
  ranked_table$name <- factor(ranked_table$name,levels=ranked_table$name,ordered=TRUE)

  posterior <- npmle_res$post_mu[ord,]

  return(list(ranked_table=ranked_table,posterior=posterior,mu=npmle_res$mu,pr_mu=npmle_res$p_mu))
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

plot_rt <- function(rc,xlab="Proportion") {
  #rc is the result of rank_cluster.norm()
  post_df <- melt(rc$posterior)
  post_df$group <- rc$ranked_table$group[match(post_df$Var1,rc$ranked_table$name)]
  post_df$p_grp <- rc$ranked_table$p_grp[match(post_df$Var1,rc$ranked_table$name)]

  return(ggplot(rc$ranked_table,aes(y=name,x=p,color=group,alpha=p_grp))+
           geom_point(pch=3)+
           geom_point(aes(x=pm),pch=4)+
           geom_point(data=post_df,aes(y=Var1,x=as.numeric(Var2),color=group,size=value,alpha=value))+
           geom_errorbarh(aes(xmin=p_LCL,xmax=p_UCL),height=0)+
           scale_y_discrete("",limits=rev(levels(rc$ranked_table$name)))+
           scale_x_continuous(xlab,breaks=rc$mu[!duplicated(round(rc$mu,2))],
                              labels=round(rc$mu[!duplicated(round(rc$mu,2))],3),minor_breaks=rc$mu)+
           scale_color_manual(values=rep(brewer.pal(8,"Dark2"),1+floor(length(levels(rc$ranked_table$group))/8)))+
           scale_size_area(max_size=5)+
           scale_alpha(limits=c(0,1),range=c(0,1))+
           theme_bw()+
           guides(color=FALSE,size=FALSE,alpha=FALSE))
}

