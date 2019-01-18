#functions for grouped ranking
#use this instead? devtools::use_package("dplyr") # Defaults to imports
library(ggplot2)
#library(coda) #what does coda do? (which command)
library(reshape2)
library(clue)
library(Hmisc)
library(RColorBrewer)
#TODO binomial data should be something else. use LBW for poisson
npmle.bin <- function(y,n,k=NULL,n.iter=1000,row_names=NULL) {
  #k is number of initial clusters
#  theta <- quantile(y/n,probs=(0:(k-1))/(k-1))
  if (is.null(k)) {
    theta<-sort(y/n) #sorted probabilities
    k<-length(theta) #k = number of units to rank
  } else {
    theta <- seq(min(y/n),max(y/n),length=k) #mass points of F We're estimating these
  }
  p_theta <- rep(1/k,k) #probabilities of each mass point. We're estimating these

  E_z <- matrix(NA,length(y),k) #expeted value of the probability that youre in each of hte groups
  #calculating the p that zij is equal to theta star j
  for (j in 1:n.iter) {
    for (i in 1:k) {
      #numerator
      E_z[,i] <- log(p_theta[i])+dbinom(y,n,theta[i],log=TRUE) #TODO E_z is log(pr associated with group) + log(quantile)?
    }
    E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x))))) #normalizes (pseudo zs). This is the E part of EM alg
    p_theta <- apply(E_z,2,mean) #M step: means over the
    theta <- y%*%E_z/n%*%E_z #theta for each group
  }
  #this reduces down to needed number of groups k
  ord<-order(theta)
  theta<-c(theta[ord]) #sorts
  p_theta<-p_theta[ord] #sorts

  p_theta <- tapply(p_theta,cumsum(!duplicated(round(theta,8))),sum)
  #cumsum numbers groups is ascending order. sums the pthetas that goes with each group. See pictures
  theta <- theta[!duplicated(round(theta,8))]

  E_z <- matrix(NA,length(y),length(theta))
  #final posterior probabilties for each county. Pr(county in group i)
  for (i in 1:length(theta)) {
    E_z[,i] <- log(p_theta[i])+dbinom(y,n,theta[i],log=TRUE)
  }
  E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x))))) #normalizes probabilities. subtract max to avoid underflow

  rownames(E_z)<-row_names
  colnames(E_z)<-signif(theta,3) #group names are rounded

  return(list(theta=theta, p_theta=p_theta, post_theta=E_z))
  #(priors, posterior)
}

rank_cluster.bin <- function(y,n,k=NULL,scale=identity,weighted=TRUE,n.iter=1000,n.samp=10000,row_names=NULL) {
  #assigns ranks then clusters to each item in a list for binomial data
  N <- length(y)

  npmle_res <- npmle.bin(y,n,k,n.iter,row_names)

  smp <- apply(npmle_res$post_theta,1,
               function(x,theta,n.samp)
                 sample(theta,n.samp,replace=TRUE,prob=x),
               theta=scale(npmle_res$theta),n.samp=n.samp)
  smp <- t(smp)
  smp.ord <- apply(smp,2,sort)
  if (weighted) #inverse variance weighting
    wgt <- 1/pmax(.Machine$double.eps,apply(smp,1,var)) #if variance is zero, uses v small value to make it impossible to reassign to new group
  else wgt <- rep(1,N)

  loss <- matrix(NA,N,N)
  for (i in 1:N) {
    for (j in 1:N) {
      loss[i,j] <- wgt[i] * mean((smp[i,]-smp.ord[j,])^2)
    }
  }

  rnk <- as.numeric(clue::solve_LSAP(loss))
  grp <- match(apply(smp.ord,1,getmode),scale(npmle_res$theta))[rnk] #matches rank positions to groups using mode. The mode version minimizes indicator (see pic)
  #^ We could replace this with squared error diff. See pic.
  grp <- factor(grp)
  p_grp <- npmle_res$post_theta[cbind(1:N,as.numeric(grp))]
  levels(grp) <- signif(npmle_res$theta,3) #labels

  ord <- order(rnk)

  CI <- Hmisc::binconf(y,n) #creating confidence intervals

  ranked_table <- data_frame(name=row_names,rank=rnk,group=factor(grp),
                             y=y,n=n,p=y/n,
                             p_LCL=CI[,2],p_UCL=CI[,3],
                             pm=c(npmle_res$post_theta%*%npmle_res$theta),
                             p_grp=p_grp)
  ranked_table <- ranked_table[ord,]
  ranked_table$name <- factor(ranked_table$name,levels=ranked_table$name,ordered=TRUE)

  posterior <- npmle_res$post_theta[ord,]

  return(list(ranked_table=ranked_table,posterior=posterior,theta=npmle_res$theta,pr_theta=npmle_res$p_theta))
}

getmode <- function(v) {
  #retrieves mode from list v
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

plot_rt <- function(rc,xlab="Proportion") {
  post_df <- reshape2::melt(rc$posterior)
  post_df$group <- rc$ranked_table$group[match(post_df$Var1,rc$ranked_table$name)]
  post_df$p_grp <- rc$ranked_table$p_grp[match(post_df$Var1,rc$ranked_table$name)]

  return(ggplot2::ggplot(rc$ranked_table,aes(y=name,x=p,color=group,alpha=p_grp))+
    ggplot2::geom_point(pch=3)+
    ggplot2::geom_point(aes(x=pm),pch=4)+
    ggplot2::geom_point(data=post_df,aes(y=Var1,x=as.numeric(Var2),color=group,size=value,alpha=value))+
    ggplot2::geom_errorbarh(aes(xmin=p_LCL,xmax=p_UCL),height=0)+
    ggplot2::scale_y_discrete("",limits=rev(levels(rc$ranked_table$name)))+
    ggplot2::scale_x_continuous(xlab,breaks=rc$theta[!duplicated(round(rc$theta,2))],
                     labels=round(rc$theta[!duplicated(round(rc$theta,2))],3),minor_breaks=rc$theta)+
    ggplot2::scale_color_manual(values=rep(RColorBrewer::brewer.pal(8,"Dark2"),1+floor(length(levels(rc$ranked_table$group))/8)))+
    ggplot2::scale_size_area(max_size=5)+
    ggplot2::scale_alpha(limits=c(0,1),range=c(0,1))+
    ggplot2::theme_bw()+
    ggplot2::guides(color=FALSE,size=FALSE,alpha=FALSE))
}

