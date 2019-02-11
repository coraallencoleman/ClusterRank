#functions for ClusterRank

library(ggplot2)
library(reshape2)
library(clue)
library(Hmisc)
library(RColorBrewer)

ClusterRankBin <- function(y,n=NULL,k=NULL,
                        scale=identity,weighted=TRUE,n.iter=1000,n.samp=10000,row_names=NULL) {
  # Assigns ranks then clusters to each item in a list based on Binomial data. Calls npmleBin()
  #
  # Args:
  #   y: number of binomial distributed events
  #   n: number of attempts
  #   k: number of starting clusters desired. Defaults to length(y)
  #   scale: scale on which to do ranking
  #   weighted: boolean indicating if inverse variance weighted is used
  #   n.iter: iterations used in EM algorithm
  #   n.samp: number of samples from posterior distribution
  #   row_names: optional row names argument
  #
  # Returns:
  #     list including ranked_table, posterior, cluster theta, pr_theta
  #
  N <- length(y)
  if(missing(n)){
      stop("n required for binomial data")
  }
  npmle_res <- npmleBin(y=y,n=n,k=k,n.iter=n.iter,row_names=row_names)
  # samples from posterior distribution. Samples from cluster thetas with prob x = post_theta
  smp <- apply(npmle_res$post_theta,1,
               function(x,theta,n.samp)
                 sample(theta,n.samp,replace=TRUE,prob=x),
               theta=scale(npmle_res$theta),n.samp=n.samp)
  smp <- t(smp) #transposes
  smp.ord <- apply(smp,2,sort) #sorts samples by ??
  if (weighted){ #inverse variance weighting
    wgt <- 1/pmax(.Machine$double.eps,apply(smp,1,var)) #if variance is zero, uses v small value to weight,
    #making it impossible to reassign a low variance estimate to wrong group
  }
  else {
    wgt <- rep(1,N)
  }
  #weighted square error loss rank optimization
  lossRank <- matrix(NA,N,N)
  for (i in 1:N) {
    for (j in 1:N) {
      lossRank[i,j] <- wgt[i] * mean((smp[i,]-smp.ord[j,])^2)
    }
  }
  totalRankLoss = sum(lossRank)
  rnk <- as.numeric(clue::solve_LSAP(lossRank))

  # square error loss cluster optimization
  lossCluster <- matrix(NA,N,length(npmle_res$theta))
  for (i in 1:N) {
    for (j in 1:length(npmle_res$theta)) {
      lossCluster[i,j] <- mean((smp[i,]-c(scale(npmle_res$theta)[j]))^2)
    }
  }
  totalClusterLoss <- sum(lossCluster)
  cluster <- apply(lossCluster, 1, which.min) #TODO check
  cluster <- factor(cluster)
  p_cluster <- npmle_res$post_theta[cbind(1:N,as.numeric(cluster))]
  levels(cluster) <- signif(npmle_res$theta,3) #labels

  ord <- order(rnk)
  CI <- matrix(ncol = 3, nrow = N)
  colnames(CI) = c("PointEst", "Lower", "Upper")
  rownames(CI) = c(rep("", times = N))
  CI <- Hmisc::binconf(y,n) #creating confidence intervals
  ranked_table <- data.frame(name=row_names,rank=rnk,cluster=factor(cluster),
                             y=y,n=n,est = CI[,1],
                             p_LCL=CI[,2],p_UCL=CI[,3],
                             posteriorMean=c(npmle_res$post_theta%*%npmle_res$theta),
                             p_cluster=p_cluster)
  ranked_table <- ranked_table[ord,]
  ranked_table$name <- factor(ranked_table$name,levels=ranked_table$name,ordered=TRUE)
  posterior <- npmle_res$post_theta[ord,]
  return(list(ranked_table=ranked_table,posterior=posterior,theta=npmle_res$theta, pr_theta=npmle_res$p_theta))
}

ClusterRankPois <- function(y,ti=rep(1,length(y)),k=NULL,
                        scale=identity,weighted=TRUE,n.iter=1000,n.samp=10000,row_names=NULL) {
  # Assigns ranks then clusters to each item in a list based on Poisson data. Calls npmlePois()
  #
  # Args:
  #   y: number of poisson distributed events
  #   ti: time vector. Length of time.
  #   k: number of starting clusters desired. Defaults to length(y)
  #   scale: scale for ranking
  #   weighted: boolean indicating if inverse variance weighted is used
  #   n.iter: iterations used in EM algorithm
  #   n.samp: number of samples from posterior distribution
  #   row_names: optional row names argument
  #
  # Returns:
  #     list including ranked_table, posterior, theta, pr_theta
  #
  N <- length(y)
  npmle_res <- npmlePois(y=y,ti=ti,k=k,n.iter=n.iter,row_names=row_names)
  smp <- apply(npmle_res$post_theta,1, #samples from a centered version of npmle_res$theta with pr = post_theta
               function(x,theta,n.samp)
                 sample(theta,n.samp,replace=TRUE,prob=x),
               theta=scale(npmle_res$theta),n.samp=n.samp) #why centered? to prevent underflow problems?
  smp <- t(smp)
  smp.ord <- apply(smp,2,sort)
  if (weighted) { #inverse variance weighting
    wgt <- 1/pmax(.Machine$double.eps,apply(smp,1,var)) #if variance is zero, uses v small value to weight,
    #making it impossible to reassign a low variance estimate to wrong group
  }
  else {
    wgt <- rep(1,N)
  }
  loss <- matrix(NA,N,N)
  for (i in 1:N) {
    for (j in 1:N) {
      loss[i,j] <- wgt[i] * mean((smp[i,]-smp.ord[j,])^2)
    }
  }
  totalRankLoss = sum(lossRank)
  rnk <- as.numeric(clue::solve_LSAP(loss))
  # square error loss cluster optimization
  lossCluster <- matrix(NA,N,length(npmle_res$theta))
  for (i in 1:N) {
    for (j in 1:length(npmle_res$theta)) {
      lossCluster[i,j] <- mean((smp[i,]-c(scale(npmle_res$theta)[j]))^2)
    }
  }
  totalClusterLoss <- sum(lossCluster)
  cluster <- apply(lossCluster, 1, which.min) #TODO check
  cluster <- factor(cluster)
  p_cluster <- npmle_res$post_theta[cbind(1:N,as.numeric(cluster))]
  levels(cluster) <- signif(npmle_res$theta,3) #labels

  ord <- order(rnk)
  CI <- matrix(ncol = 3, nrow = N)
  colnames(CI) = c("PointEst", "Lower", "Upper")
  rownames(CI) = c(rep("", times = N))
  ests <- WilsonHilfertyPoiCI(y, ti, conf.level=0.95)
  CI[, "PointEst"] <- ests[1]
  CI[, "Lower"] <- ests[2]
  CI[, "Upper"] <- ests[3]
  ranked_table <- data.frame(name=row_names,rank=rnk,group=factor(cluster),
                             y=y,ti=ti,est = CI[,1],
                             p_LCL=CI[,2],p_UCL=CI[,3],
                             posteriorMean=c(npmle_res$post_theta%*%npmle_res$theta),
                             p_cluster=p_cluster)
  ranked_table <- ranked_table[ord,]
  ranked_table$name <- factor(ranked_table$name,levels=ranked_table$name,ordered=TRUE)

  posterior <- npmle_res$post_theta[ord,]

  return(list(ranked_table=ranked_table,posterior=posterior,theta=npmle_res$theta,pr_theta=npmle_res$p_theta))
}

ClusterRankNorm <- function(y,n=NULL,se,k=NULL, scale=identity,
                            weighted=TRUE,n.iter=1000,n.samp=10000,row_names=NULL) {
  # Assigns ranks then clusters to each item in a list based on Normal data. Calls npmleNorm()
  #
  # Args:
  #   y: means for each item to be ranked
  #   se: standard errors for each mean y
  #   k: number of starting clusters desired. Defaults to length(y)
  #   scale: scale for ranking
  #   weighted: boolean indicating if inverse variance weighted is used
  #   n.iter: iterations used in EM algorithm
  #   n.samp: number of samples from posterior distribution
  #   row_names: optional row names argument
  #
  # Returns:
  #     list including ranked_table, posterior, theta, pr_theta
  #
  N <- length(y)
  if(c(missing(se))) {
    stop("se required for normal data")
  }
  npmle_res <- npmleNorm(y=y, se=se, k=k,n.iter=n.iter,row_names=row_names)

  smp <- apply(npmle_res$post_theta,1,
               function(x,theta,n.samp)
                 sample(theta,n.samp,replace=TRUE,prob=x),
               theta=scale(npmle_res$theta),n.samp=n.samp)
  smp <- t(smp)
  smp.ord <- apply(smp,2,sort)
  if (weighted) { #inverse variance weighting
    wgt <- 1/pmax(.Machine$double.eps,apply(smp,1,var)) #if variance is zero, uses v small value to weight,
                            #making it impossible to reassign a low variance estimate to wrong group
  }
  else {
    wgt <- rep(1,N)
  }
  loss <- matrix(NA,N,N)
  for (i in 1:N) {
    for (j in 1:N) {
      loss[i,j] <- wgt[i] * mean((smp[i,]-smp.ord[j,])^2)
    }
  }
  totalRankLoss = sum(lossRank)
  rnk <- as.numeric(clue::solve_LSAP(loss))

  # square error loss cluster optimization
  lossCluster <- matrix(NA,N,length(npmle_res$theta))
  for (i in 1:N) {
    for (j in 1:length(npmle_res$theta)) {
      lossCluster[i,j] <- mean((smp[i,]-c(scale(npmle_res$theta)[j]))^2)
    }
  }
  totalClusterLoss <- sum(lossCluster)
  cluster <- apply(lossCluster, 1, which.min) #TODO check
  cluster <- factor(cluster)
  p_cluster <- npmle_res$post_theta[cbind(1:N,as.numeric(cluster))]
  levels(cluster) <- signif(npmle_res$theta,3) #labels for clusters

  ord <- order(rnk)
  CI <- matrix(ncol = 3, nrow = N)
  colnames(CI) = c("PointEst", "Lower", "Upper")
  rownames(CI) = c(rep("", times = N))
  CI[, "PointEst"] <- y
  CI[, "Lower"] <- y - 1.96*se
  CI[, "Upper"] <- y + 1.96*se
  ranked_table <- data.frame(name=row_names,rank=rnk,group=factor(cluster),
                             y=y, se = se, est = CI[,1],
                             p_LCL=CI[,2],p_UCL=CI[,3],
                             posteriorMean=c(npmle_res$post_theta%*%npmle_res$theta),
                             p_cluster=p_cluster)
  ranked_table <- ranked_table[ord,]
  ranked_table$name <- factor(ranked_table$name,levels=ranked_table$name,ordered=TRUE)

  posterior <- npmle_res$post_theta[ord,]

  return(list(ranked_table=ranked_table,posterior=posterior,theta=npmle_res$theta, pr_theta=npmle_res$p_theta))
}

npmleBin <- function(y,n,k=NULL,n.iter=1000,row_names=NULL) {
  # Estimates clusters nonparametrically using an EM algorithm. Calculates the
  # probability each item will be assigned to each cluster.
  # Called by ClusterRankBin()
  #
  # Args:
  #   y: number of binomial distributed events
  #   n: number of attempts
  #   k: initial number of clusters. Defaults to length(y)
  #   n.iter: iterations used in EM algorithm
  #   row_names: optional row names argument
  #
  # Returns:
  #     list including clusters thetas, prior distribution for each p_theta,
  #                   posterior probabilities for each item's assignment to each cluster
  #
  if (is.null(k)) {
    theta<-sort(y/n) #sorted probabilities
    k<-length(theta) #k = number of units to rank
  } else {
    theta <- seq(min(y/n),max(y/n),length=k) #starting mass points of F. We're estimating these, along with p_theta
  }
  p_theta <- rep(1/k,k) #probabilities of each mass point.

  E_z <- matrix(NA,length(y),k) #expected value of the probability that you're in each of the k theta groups
  #calculating the p that z_{ij} is equal to theta star j
  for (j in 1:n.iter) {
    for (i in 1:k) {
      #numerator
      E_z[,i] <- log(p_theta[i])+dbinom(y,n,theta[i],log=TRUE)
    }
    E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x))))) #normalizes (pseudo zs). This is the E part of EM alg
    p_theta <- apply(E_z,2,mean) #M-step: means over the matrix
    theta <- y%*%E_z/n%*%E_z #calculates optimal theta for each group
  }
  #this reduces down to needed number of clusters (<=k)
  ord<-order(theta)
  theta<-c(theta[ord]) #sorts
  p_theta<-p_theta[ord] #sorts

  p_theta <- tapply(p_theta,cumsum(!duplicated(round(theta,8))),sum) #cumsum numbers clusters is ascending order. sums the pthetas that goes with each cluster. See pictures
  theta <- theta[!duplicated(round(theta,8))] #removes duplicate thetas

  E_z <- matrix(NA,length(y),length(theta))
  #final posterior probabilties for each county. Pr(county in cluster i)
  for (i in 1:length(theta)) {
    E_z[,i] <- log(p_theta[i])+dbinom(y,n,theta[i],log=TRUE)
  }
  E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x))))) #normalizes probabilities. subtracts max to avoid underflow

  rownames(E_z)<-row_names
  colnames(E_z)<-signif(theta,3) #cluster names are rounded

  return(list(theta=theta, p_theta=p_theta, post_theta=E_z))
}

npmlePois <- function(y,ti=rep(1,length(y)),k=NULL,n.iter=1000,row_names=NULL) {
  # Estimates clusters nonparametrically using an EM algorithm. Calculates the
  # probability each item will be assigned to each cluster.
  # Called by ClusterRankPois()
  #
  # Args:
  #   y: number of poisson distributed events
  #   ti: length of time. defaults to 1
  #   k: initial number of clusters. Defaults to length(y)
  #   n.iter: iterations used in EM algorithm
  #   row_names: optional row names argument
  #
  # Returns:
  #     list including prior for theta, prior distribution for each p_theta,
  #                   posterior probabilities for each item's assignment to each cluster
  #
  if (is.null(k)) {
    theta<-sort(y/ti) #sorted probabilities.
    k<-length(theta) #number of clusters to start
  } else {
    theta <- seq(min(y/ti),max(y/ti),length=k)
  }
  p_theta <- rep(1/k,k) #evenly spaced probabilities between 0 and 1 for clusters?

  E_z <- matrix(NA,length(y),k)
  for (j in 1:n.iter) {
    for (i in 1:k) {
      E_z[,i] <- log(p_theta[i])+dpois(y, ti*theta[i],log=TRUE) #E-step
    }
    E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x))))) #E_z will be
    p_theta <- apply(E_z,2,mean) #M-step
    theta <- y%*%E_z/ti%*%E_z
  }

  ord<-order(theta)
  theta<-c(theta[ord])
  p_theta<-p_theta[ord]

  p_theta <- tapply(p_theta,cumsum(!duplicated(round(theta,8))),sum)
  theta <- theta[!duplicated(round(theta,8))]

  E_z <- matrix(NA,length(y),length(theta))
  for (i in 1:length(theta)) {
    E_z[,i] <- log(p_theta[i])+dpois(y,ti*theta[i],log=TRUE)
  }
  E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x)))))
  rownames(E_z)<-row_names
  colnames(E_z)<-signif(theta,3)

  return(list(theta=theta, p_theta=p_theta, post_theta=E_z))
}

#key point with normal version:
#comes in as y, se, unknown theta, p_theta. We assume se is known here.
#theta i hat = see pics
npmleNorm <- function(y, se, k=NULL,n.iter=1000,row_names=NULL) {
  # Estimates clusters nonparametrically using an EM algorithm. Calculates the
  # probability each item will be assigned to each cluster.
  # Called by ClusterRankNorm()
  #
  # Args:
  #   y: mean for each item
  #   se: standard error for each mean y
  #   k: initial number of clusters. Defaults to length(y)
  #   n.iter: iterations used in EM algorithm
  #   row_names: optional row names argument
  #
  # Returns:
  #     list including prior for theta, prior distribution for each p_theta,
  #                   posterior probabilities for each item's assignment to each cluster
  #
  if (is.null(k)) {
    theta<-sort(y)
    k<-length(theta) #k = # clusters
  } else {
    theta <- seq(min(y),max(y),length=k)
  }
  p_theta <- rep(1/k, k)

  E_z <- matrix(NA,length(y),k)
  for (j in 1:n.iter) {
    for (i in 1:k) {
      E_z[,i] <- log(p_theta[i])+dnorm(y, theta[i], se, log=TRUE) #uses the se that comes in with data
    }
    E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x)))))
    p_theta <- apply(E_z,2,mean)
    theta <- ((y/se^2)%*%E_z)/((1/se^2)%*%E_z)
  }

  ord<-order(theta)
  theta<-c(theta[ord])
  p_theta<-p_theta[ord]

  p_theta <- tapply(p_theta,cumsum(!duplicated(round(theta,8))),sum)
  theta <- theta[!duplicated(round(theta,8))]

  E_z <- matrix(NA,length(y),length(theta))
  for (i in 1:length(theta)) {
    E_z[,i] <- log(p_theta[i])+dnorm(y, theta[i], se, log=TRUE)
  }
  E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x)))))

  rownames(E_z)<-row_names
  colnames(E_z)<-signif(theta,3)

  return(list(theta=theta, p_theta=p_theta, post_theta=E_z))
}

getmode <- function(v) {
  #returns mode from list v
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#data type agnostic
PlotClusterRank <- function(ClusterRank,xlab=NULL, maintitle=NULL) {
  # Creates a plot for a ClusterRank object
  #
  # Args:
  #   ClusterRank: a ClusterRank object, the output of ClusterRankBin, ClusterRankPois, ClusterRankNorm
  #   xlab: label for x axis
  #   maintitle: title for plot
  #
  # Returns:
  #     a visualization using the result of ClusterRankBin, ClusterRankPois or ClusterRankNorm.
  #     Shows ranks with clusters and confidence intervals of ranks.
  #
  post_df <- reshape2::melt(ClusterRank$posterior)
  post_df$cluster <- ClusterRank$ranked_table$cluster[match(post_df$Var1,ClusterRank$ranked_table$name)]
  post_df$p_cluster <- ClusterRank$ranked_table$p_cluster[match(post_df$Var1,ClusterRank$ranked_table$name)]

  return(ggplot2::ggplot(ClusterRank$ranked_table,aes(y=name,x=est,color=cluster,alpha=p_cluster))+
    ggplot2::geom_point(pch=3)+
    ggplot2::geom_point(aes(x=posteriorMean),pch=4)+
    ggplot2::geom_point(data=post_df,aes(y=Var1,x=as.numeric(Var2),color=cluster,size=value,alpha=value))+
    ggplot2::geom_errorbarh(aes(xmin=p_LCL,xmax=p_UCL),height=0)+
    ggplot2::scale_y_discrete("",limits=rev(levels(ClusterRank$ranked_table$name)))+
    ggplot2::scale_x_continuous(xlab,breaks=ClusterRank$theta[!duplicated(round(ClusterRank$theta,2))],
                     labels=round(ClusterRank$theta[!duplicated(round(ClusterRank$theta,2))],3),minor_breaks=ClusterRank$theta)+
    ggplot2::scale_color_manual(values=rep(RColorBrewer::brewer.pal(8,"Dark2"),1+floor(length(levels(ClusterRank$ranked_table$cluster))/8)))+
    ggplot2::scale_size_area(max_size=5)+
    ggplot2::scale_alpha(limits=c(0,1),range=c(0,1))+
    ggplot2::theme_bw()+
    ggplot2::guides(color=FALSE,size=FALSE,alpha=FALSE))
}

#change this to analog to Wilson CI:  https://www.ine.pt/revstat/pdf/rs120203.pdf
#(uses phat in denom) y/t = lambda hat
#var(lambda hat) = lT/T^2 = lambda/T. See binconf code for solving
WilsonHilfertyPoiCI <- function (x, ti, conf.level=0.95) {
  # Calculates approximate Poisson CI by Wilson & Hilferty (1931)
  #
  # Args:
  #   x: number of Poisson-distributed events
  #   ti: length of time. defaults to 1
  #   conf.level: desired converage for confidence intervals
  #
  # Returns:
  #     vector containing estimate, lower bound, upper bound
  #
  est = x/ti
  est1 = est+1
  alpha = 1 - conf.level
  alpha1 = 1-alpha/2
  Z_alpha1 = qnorm(alpha1)
  lower = est*(1-(1/(9*est)) - (Z_alpha1/(3*sqrt(est))))^3
  upper = est1*(1-(1/(9*est1)) + (Z_alpha1/(3*sqrt(est1))))^3
  return(c(est,lower, upper))
}


