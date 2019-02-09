#functions for grouped ranking
#use this instead? devtools::use_package("dplyr") # Defaults to imports
library(ggplot2)
library(reshape2)
library(clue)
library(Hmisc)
library(RColorBrewer)
#To fix input problems
# ClusterRankBin() then within it calls
# ClusterRankPois()
# ClusterRankNorm()
ClusterRank <- function(y,n=NULL,se=NULL,ti=rep(1,length(y)),k=NULL,datatype,
        scale=identity,weighted=TRUE,n.iter=1000,n.samp=10000,row_names=NULL) {

            #take df for y,n OR y,se OR y,optional ti instead?
  # assigns ranks then clusters to each item in a list
  N <- length(y)
  if (datatype == "binomial"){
      if(missing(n)) {
        stop("n required for binomial data")
      }
      npmle_res <- npmle.bin(y=y,n=n,k=k,n.iter=n.iter,row_names=row_names)
  } else if (datatype == "poisson") {
      npmle_res <- npmle.pois(y=y,ti=ti,k=k,n.iter=n.iter,row_names=row_names)
  } else if (datatype == "normal"){
        if(c(missing(se))) {
        stop("se required for binomial data")
        }
      npmle_res <- npmle.norm(y=y, se=se, n=n, k=k,n.iter=n.iter,row_names=row_names)
  } else {
    stop("datatype must be binomial, poisson, or normal")
  }

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

  rnk <- as.numeric(clue::solve_LSAP(loss))
  grp <- match(apply(smp.ord,1,getmode),scale(npmle_res$theta))[rnk]
  #matches rank positions to groups using mode. The mode version minimizes indicator (see pic)
  #^ We could replace this with squared error diff to make things more consistent. See pic. TODO
  grp <- factor(grp)
  p_grp <- npmle_res$post_theta[cbind(1:N,as.numeric(grp))]
  levels(grp) <- signif(npmle_res$theta,3) #labels

  ord <- order(rnk)
  CI <- matrix(ncol = 3, nrow = N)
  colnames(CI) = c("PointEst", "Lower", "Upper")
  rownames(CI) = c(rep("", times = N))
  if (datatype == "binomial"){
    CI <- Hmisc::binconf(y,n) #creating confidence intervals
  } else if (datatype == "poisson") { #TODO update to score formula for better inference
    ests <- exactPoiCI(y, ti, conf.level=0.95)
    CI[, "PointEst"] <- ests[1] #as.numeric(poisson.test(y, T=ti, conf.level = 0.95)$estimate)
    CI[, "Lower"] <- ests[2] #poisson.test(y, T=ti, conf.level = 0.95)$conf.int[1]
    CI[, "Upper"] <- ests[3] #poisson.test(y, T=ti, conf.level = 0.95)$conf.int[2]
  } else if (datatype == "normal"){
    CI[, "PointEst"] <- mean(y)
    CI[, "Lower"] <- mean(y) - 1.96*se
    CI[, "Upper"] <- mean(y) + 1.96*se
  } else {
    stop("datatype must be binomial, poisson, or normal")
  }
    #TODO update this for normal make a separate function for each data type
  ranked_table <- data.frame(name=row_names,rank=rnk,group=factor(grp),
                             y=y,n=n,est = CI[,1], #p=y/n,
                             p_LCL=CI[,2],p_UCL=CI[,3],
                             pm=c(npmle_res$post_theta%*%npmle_res$theta),
                             p_grp=p_grp)
  ranked_table <- ranked_table[ord,]
  ranked_table$name <- factor(ranked_table$name,levels=ranked_table$name,ordered=TRUE)

  posterior <- npmle_res$post_theta[ord,]

  return(list(ranked_table=ranked_table,posterior=posterior,theta=npmle_res$theta,pr_theta=npmle_res$p_theta))
}

npmle.bin <- function(y,n,k=NULL,n.iter=1000,row_names=NULL) {
  #k is number of initial clusters
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
  #this reduces down to needed number of groups (<=k)
  ord<-order(theta)
  theta<-c(theta[ord]) #sorts
  p_theta<-p_theta[ord] #sorts

  p_theta <- tapply(p_theta,cumsum(!duplicated(round(theta,8))),sum) #cumsum numbers groups is ascending order. sums the pthetas that goes with each group. See pictures
  theta <- theta[!duplicated(round(theta,8))] #removes duplicate thetas

  E_z <- matrix(NA,length(y),length(theta))
  #final posterior probabilties for each county. Pr(county in group i)
  for (i in 1:length(theta)) {
    E_z[,i] <- log(p_theta[i])+dbinom(y,n,theta[i],log=TRUE)
  }
  E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x))))) #normalizes probabilities. subtracts max to avoid underflow

  rownames(E_z)<-row_names
  colnames(E_z)<-signif(theta,3) #group names are rounded

  return(list(theta=theta, p_theta=p_theta, post_theta=E_z))
  #return(prior for theta, prior for p_theta, posterior)
}

npmle.pois <- function(y,ti=rep(1,length(y)),k=NULL,n.iter=1000,row_names=NULL) {
  #y, persontime ti
  if (is.null(k)) {
    theta<-sort(y/ti) #sorted probabilities.
    k<-length(theta) #number of groups to start
  } else {
    theta <- seq(min(y/ti),max(y/ti),length=k)
  }
  p_theta <- rep(1/k,k) #evenly spaced probabilities between 0 and 1 for groups?

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
  print(dim(E_z))
  rownames(E_z)<-row_names
  colnames(E_z)<-signif(theta,3)

  return(list(theta=theta, p_theta=p_theta, post_theta=E_z))
}

#key point with normal version:
#comes in as y, se, unknown theta, p_theta. We assume se is known here.
#theta i hat = see pics
npmle.norm <- function(y, se, n, k=NULL,n.iter=1000,row_names=NULL) {
  if (is.null(k)) {
    theta<-sort(y)
    k<-length(theta) #k = # groups
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

#data type agnostic
getmode <- function(v) {
  #retrieves mode from list v
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#data type agnostic
PlotClusterRank <- function(ClusterRank,xlab=NULL, maintitle=NULL) {
  post_df <- reshape2::melt(ClusterRank$posterior)
  post_df$group <- ClusterRank$ranked_table$group[match(post_df$Var1,ClusterRank$ranked_table$name)]
  post_df$p_grp <- ClusterRank$ranked_table$p_grp[match(post_df$Var1,ClusterRank$ranked_table$name)]

  return(ggplot2::ggplot(ClusterRank$ranked_table,aes(y=name,x=p,color=group,alpha=p_grp))+
    ggplot2::geom_point(pch=3)+
    ggplot2::geom_point(aes(x=pm),pch=4)+
    ggplot2::geom_point(data=post_df,aes(y=Var1,x=as.numeric(Var2),color=group,size=value,alpha=value))+
    ggplot2::geom_errorbarh(aes(xmin=p_LCL,xmax=p_UCL),height=0)+
    ggplot2::scale_y_discrete("",limits=rev(levels(ClusterRank$ranked_table$name)))+
    ggplot2::scale_x_continuous(xlab,breaks=ClusterRank$theta[!duplicated(round(ClusterRank$theta,2))],
                     labels=round(ClusterRank$theta[!duplicated(round(ClusterRank$theta,2))],3),minor_breaks=ClusterRank$theta)+
    ggplot2::scale_color_manual(values=rep(RColorBrewer::brewer.pal(8,"Dark2"),1+floor(length(levels(ClusterRank$ranked_table$group))/8)))+
    ggplot2::scale_size_area(max_size=5)+
    ggplot2::scale_alpha(limits=c(0,1),range=c(0,1))+
    ggplot2::theme_bw()+
    ggplot2::guides(color=FALSE,size=FALSE,alpha=FALSE))
}

#change this to analog to Wilson interval for binomial
#(uses phat in denom) y/t = lambda hat
#var(lambda hat) = lT/T^2 = lambda/T. See binconf code for solving
exactPoiCI <- function (y, ti, conf.level=0.95) {
    alpha = 1 - conf.level
    est = y/ti
    upper <- 0.5 * qchisq(1-alpha/2, 2*(y/ti)+2)
    lower <- 0.5 * qchisq(alpha/2, 2*(y/ti))
    return(c(est,lower, upper))
}


