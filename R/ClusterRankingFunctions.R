#functions for ClusterRank

library(ggplot2)
library(reshape2)
library(clue)
library(Hmisc)
library(RColorBrewer)

ClusterRank <- function(ranked.table="character", data.type = "string", posterior="matrix",
                        cluster.thetas = "vector", thetas = "vector", smp="matrix",
                        smp.ord="matrix", post.dist.theta="matrix"){
  # Creates a ClusterRank class
  #
  # Args:
  #   ranked.table
  #   data.type
  #   posterior="matrix",
  #  cluster.thetas = "vector", thetas = "vector", smp="matrix",
  #  smp.ord="matrix", post.dist.theta="matrix"
  #
  me <- list(ranked.table=ranked.table, data.type = data.type, posterior=posterior,
             cluster.thetas = cluster.thetas, thetas = thetas,
             smp = smp, smp.ord = smp.ord, post.dist.theta = post.dist.theta)
  ## Set the name for the class
  class(me) <- append(class(me),"ClusterRank")
  return(me)
}


ClusterRankBin <- function(y,n,k=NULL, scale=identity,weighted=FALSE,n.iter=1000,n.samp=10000,row.names=NULL,
                           sig.digits=6, return.post=FALSE) {
  #in future, we could ask which kind of tie breaking we'd like, similarly to rank() function
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
  #   row.names: optional row names argument
  #   sig.digits: optional significant digits argument
  #   return.post: optional boolean for returning posteriors todo remove?
  #
  # Returns:
  #     ClusterRank object including ranked.table, posterior, cluster theta, thetas,
  #     smp, smp.ord, post_dist.theta
  #
  N <- length(y)
  if(missing(n)){
      stop("n required for binomial data")
  }
  if (is.null(row.names)){
    row.names <- as.character(seq(1:N))
  } else{
    row.names <- as.character(row.names)
  }

  if (!all.equal(length(y), length(n)) & !all.equal(length(y),length(row.names))){
    stop("y, n, and row.names must be vectors of the same length")
  }
  if (length(unique(y/n)) == 1){
    stop("All units have identical point mass at ", unique(y/n), " so these units cannot be clustered and ranked.")
  }

  npmle.res <- npmleBin(y=y,n=n,k=k,n.iter=n.iter,row.names=row.names, sig.digits=sig.digits)
  # samples from posterior distribution. Samples from cluster thetas with prob x = post.theta

  smp <- apply(npmle.res$post.theta,1,
               function(x,theta,n.samp)
                 sample(theta,n.samp,replace=TRUE,prob=x),
               theta=scale(npmle.res$theta),n.samp=n.samp)
  smp <- t(smp) #transposes
  smp.ord <- apply(smp,2,sort)
  #posterior distribution of theta
  post.dist.theta <- t(apply(round(smp.ord,sig.digits), 1, function(x, levels) table(factor(x, levels = levels))/length(x), levels=round(scale(npmle.res$theta),sig.digits)))
  if (weighted){ #inverse variance weighting
    wgt <- 1/pmax(.Machine$double.eps,apply(smp,1,var)) #if variance is zero, uses v small value to weight,
    #making it impossible to reassign a low variance estimate to wrong rank?
  }else {
    wgt <- rep(1,N)
  }
  #weighted square error loss rank optimization

  #RANKING ASSIGNMENT
  lossRank <- matrix(NA,N,N)
  for (i in 1:N) {
    for (j in 1:N) {
      lossRank[i,j] <- wgt[i] * mean((smp[i,]-smp.ord[j,])^2)
    }
  }
  totalRankLoss = sum(diag(lossRank)) #todo change for others
  rnk <- as.numeric(clue::solve_LSAP(lossRank))

  CI <- matrix(ncol = 3, nrow = N)
  colnames(CI) = c("PointEst", "Lower", "Upper")
  rownames(CI) = c(rep("", times = N))
  CI <- Hmisc::binconf(y,n) #creating confidence intervals

  #TODO check here for cluster ties, then add to pois, norm versions
  #should we actually just look at post.dist.theta?
  tab <- t(apply(smp.ord, 1, function(x, levels) table(factor(x, levels = levels)), levels = sort(unique(c(smp.ord)))))
  which(duplicated(tab, MARGIN = 1))
  if (anyDuplicated(tab, MARGIN = 1) != 0){
  #TODO flag this, give warning
    warning("Rankings ", duplicated(tab, MARGIN = 1), "are tied so these assignments are arbitrary.")
    for (tie in which(duplicated(tab, MARGIN = 1))){
      #12, 13, 14 replcae with 12, 12, 12 or with min, max, mean, random
      #the ord <- order(rnk, row.names)


      #we need to check if any of the ties are next to each other
      #(it should give two duplicate if there are more than one)
      #find the range of duplicates, then rearrange based on the sort of their CI[,1] p estimates
      stop("Ties exist in cluster assignment between ", tie, " and ", tie-1) #TODO
      #tie breaker using posterior means rnk must respect posterior means
      if ((CI[tie,1] > CI[tie-1,1]) && (rnk[tie] < rnk[tie-1])){ #todo does this catch all cases
        print(paste("Switching rank assignments of ", tie, " and ", tie-1))
        tmp <- rnk[tie]
        rnk[tie] <- rnk[tie-1]
        rnk[tie-1] <- tmp
      }
    }
  }

  #CLUSTER ASSIGNMENTS
  lossCluster <- matrix(NA,N,length(npmle.res$theta)) # square error loss cluster optimization
  for (i in 1:N) {
    for (j in 1:length(npmle.res$theta)) {
      lossCluster[i,j] <- mean((smp[i,]-c(scale(npmle.res$theta)[j]))^2)
    }
  }
  totalClusterLoss <- sum(lossCluster)
  cluster <- apply(lossCluster, 1, which.min)
  cluster <- factor(cluster)
  p_cluster <- npmle.res$post.theta[cbind(1:N,as.numeric(cluster))]
  levels(cluster) <- signif(npmle.res$theta,sig.digits) #cluster labels

  #DATA TABLE
  ord <- order(rnk)
  ranked_table <- data.frame(name=row.names,rank=rnk,cluster=factor(cluster),
                             y=y,n=n,est = CI[,1],
                             p_LCL=CI[,2],p_UCL=CI[,3],
                             posteriorMean=c(npmle.res$post.theta%*%npmle.res$theta),
                             p_cluster=p_cluster)
  ranked_table <- ranked_table[ord,]
  ranked_table$name <- factor(ranked_table$name,levels=ranked_table$name,ordered=TRUE)
  posterior <- npmle.res$post.theta[ord,]

  # if (return.post) {
  #   print(list(ranked_table=ranked_table,posterior=posterior,cluster.thetas=npmle.res$theta,thetas=npmle.res$p.theta, smp=smp,smp.ord=smp.ord, post.dist.theta=post.dist.theta))
  # } else{
  #   print(list(ranked_table=ranked_table,posterior=posterior,cluster.thetas=npmle.res$theta, thetas=npmle.res$p.theta))
  # }
  #TODO check this
  obj <- ClusterRank(ranked.table=ranked_table, data.type = "Binomial", posterior=posterior,
                     cluster.thetas=npmle.res$theta, thetas=npmle.res$p.theta,
                     smp=smp,smp.ord=smp.ord, post.dist.theta=post.dist.theta)
  #todo take posteriors out of main object. Have a method that takes the object and generates the posteriors
  return(obj)
}

ClusterRankPois <- function(y,ti=rep(1,length(y)),k=NULL,
                        scale=identity,weighted=FALSE,n.iter=1000,n.samp=10000,row.names=NULL, sig.digits=6, return.post=FALSE) {
  # Assigns ranks then clusters to each item in a list based on Poisson data. Calls npmlePois()
  #
  # Args:
  #   y: number of Poisson distributed events
  #   ti: time vector. Length of time. (can be person-time)
  #   k: number of starting clusters desired. Defaults to length(y)
  #   scale: scale for ranking
  #   weighted: boolean indicating if inverse variance weighted is used
  #   n.iter: iterations used in EM algorithm
  #   n.samp: number of samples from posterior distribution
  #   row.names: optional row names argument
  #
  # Returns:
  #     list including ranked_table, posterior, cluster thetas, thetas
  #
  N <- length(y)
  if (is.null(row.names)){
    row.names <- as.character(seq(1:N))
  } else{
    row.names <- as.character(row.names)
  }
  if (!all.equal(length(y), length(ti)) & !all.equal(length(y),length(row.names))){ #TODO there's probably a more elegant way
    stop("y, ti, and row.names must be vectors of the same length")
  }
  if (length(unique(y)) == 1){
    stop("All units have identical point mass at y = ", unique(y), " so these units cannot be clustered and ranked.")
  }

  npmle.res <- npmlePois(y=y,ti=ti,k=k,n.iter=n.iter,row.names=row.names, sig.digits=sig.digits)
  smp <- apply(npmle.res$post.theta,1, #samples from a centered version of npmle.res$theta with pr = post.theta
               function(x,theta,n.samp)
                 sample(theta,n.samp,replace=TRUE,prob=x),
               theta=scale(npmle.res$theta),n.samp=n.samp)
  smp <- t(smp)
  smp.ord <- apply(smp,2,sort)
  #posterior distribution of theta
  post.dist.theta <- t(apply(round(smp.ord,sig.digits), 1, function(x, levels) table(factor(x, levels = levels))/length(x), levels=round(scale(npmle.res$theta),sig.digits)))
  if (weighted) { #inverse variance weighting
    wgt <- 1/pmax(.Machine$double.eps,apply(smp,1,var)) #if variance is zero, uses v small value to weight,
    #making it impossible to reassign a low variance estimate to wrong rank, right?
  }
  else {
    wgt <- rep(1,N)
  }
  lossRank <- matrix(NA,N,N)
  for (i in 1:N) {
    for (j in 1:N) {
      lossRank[i,j] <- wgt[i] * mean((smp[i,]-smp.ord[j,])^2)
    }
  }
 #totalRankLoss = sum(diag(lossRank)) #todo check this calculation
  rnk <- as.numeric(clue::solve_LSAP(lossRank))
    # square error loss cluster optimization
  lossCluster <- matrix(NA,N,length(npmle.res$theta))
  for (i in 1:N) {
    for (j in 1:length(npmle.res$theta)) {
      lossCluster[i,j] <- mean((smp[i,]-c(scale(npmle.res$theta)[j]))^2)
    }
  }
  #totalClusterLoss <- sum(lossCluster)
  cluster <- apply(lossCluster, 1, which.min)
  cluster <- factor(cluster)
  p_cluster <- npmle.res$post.theta[cbind(1:N,as.numeric(cluster))]
  levels(cluster) <- signif(npmle.res$theta,sig.digits) #labels

  ord <- order(rnk)
  CI <- matrix(ncol = 3, nrow = N)
  colnames(CI) = c("PointEst", "Lower", "Upper")
  rownames(CI) = c(rep("", times = N))
  ests <- WilsonHilfertyPoiCI(y, ti, conf.level=0.95)
  CI[, "PointEst"] <- ests[1]
  CI[, "Lower"] <- ests[2]
  CI[, "Upper"] <- ests[3]
  ranked_table <- data.frame(name=row.names,rank=rnk,cluster=factor(cluster),
                             y=y,ti=ti,est = CI[,1],
                             p_LCL=CI[,2],p_UCL=CI[,3],
                             posteriorMean=c(npmle.res$post.theta%*%npmle.res$theta),
                             p_cluster=p_cluster)
  ranked_table <- ranked_table[ord,]
  ranked_table$name <- factor(ranked_table$name,levels=ranked_table$name,ordered=TRUE)
  posterior <- npmle.res$post.theta[ord,]

  # if (return.post) {
  #   return(list(ranked_table=ranked_table,posterior=posterior,theta=npmle.res$theta,thetas=npmle.res$p.theta, smp=smp,smp.ord=smp.ord, ord=ord))
  # } else{
  #   return(list(ranked_table=ranked_table,posterior=posterior,theta=npmle.res$theta, thetas=npmle.res$p.theta))
  # }
  obj <- ClusterRank(ranked.table=ranked_table, data.type = "Poisson", posterior=posterior,
                     cluster.thetas=npmle.res$theta, thetas=npmle.res$p.theta,
                     smp=smp,smp.ord=smp.ord, post.dist.theta=post.dist.theta)
  return(obj)
}

ClusterRankNorm <- function(y, se, k=NULL, scale=identity,
                            weighted=FALSE,n.iter=1000,n.samp=10000,row.names=NULL, sig.digits=6, return.post=FALSE) {
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
  #   row.names: optional row names argument
  #
  # Returns:
  #     list including ranked_table, posterior, theta, thetas
  #
  N <- length(y)
  if(c(missing(se))) {
    stop("se required for normal data")
  }
  if (is.null(row.names)){
    row.names <- as.character(seq(1:N))
  } else{
    row.names <- as.character(row.names)
  }
  if (!all.equal(length(y), length(se)) & !all.equal(length(y),length(row.names))){
    stop("y, se, and row.names must be vectors of the same length")
  }

  if (length(unique(y/se)) == 1){
    stop("All units have identical point mass at ", unique(y/se), " so these units are all in the same cluster and cannot be sensibly ranked.")
  }

  npmle.res <- npmleNorm(y=y, se=se, k=k,n.iter=n.iter,row.names=row.names, sig.digits=sig.digits)
  # samples from posterior distribution. Samples from cluster thetas with prob x = post.theta
  smp <- apply(npmle.res$post.theta,1,
               function(x,theta,n.samp)
                 sample(theta,n.samp,replace=TRUE,prob=x),
               theta=scale(npmle.res$theta),n.samp=n.samp)
  smp <- t(smp)
  smp.ord <- apply(smp,2,sort)
  #posterior distribution of theta
  post.dist.theta <- t(apply(round(smp.ord,sig.digits), 1, function(x, levels) table(factor(x, levels = levels))/length(x), levels=round(scale(npmle.res$theta),sig.digits)))
  if (weighted) { #inverse variance weighting
    wgt <- 1/pmax(.Machine$double.eps,apply(smp,1,var)) #if variance is zero, uses v small value to weight,
                            #making it impossible to reassign a low variance estimate to wrong cluster
  }
  else {
    wgt <- rep(1,N)
  }
  lossRank <- matrix(NA,N,N) #weighted square error loss rank optimization
  for (i in 1:N) {
    for (j in 1:N) {
      lossRank[i,j] <- wgt[i] * mean((smp[i,]-smp.ord[j,])^2)
    }
  }
  totalRankLoss = sum(diag(lossRank)) #TODO should be diag everywhere
  rnk <- as.numeric(clue::solve_LSAP(lossRank))
  # square error loss cluster optimization
  lossCluster <- matrix(NA,N,length(npmle.res$theta))
  for (i in 1:N) {
    for (j in 1:length(npmle.res$theta)) {
      lossCluster[i,j] <- mean((smp[i,]-c(scale(npmle.res$theta)[j]))^2)
    }
  }
  totalClusterLoss <- sum(lossCluster) #TODO check with Ron
  cluster <- apply(lossCluster, 1, which.min)
  cluster <- factor(cluster)
  p_cluster <- npmle.res$post.theta[cbind(1:N,as.numeric(cluster))]
  levels(cluster) <- signif(npmle.res$theta,sig.digits) #labels for clusters

  ord <- order(rnk)
  CI <- matrix(ncol = 3, nrow = N)
  colnames(CI) = c("PointEst", "Lower", "Upper")
  rownames(CI) = c(rep("", times = N))
  CI[, "PointEst"] <- y
  CI[, "Lower"] <- y - 1.96*se
  CI[, "Upper"] <- y + 1.96*se
  ranked_table <- data.frame(name=row.names,rank=rnk,cluster=factor(cluster),
                             y=y, se = se, est = CI[,1],
                             p_LCL=CI[,2],p_UCL=CI[,3],
                             posteriorMean=c(npmle.res$post.theta%*%npmle.res$theta),
                             p_cluster=p_cluster)
  ranked_table <- ranked_table[ord,]
  ranked_table$name <- factor(ranked_table$name,levels=ranked_table$name,ordered=TRUE)
  posterior <- npmle.res$post.theta[ord,]
  # if (return.post) {
  #   return(list(ranked_table=ranked_table,posterior=posterior,theta=npmle.res$theta,thetas=npmle.res$p.theta, smp=smp,smp.ord=smp.ord, post.dist.theta=post.dist.theta))
  # } else{
  #   return(list(ranked_table=ranked_table,posterior=posterior,theta=npmle.res$theta, thetas=npmle.res$p.theta))
  # }
  obj <- ClusterRank(ranked.table=ranked_table, data.type = "Normal", posterior=posterior,
                     cluster.thetas=npmle.res$theta, thetas=npmle.res$p.theta,
                     smp=smp,smp.ord=smp.ord, post.dist.theta=post.dist.theta)
  return(obj)
}

#NOTE: The npmle functions should only be called by the ClusterRank functions

npmleBin <- function(y,n,k=NULL,n.iter=1000,row.names,sig.digits) {
  # Estimates clusters nonparametrically using an EM algorithm. Calculates the
  # probability each item will be assigned to each cluster.
  # Called by ClusterRankBin()
  #
  # Args:
  #   y: number of binomial distributed events
  #   n: number of attempts
  #   k: initial number of clusters. Defaults to length(y)
  #   n.iter: iterations used in EM algorithm
  #   row.names: optional row names argument
  #
  # Returns:
  #     list including clusters thetas, prior distribution for each p.theta,
  #                   posterior probabilities for each item's assignment to each cluster
  #
  if (is.null(k)) {
    theta<-sort(y/n) #sorted probabilities
    k<-length(theta) #k = number of units to rank
  } else {
    theta <- seq(min(y/n),max(y/n),length=k) #starting mass points of F. We're estimating these, along with p.theta
  }
  p.theta <- rep(1/k,k) #probabilities of each mass point.

  if (length(unique(theta)) > 1){ #checks that thetas have > 1 point mass. If they are, skips EM step and proceeds
    E_z <- matrix(NA,length(y),k) #expected value of the probability that a county is in each of the k theta clusters
    #calculating the p that z_{ij} is equal to theta star j
    for (j in 1:n.iter) {
      for (i in 1:k) {
        #numerator
        E_z[,i] <- log(p.theta[i])+dbinom(y,n,theta[i],log=TRUE)
      }
      E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x))))) #normalizes (pseudo zs). This is the E part of EM alg
      p.theta <- apply(E_z,2,mean) #M-step: means over the matrix
      theta <- y%*%E_z/n%*%E_z #calculates optimal theta for each cluster
      #^this shouldnt create problems for a single point mass
    }
  } #end of length(y) > 1 cases
  if (length(unique(theta)) == 1){ #special case when nclusters = 1
    E_z <- rep(1, times = length(y))
  }

  #this reduces down to needed number of clusters (<=k)
  ord<-order(theta)
  theta<-c(theta[ord]) #sorts
  p.theta<-p.theta[ord] #sorts

  p.theta <- tapply(p.theta,cumsum(!duplicated(round(theta,sig.digits))),sum) #cumsum numbers clusters is ascending order. sums the pthetas that goes with each cluster. See pictures
  theta <- theta[!duplicated(round(theta,sig.digits))] #removes duplicate thetas

  E_z <- matrix(NA,length(y),length(theta))
  #final posterior probabilties for each county. Pr(county in cluster i)
  for (i in 1:length(theta)) {
    E_z[,i] <- log(p.theta[i])+dbinom(y,n,theta[i],log=TRUE)
  }

  #Creates post.theta data frame
  if (length(unique(theta)) == 1){ #special case when nclusters = 1
    E_z <- as.data.frame(exp(E_z-max(E_z))/sum(exp(E_z-max(E_z)))) #normalizes probabilities to avoid underflow
  } else {
    E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x))))) #normalizes probabilities to avoid underflow
  }
  rownames(E_z) <- row.names
  #TODO drop first element of vector
  colnames(E_z)<-signif(theta,sig.digits) #cluster names are rounded

  return(list(theta=theta, p.theta=p.theta, post.theta=E_z))
}

npmlePois <- function(y,ti,k=NULL,n.iter=1000,row.names, sig.digits) {
  # Estimates clusters nonparametrically using an EM algorithm. Calculates the
  # probability each item will be assigned to each lambda cluster.
  # Called by ClusterRankPois()
  #
  # Args:
  #   y: number of poisson distributed events
  #   n: number of people
  #   ti: length of time. defaults to 1
  #   k: initial number of clusters. Defaults to length(y)
  #   n.iter: iterations used in EM algorithm
  #   row.names: optional row names argument
  #
  # Returns:
  #     list including prior for theta, prior distribution for each p.theta,
  #                   posterior probabilities for each item's assignment to each cluster
  #
  if (is.null(k)) {
    theta<-sort(y/ti) #sorted probabilities.
    k<-length(theta) #number of clusters to start
  } else {
    theta <- seq(min(y/ti),max(y/ti),length=k)
  }
  p.theta <- rep(1/k,k) #evenly spaced probabilities between 0 and 1 for clusters

  E_z <- matrix(NA,length(y),k)
  for (j in 1:n.iter) {
    for (i in 1:k) {
      E_z[,i] <- log(p.theta[i])+dpois(y, ti*theta[i],log=TRUE) #E-step
    }
    E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x))))) #E_z will be
    p.theta <- apply(E_z,2,mean) #M-step
    theta <- y%*%E_z/ti%*%E_z
  }

  ord<-order(theta)
  theta<-c(theta[ord])
  p.theta<-p.theta[ord]

  p.theta <- tapply(p.theta,cumsum(!duplicated(round(theta,sig.digits))),sum)
  theta <- theta[!duplicated(round(theta, sig.digits))]

  E_z <- matrix(NA,length(y),length(theta))
  for (i in 1:length(theta)) {
    E_z[,i] <- log(p.theta[i])+dpois(y,ti*theta[i],log=TRUE)
  }
  if (length(theta) == 1){ #special case when nclusters = 1 when nclusters = 1, the matrix is transposed: #theta (rows) x items (col) TODO
    E_z <- apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x)))) #normalizes probabilities to avoid underflow
  } else {
    E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x))))) #normalizes probabilities to avoid underflow
  }
  rownames(E_z)<-row.names
  colnames(E_z)<-signif(theta,sig.digits)

  return(list(theta=theta, p.theta=p.theta, post.theta=E_z))
}

#key point with normal version:
#comes in as y, se, unknown theta, p.theta. We assume se is known here.
#theta i hat = see pics
npmleNorm <- function(y, se, k=NULL,n.iter=1000,row.names, sig.digits) {
  # Estimates clusters nonparametrically using an EM algorithm. Calculates the
  # probability each item will be assigned to each cluster.
  # Called by ClusterRankNorm()
  #
  # Args:
  #   y: mean for each item
  #   se: standard error for each mean y
  #   k: initial number of clusters. Defaults to length(y)
  #   n.iter: iterations used in EM algorithm
  #   row.names: optional row names argument
  #
  # Returns:
  #     list including prior for theta, prior distribution for each p.theta,
  #                   posterior probabilities for each item's assignment to each cluster
  #
  if (is.null(k)) {
    theta<-sort(y)
    k<-length(theta) #k = # clusters
  } else {
    theta <- seq(min(y),max(y),length=k)
  }
  p.theta <- rep(1/k, k)

  E_z <- matrix(NA,length(y),k)
  for (j in 1:n.iter) {
    for (i in 1:k) {
      E_z[,i] <- log(p.theta[i])+dnorm(y, theta[i], se, log=TRUE) #uses the se that comes in with data
    }
    E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x)))))
    p.theta <- apply(E_z,2,mean)
    theta <- ((y/se^2)%*%E_z)/((1/se^2)%*%E_z)
  }

  ord<-order(theta)
  theta<-c(theta[ord])
  p.theta<-p.theta[ord]

  p.theta <- tapply(p.theta,cumsum(!duplicated(round(theta,sig.digits))),sum)
  theta <- theta[!duplicated(round(theta,sig.digits))]

  E_z <- matrix(NA,length(y),length(theta))
  for (i in 1:length(theta)) {
    E_z[,i] <- log(p.theta[i])+dnorm(y, theta[i], se, log=TRUE)
  }
  if (length(theta) == 1){ #special case when nclusters = 1 when nclusters = 1, the matrix is transposed: #theta (rows) x items (col)
    E_z <- apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x)))) #normalizes probabilities to avoid underflow
  } else {
    E_z <- t(apply(E_z,1,function(x) exp(x-max(x))/sum(exp(x-max(x))))) #normalizes probabilities to avoid underflow
  }

  rownames(E_z)<-row.names
  colnames(E_z)<-signif(theta,sig.digits)

  return(list(theta=theta, p.theta=p.theta, post.theta=E_z))
}

getmode <- function(v) {
  #returns mode from list v
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#data type agnostic
#TODO make this take a ClusterRank class object only
PlotClusterRank <- function(ClusterRank,xlab="Ranking Measure", maintitle="Clustered Rankings") {
  # Creates a plot for a ClusterRank class object
  #
  # Args:
  #   ClusterRank: a ClusterRank class object,
  #     the output of ClusterRankBin, ClusterRankPois, ClusterRankNorm
  #   xlab: label for x axis
  #   maintitle: title for plot
  #
  # Returns:
  #     a visualization using the result of ClusterRankBin, ClusterRankPois or ClusterRankNorm.
  #     Shows ranks with clusters and confidence intervals of ranks.
  #

  post_df <- reshape2::melt(ClusterRank$posterior, as.is=TRUE)
  #name is an ordered factor, so the compare here isnt working
  post_df$cluster <- ClusterRank$ranked_table$cluster[match(post_df$Var1,ClusterRank$ranked_table$name)]
  post_df$p_cluster <- ClusterRank$ranked_table$p_cluster[match(post_df$Var1,ClusterRank$ranked_table$name)]

  #remove ggplot2:: here becuse we added ggplot2 to DESCRIPTION
  return(ggplot2::ggplot(ClusterRank$ranked_table,aes(y=name,x=est,color=cluster,alpha=p_cluster))+
    ggplot2::geom_point(pch=3)+
    ggplot2::geom_point(aes(x=posteriorMean),pch=4)+
    ggplot2::geom_point(data=post_df,aes(y=Var1,x=as.numeric(Var2),color=cluster, size=value,alpha=value))+
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


