# Statistical functions and simulations for Bernoulli and Poisson outcomes.
# Requires base R only. Run: Rscript non_gaussian_simulations.R
# Results are created in a results/ folder beside this script.
# Exact binary and count outcomes are used without normal approximation.
# All choices below are fixed before the evaluation simulations are generated.
args <- commandArgs(FALSE)
filearg <- sub("^--file=", "", args[grepl("^--file=", args)])
if(length(filearg)!=1L) stop("Run this script with Rscript.")
root <- dirname(normalizePath(gsub("~+~", " ", filearg, fixed=TRUE), mustWork=TRUE))
out <- file.path(root,"results")
dir.create(out,showWarnings=FALSE)
started <- Sys.time()
K <- 4L; alpha <- .05; M <- 50000L; N <- 100L; batch <- 10L; looks <- N/batch
models <- list(
  Bernoulli=list(null=.20,sparse=.30,dense=.25,single=.30,grid=c(.24,.28,.36,.52)),
  Poisson=list(null=.50,sparse=.70,dense=.60,single=.70,grid=c(.575,.65,.80,1.10)))
scenarios <- expand.grid(scenario=c("null","sparse","dense"),
  dependence=c("indep","shared"),stringsAsFactors=FALSE)

mass <- function(x,n,theta,model) {
  if(model=="Bernoulli") dbinom(x,n,theta) else dpois(x,n*theta)
}
cdf <- function(x,n,theta,model) {
  if(model=="Bernoulli") pbinom(x,n,theta) else ppois(x,n*theta)
}
tailprob <- function(x,n,theta,model) {
  if(model=="Bernoulli") pbinom(x-1,n,theta,lower.tail=FALSE) else
    ppois(x-1,n*theta,lower.tail=FALSE)
}
log_lr <- function(s,n,b,null,model) {
  if(model=="Bernoulli") s*log(b/null)+(n-s)*log((1-b)/(1-null)) else
    s*log(b/null)-n*(b-null)
}
uniform_cal <- function(p) {
  L <- -log(p); ans <- L; small <- L<1e-3
  ans[small] <- .5+L[small]/6+L[small]^2/24+L[small]^3/120+L[small]^4/720
  ans[!small] <- (expm1(L[!small])-L[!small])/L[!small]^2
  ans[p==0] <- Inf
  ans
}
e_mean <- function(s,n,bet,null,model) rowMeans(exp(log_lr(s,n,bet,null,model)))
e_mix <- function(s,n,bets,null,model) {
  ans <- numeric(nrow(s))
  for(b in bets) ans <- ans+e_mean(s,n,b,null,model)
  ans/length(bets)
}

# Exact integer boundaries: keep the total probability of having survived all
# looks at least 1-alpha*t/(K*T). Only surviving counts need to be stored.
# The Poisson cdf accounts for the whole upper tail; no tail is truncated away.
make_boundaries <- function(model,d) {
  previous <- 1; rows <- list()
  for(t in seq_len(looks)) {
    budget <- alpha*t/(K*looks)
    survive <- function(b) sum(previous*cdf(b-1-(seq_along(previous)-1),batch,d$null,model))
    b <- 0L
    while(survive(b)<1-budget) b <- b+1L
    counts <- 0:(b-1L)
    next_prob <- vapply(counts,function(s)
      sum(previous*mass(s-(seq_along(previous)-1),batch,d$null,model)),0.0)
    stopifnot(abs(sum(next_prob)-survive(b))<1e-12,
      1-sum(next_prob)<=budget+1e-12, b==0L || survive(b-1L)<1-budget)
    rows[[t]] <- data.frame(model=model,look=t,observations=t*batch,boundary=b,
      target_cumulative_error=budget,exact_cumulative_error=1-sum(next_prob))
    previous <- next_prob
  }
  do.call(rbind,rows)
}
boundaries <- do.call(rbind,lapply(names(models),function(m) make_boundaries(m,models[[m]])))
write.csv(boundaries,file.path(out,"non_gaussian_boundaries.csv"),row.names=FALSE)

simulate <- function(model,scenario,dependence,seed) {
  d <- models[[model]]; rho <- if(dependence=="indep") 0 else .5
  theta <- switch(scenario,null=rep(d$null,K),sparse=c(d$sparse,rep(d$null,K-1)),
    dense=rep(d$dense,K))
  set.seed(seed)
  S <- matrix(0,M,K)
  methods <- c("Bonf_looks","Exact_spending","E_LR_planned","E_mixture_planned",
    "E_LR_every","E_mixture_every")
  first <- matrix(N+1L,M,length(methods),dimnames=list(NULL,methods))
  bds <- boundaries$boundary[boundaries$model==model]
  for(n in seq_len(N)) {
    if(model=="Bernoulli") {
      U <- matrix(runif(M*K),M,K)
      if(rho>0) {
        shared <- runif(M)<rho; common <- runif(M)
        U[shared,] <- rep(common[shared],K)
      }
      X <- sweep(U,2,theta,"<=")
    } else {
      common <- if(rho>0) rpois(M,rho*d$null) else numeric(M)
      X <- matrix(rpois(M*K,rep(theta-rho*d$null,each=M)),M,K)+common
    }
    S <- S+X
    single <- e_mean(S,n,d$single,d$null,model)
    mixture <- e_mix(S,n,d$grid,d$null,model)
    for(m in c("E_LR_every","E_mixture_every")) {
      hit <- if(m=="E_LR_every") single>=1/alpha else mixture>=1/alpha
      first[hit & first[,m]==N+1L,m] <- n
    }
    if(n%%batch==0) {
      # Lookup over integer sufficient statistics avoids repeated cdf calls.
      p_lookup <- tailprob(0:max(S),n,d$null,model)
      P <- matrix(p_lookup[S+1L],M,K)
      hit <- cbind(Bonf_looks=apply(P,1,min)<=alpha/(K*looks),
        Exact_spending=apply(S,1,max)>=bds[n/batch],
        E_LR_planned=single>=1/alpha,E_mixture_planned=mixture>=1/alpha)
      for(m in colnames(hit)) first[hit[,m] & first[,m]==N+1L,m] <- n
    }
  }
  fixed <- cbind(Bonferroni=apply(P,1,min)<=alpha/K,
    E_half=rowMeans(.5/sqrt(P))>=1/alpha,
    E_uniform=rowMeans(uniform_cal(P))>=1/alpha,
    E_tuned=rowMeans(exp(log(alpha)+sqrt(2*log(1/alpha))*qnorm(P,lower.tail=FALSE)))>=1/alpha,
    E_LR=single>=1/alpha,E_mixture=mixture>=1/alpha)
  stopifnot(identical(unname(fixed[,"Bonferroni"]),
    rowMeans((P<=alpha/K)/(alpha/K))>=1/alpha))
  for(m in c("E_LR","E_mixture")) {
    stopifnot(all(first[,paste0(m,"_every")]<=first[,paste0(m,"_planned")]),
      all(!fixed[,m] | first[,paste0(m,"_planned")]<=N))
  }
  summarize <- function(rej,used,analysis) {
    count <- colSums(rej); rate <- count/M
    data.frame(model=model,analysis=analysis,scenario=scenario,dependence=dependence,
      method=colnames(rej),rejections=as.integer(count),replicates=M,rate=rate,
      mcse=sqrt(rate*(1-rate)/M),mean_observations=colMeans(used),
      mcse_observations=apply(used,2,sd)/sqrt(M),sum_observations=colSums(used),
      sum_squared_observations=colSums(used^2),row.names=NULL)
  }
  results <- rbind(summarize(fixed,matrix(N,M,ncol(fixed)),"fixed"),
    summarize(first<=N,pmin(first,N),"sequential"))
  design <- data.frame(model=model,scenario=scenario,dependence=dependence,rho=rho,
    null=d$null,theta1=theta[1],theta2=theta[2],theta3=theta[3],theta4=theta[4],
    single=d$single,grid1=d$grid[1],grid2=d$grid[2],grid3=d$grid[3],grid4=d$grid[4],
    replicates=M,observations=N,batch=batch,looks=looks,seed=seed)
  centered <- sweep(S,2,N*theta,"-")
  moments <- list()
  for(k in seq_len(K)) {
    moments[[length(moments)+1L]] <- data.frame(model=model,scenario=scenario,
      dependence=dependence,moment="mean",component1=k,component2=k,
      estimate=mean(S[,k]),expected=N*theta[k],mcse=sd(S[,k])/sqrt(M))
    for(j in k:K) {
      expected <- if(k==j) N*theta[k]*(if(model=="Bernoulli") 1-theta[k] else 1) else
        N*rho*(if(model=="Bernoulli") min(theta[k],theta[j])-theta[k]*theta[j] else d$null)
      product <- centered[,k]*centered[,j]
      moments[[length(moments)+1L]] <- data.frame(model=model,scenario=scenario,
        dependence=dependence,moment="central_second",component1=k,component2=j,
        estimate=mean(product),expected=expected,mcse=sd(product)/sqrt(M))
    }
  }
  list(results=results,design=design,moments=do.call(rbind,moments))
}
results <- list(); designs <- list(); moments <- list(); j <- 0L
for(model in names(models)) for(i in seq_len(nrow(scenarios))) {
  j <- j+1L; sc <- scenarios[i,]
  cat(model,sc$scenario,sc$dependence,"\n")
  ans <- simulate(model,sc$scenario,sc$dependence,202614L+j)
  results[[j]] <- ans$results; designs[[j]] <- ans$design; moments[[j]] <- ans$moments
}
results <- do.call(rbind,results); designs <- do.call(rbind,designs)
write.csv(results,file.path(out,"non_gaussian_results.csv"),row.names=FALSE)
write.csv(designs,file.path(out,"non_gaussian_design.csv"),row.names=FALSE)
write.csv(do.call(rbind,moments),file.path(out,"non_gaussian_moments.csv"),row.names=FALSE)
capture.output(print(results[,1:9],row.names=FALSE),
  file=file.path(out,"non_gaussian_summary.txt"))
writeLines(c(paste("Started:",format(started)),paste("Finished:",format(Sys.time())),
  "Fixed design: two models; null, sparse and dense alternatives; independent and shared-component observations.",
  "Each of twelve scenarios uses 50,000 paths. All methods within a scenario share observations.",
  "No normal approximations or Monte Carlo-calibrated rejection thresholds are used.",
  capture.output(sessionInfo())),file.path(out,"non_gaussian_run.txt"))
cat("Non-Gaussian simulations completed.\n")
