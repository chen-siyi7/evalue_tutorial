
if (!requireNamespace("mvtnorm", quietly=TRUE)) stop("mvtnorm is required")
args <- commandArgs(trailingOnly=FALSE)
filearg <- sub("^--file=", "", args[grepl("^--file=", args)])
if(length(filearg)!=1L) stop("Run this script with Rscript.")
root <- dirname(normalizePath(gsub("~+~", " ", filearg, fixed=TRUE), mustWork=TRUE))
out <- file.path(root, "results")
dir.create(out, showWarnings=FALSE)
K <- 4L; alpha <- 0.05; grid <- c(.1,.3,.5,.7,.9)
Mfixed <- 100000L; Mseq <- 50000L
Tlooks <- 10L; nbatch <- 20L; Nmax <- Tlooks*nbatch
mu_tuned <- sqrt(2*log(1/alpha)/Nmax)
# Prespecified planning mean .10, with equal weights over a factor-of-four range.
# The alternative means .07 and .15 are not points on this grid.
stream_grid <- .10*2^(-2:2)
stream_mixture <- function(S,n) {
  ans <- numeric(nrow(S))
  for(b in stream_grid) ans <- ans+rowMeans(exp(b*S-n*b^2/2))
  ans/length(stream_grid)
}
# Exact one-sided group-sequential boundaries for equally spaced looks.
# Cumulative z-statistics at looks 1..T have correlation sqrt(i/j) under the null.
# Miwa's algorithm is deterministic, so this step leaves the RNG stream unchanged.
Rlooks <- outer(1:Tlooks,1:Tlooks,function(i,j) sqrt(pmin(i,j)/pmax(i,j)))
gs_const <- function(shape,level) {
  cross <- function(c) 1-mvtnorm::pmvnorm(upper=c*shape,corr=Rlooks,
                                          algorithm=mvtnorm::Miwa(steps=256))[1]
  uniroot(function(c) cross(c)-level,c(1,5),tol=1e-9)$root
}
pocock_shape <- rep(1,Tlooks); obf_shape <- sqrt(Tlooks/(1:Tlooks))
gs_bounds <- list(
  stouffer_pocock=gs_const(pocock_shape,alpha)*pocock_shape,
  stouffer_obf=gs_const(obf_shape,alpha)*obf_shape,
  bonf_pocock=gs_const(pocock_shape,alpha/K)*pocock_shape,
  bonf_obf=gs_const(obf_shape,alpha/K)*obf_shape)
sigma <- function(kind) {
  if (kind=="indep") return(diag(K))
  if (kind=="equicorr") return(.5*diag(K)+.5)
  s <- diag(K); s[1,2] <- s[2,1] <- s[3,4] <- s[4,3] <- .5; s
}
uniform_cal <- function(p) {
  L <- -log(p)
  ans <- L
  small <- abs(L)<1e-3
  ans[small] <- .5+L[small]/6+L[small]^2/24+L[small]^3/120+L[small]^4/720
  ans[!small] <- (expm1(L[!small])-L[!small])/L[!small]^2
  ans[p==0] <- Inf
  ans
}
summary_rows <- function(rej, scenario, dep) {
  count <- colSums(rej); n <- nrow(rej); rate <- count/n
  data.frame(scenario=scenario,dependence=dep,method=colnames(rej),
             rejections=as.integer(count),replicates=n,rate=rate,
             mcse=sqrt(rate*(1-rate)/n),row.names=NULL)
}
fixed_decisions <- function(Z,S) {
  P <- pnorm(Z,lower.tail=FALSE)
  Ehalf <- rowMeans(.5/sqrt(P))
  Egrid <- rowMeans(sapply(grid,function(k) rowMeans(k*P^(k-1))))
  Euniform <- rowMeans(uniform_cal(P))
  Estep <- rowMeans((P<=alpha)/alpha)
  EstepK <- rowMeans((P<=alpha/K)/(alpha/K))
  ELR <- rowMeans(exp(sqrt(2*log(1/alpha))*Z-log(1/alpha)))
  C <- rowMeans(tan((.5-pmin(pmax(P,1e-15),1-1e-15))*pi))
  cbind(Fisher=pchisq(-2*rowSums(log(P)),2*K,lower.tail=FALSE)<=alpha,
        Stouffer_IND=pnorm(rowSums(Z)/sqrt(K),lower.tail=FALSE)<=alpha,
        Stouffer_COR=pnorm(rowSums(Z)/sqrt(sum(S)),lower.tail=FALSE)<=alpha,
        Bonferroni=apply(P,1,min)<=alpha/K,
        Cauchy=.5-atan(C)/pi<=alpha,
        E_half=Ehalf>=1/alpha,E_grid=Egrid>=1/alpha,
        E_uniform=Euniform>=1/alpha,E_step_alpha=Estep>=1/alpha,
        E_step_alphaK=EstepK>=1/alpha,E_LR_alpha=ELR>=1/alpha)
}
# Deterministic checks for calibration, rejection-rule equivalence, and the worked example.
check_p <- c(.00001,.0034,.05,.5,.999999,1)
quadrature <- vapply(check_p,function(p) integrate(function(k) k*p^(k-1),0,1)$value,0.0)
stopifnot(max(abs(uniform_cal(check_p)-quadrature))<1e-7)
z_example <- c(2.71,2.12,1.88,1.05)
p_example <- pnorm(z_example,lower.tail=FALSE)
stopifnot(mean((p_example<=alpha)/alpha)==15)
stopifnot(abs(mean(exp(sqrt(2*log(1/alpha))*z_example-log(1/alpha)))-13.15142648)<1e-7)
set.seed(202605)
design <- list(
  list("null","indep",rep(0,K)),list("null","equicorr",rep(0,K)),
  list("null","block",rep(0,K)),list("sparse","indep",c(3,0,0,0)),
  list("sparse","equicorr",c(3,0,0,0)),list("dense","indep",rep(1,K)),
  list("dense","equicorr",rep(1,K)))
fixed <- do.call(rbind,lapply(design,function(d) {
  S <- sigma(d[[2]])
  Z <- mvtnorm::rmvnorm(Mfixed,mean=d[[3]],sigma=S)
  rej <- fixed_decisions(Z,S)
  stopifnot(identical(rej[,"Bonferroni"],rej[,"E_step_alphaK"]))
  summary_rows(rej,d[[1]],d[[2]])
}))
write.csv(fixed,file.path(out,"fixed_results.csv"),row.names=FALSE)

seq_run <- function(scenario,rho,mu) {
  S <- (1-rho)*diag(K)+rho
  L <- chol(S); cumS <- matrix(0,Mseq,K)
  log_hedge <- matrix(0,Mseq,K)
  methods <- c("Fisher_naive","Stouffer_naive","Stouffer_Bonf_looks",
               "Bonf_endpoints_looks","E_LR_020","E_LR_horizon","E_hedged_step",
               "Stouffer_Pocock","Stouffer_OBF","Bonf_Pocock","Bonf_OBF","E_LR_020_every","E_LR_mixture","E_LR_mixture_every")
  first <- matrix(Tlooks+1L,Mseq,length(methods),dimnames=list(NULL,methods))
  every_hit <- rep(FALSE,Mseq)
  first_observation <- rep(Nmax+1L,Mseq)
  mixture_hit <- rep(FALSE,Mseq)
  mixture_first <- rep(Nmax+1L,Mseq)
  for(t in seq_len(Tlooks)) {
    # Retain the original script's draws and RNG ordering exactly.
    batch <- matrix(0,Mseq,K)
    for(i in seq_len(nbatch)) {
      sample <- matrix(rnorm(Mseq*K),Mseq,K)%*%L
      batch <- batch+sweep(sample,2,mu,"+")
      # The same e-process checked after every observation (no extra draws).
      n_obs <- (t-1)*nbatch+i
      hit_now <- rowMeans(exp(.2*(cumS+batch)-n_obs*.2^2/2))>=1/alpha
      first_observation[hit_now & first_observation==Nmax+1L] <- n_obs
      every_hit <- every_hit | hit_now
      mix_now <- stream_mixture(cumS+batch,n_obs)>=1/alpha
      mixture_first[mix_now & mixture_first==Nmax+1L] <- n_obs
      mixture_hit <- mixture_hit | mix_now
    }
    cumS <- cumS+batch; n <- t*nbatch; Z <- cumS/sqrt(n)
    P <- pnorm(Z,lower.tail=FALSE)
    pf <- pchisq(-2*rowSums(log(P)),2*K,lower.tail=FALSE)
    ps <- pnorm(rowSums(Z)/sqrt(sum(S)),lower.tail=FALSE)
    E020 <- rowMeans(exp(.2*cumS-n*.2^2/2))
    Etuned <- rowMeans(exp(mu_tuned*cumS-n*mu_tuned^2/2))
    # Fresh batches, not repeated calibration of cumulative fixed-time p-values.
    pb <- pnorm(batch/sqrt(nbatch),lower.tail=FALSE)
    log_hedge <- log_hedge+log(.5+.5*(pb<=alpha)/alpha)
    Eh <- rowMeans(exp(log_hedge))
    zs <- rowSums(Z)/sqrt(sum(S)); zmax <- apply(Z,1,max)
    crossed <- cbind(pf<=alpha,ps<=alpha,ps<=alpha/Tlooks,
                     apply(P,1,min)<=alpha/(K*Tlooks),
                     E020>=1/alpha,Etuned>=1/alpha,Eh>=1/alpha,
                     zs>=gs_bounds$stouffer_pocock[t],zs>=gs_bounds$stouffer_obf[t],
                     zmax>=gs_bounds$bonf_pocock[t],zmax>=gs_bounds$bonf_obf[t],every_hit,
                     stream_mixture(cumS,n)>=1/alpha,mixture_hit)
    first[crossed & first==Tlooks+1L] <- t
  }
  rej <- first<=Tlooks
  res <- summary_rows(rej,scenario,if(rho==0) "indep" else "equicorr")
  used <- pmin(first,Tlooks)*nbatch
  used[,"E_LR_020_every"] <- pmin(first_observation,Nmax)
  stopifnot(all((first_observation<=Nmax)==rej[,"E_LR_020_every"]))
  stopifnot(all(used[,"E_LR_020_every"]<=used[,"E_LR_020"]))
  used[,"E_LR_mixture_every"] <- pmin(mixture_first,Nmax)
  stopifnot(all((mixture_first<=Nmax)==rej[,"E_LR_mixture_every"]))
  stopifnot(all(used[,"E_LR_mixture_every"]<=used[,"E_LR_mixture"]))
  res$mean_observations <- colMeans(used)
  res$mcse_observations <- apply(used,2,sd)/sqrt(Mseq)
  res$sum_observations <- colSums(used)
  res$sum_squared_observations <- colSums(used^2)
  res
}
set.seed(202606)
sequential <- rbind(seq_run("null",.5,rep(0,K)),seq_run("sparse",.5,c(.15,0,0,0)),
                    seq_run("dense",.5,rep(.07,K)),seq_run("null",0,rep(0,K)))
write.csv(sequential,file.path(out,"sequential_results.csv"),row.names=FALSE)

examples <- function(label,p,z=qnorm(p,lower.tail=FALSE)) {
  E <- cbind(half=.5/sqrt(p),grid=rowMeans(sapply(grid,function(k) k*p^(k-1))),
             uniform=uniform_cal(p),step_alpha=(p<=alpha)/alpha,
             step_alphaK=(p<=alpha/K)/(alpha/K),
             LR_alpha=exp(sqrt(2*log(1/alpha))*z-log(1/alpha)))
  data.frame(example=label,method=colnames(E),mean_e=colMeans(E),
             reciprocal_p=pmin(1,1/colMeans(E)),reject=colMeans(E)>=1/alpha)
}
worked <- rbind(examples("7.1",p_example,z_example),examples("7.2",c(1.3e-5,.020,.081,.184)))
write.csv(worked,file.path(out,"worked_examples.csv"),row.names=FALSE)
walk_z <- rbind(z_example,c(3.0,2.1,1.5,.8))
walk_e <- rowMeans(uniform_cal(pnorm(walk_z,lower.tail=FALSE)))
write.csv(data.frame(study=1:2,mean_e=walk_e,product_e=cumprod(walk_e),
                    reciprocal_p=pmin(1,1/cumprod(walk_e))),
          file.path(out,"workflow_example.csv"),row.names=FALSE)
write.csv(data.frame(z=z_example,p=p_example,step=(p_example<=alpha)/alpha,
                    LR=exp(sqrt(2*log(1/alpha))*z_example-log(1/alpha))),
          file.path(out,"example_7_1_components.csv"),row.names=FALSE)
write.csv(data.frame(look=1:Tlooks,as.data.frame(gs_bounds)),
          file.path(out,"group_sequential_boundaries.csv"),row.names=FALSE)

# Optional stopping versus a fixed-n z-test: one stream, X_i ~ N(delta,1), monitored
# after every observation; stop when the likelihood-ratio e-process reaches 1/alpha.
# The design guess is ratio*delta. "single" bets on the guess; "mixture" is a fixed,
# equally weighted mixture over guess*mix_grid (valid by the fixed-mixture theorem).
# Pilot paths estimate a cap targeting 80% power under delta.
# An independent simulation below evaluates each frozen cap.
mix_grid <- 2^(-2:2)
stop_run <- function(delta,ratio,method="single",M=50000L) {
  nNP <- ceiling(((qnorm(1-alpha)+qnorm(.8))/delta)^2); cap <- 10L*nNP
  mus <- ratio*delta*(if(method=="single") 1 else mix_grid)
  S <- numeric(M); tau <- rep(Inf,M); alive <- rep(TRUE,M)
  for(n in seq_len(cap)) {
    S[alive] <- S[alive]+rnorm(sum(alive),delta)
    if(method=="single") {
      mu <- mus
      hit <- alive & (mu*S-n*mu^2/2>=log(1/alpha))
    } else {
      Sa <- S[alive]
      Emix <- rowMeans(exp(outer(Sa,mus)-matrix(n*mus^2/2,length(Sa),length(mus),byrow=TRUE)))
      hit <- alive; hit[alive] <- Emix>=1/alpha
    }
    tau[hit] <- n; alive <- alive & !hit
    if(!any(alive)) break
  }
  stopifnot(mean(is.finite(tau))>=.8)
  nmax <- as.numeric(quantile(tau,.8,type=1)); stopped <- pmin(tau,nmax)
  # Distribution-free order-statistic interval for the population 80% quantile.
  ranks <- pmax(1L,pmin(M,qbinom(c(.025,.975),M,.8)+c(0L,1L)))
  cap_interval <- sort(tau)[ranks]
  data.frame(method=method,delta=delta,mu_over_delta=ratio,n_fixed=nNP,nmax80=nmax,
             nmax80_ratio=nmax/nNP,mean_stop=mean(stopped),mean_stop_ratio=mean(stopped)/nNP,
             mcse_ratio=sd(stopped)/sqrt(M)/nNP,p_stop_by_fixed_n=mean(tau<=nNP),
             p_reject_by_cap=mean(tau<=nmax),
             p_reject_by_horizon=mean(is.finite(tau)),simulation_horizon=cap,
             cap_ci_lower=cap_interval[1],cap_ci_upper=cap_interval[2],replicates=M)
}
set.seed(202607)
stopping <- do.call(rbind,lapply(c(.2,.1),function(d)
  do.call(rbind,lapply(c(1,.5,2),function(r) stop_run(d,r)))))
set.seed(202609)
stopping <- rbind(stopping,do.call(rbind,lapply(c(.2,.1),function(d)
  do.call(rbind,lapply(c(1,.5,2),function(r) stop_run(d,r,"mixture"))))))
stopping_design <- stopping
write.csv(stopping_design,file.path(out,"stopping_design.csv"),row.names=FALSE)

# Evaluate the estimated design caps on fresh paths, under the alternative and null.
# Standard errors below condition on the cap chosen by the independent pilot.
evaluate_cap <- function(delta,guess,method,cap,M=50000L) {
  mus <- guess*(if(method=="single") 1 else mix_grid)
  S <- numeric(M); used <- rep(as.integer(cap),M); alive <- rep(TRUE,M)
  for(n in seq_len(cap)) {
    active <- which(alive)
    if(!length(active)) break
    S[active] <- S[active]+rnorm(length(active),delta)
    if(method=="single") {
      hit <- mus*S[active]-n*mus^2/2>=log(1/alpha)
    } else {
      Emix <- rowMeans(exp(outer(S[active],mus)-
                     matrix(n*mus^2/2,length(active),length(mus),byrow=TRUE)))
      hit <- Emix>=1/alpha
    }
    crossed <- active[hit]; used[crossed] <- n; alive[crossed] <- FALSE
  }
  list(rejections=sum(!alive),mean_n=mean(used),mcse_n=sd(used)/sqrt(M),
       sum_n=sum(used),sum_n_squared=sum(used^2))
}
set.seed(202612)
stopping <- do.call(rbind,lapply(seq_len(nrow(stopping_design)),function(i) {
  d <- stopping_design[i,]; cap <- as.integer(d$nmax80); M <- 50000L
  alt <- evaluate_cap(d$delta,d$delta*d$mu_over_delta,d$method,cap,M)
  nul <- evaluate_cap(0,d$delta*d$mu_over_delta,d$method,cap,M)
  power <- alt$rejections/M; size <- nul$rejections/M
  data.frame(method=d$method,delta=d$delta,mu_over_delta=d$mu_over_delta,
    n_fixed=d$n_fixed,fixed_power=pnorm(sqrt(d$n_fixed)*d$delta-qnorm(1-alpha)),
    nmax80=cap,nmax80_ratio=cap/d$n_fixed,
    cap_ci_lower=d$cap_ci_lower,cap_ci_upper=d$cap_ci_upper,
    mean_stop=alt$mean_n,mean_stop_ratio=alt$mean_n/d$n_fixed,
    mcse_ratio=alt$mcse_n/d$n_fixed,sum_observations=alt$sum_n,
    sum_squared_observations=alt$sum_n_squared,
    rejections=alt$rejections,p_reject_by_cap=power,
    mcse_power=sqrt(power*(1-power)/M),null_rejections=nul$rejections,
    null_rate=size,mcse_null=sqrt(size*(1-size)/M),
    replicates=M,pilot_replicates=d$replicates)
}))
write.csv(stopping,file.path(out,"stopping_results.csv"),row.names=FALSE)

# Unplanned extension: the design above plans ten looks. Here every trial that has not
# rejected by look 10 is extended by ten more batches. Classical procedures keep the
# boundary value of their last planned look; the e-process keeps its threshold 1/alpha.
seq_extend <- function(scenario,rho,mu,looks=2L*Tlooks) {
  S <- (1-rho)*diag(K)+rho; L <- chol(S); cumS <- matrix(0,Mseq,K)
  methods <- c("Stouffer_Pocock","Stouffer_OBF","Bonf_Pocock","Bonf_OBF","E_LR_020","E_LR_020_every","E_LR_mixture","E_LR_mixture_every")
  first <- matrix(looks+1L,Mseq,length(methods),dimnames=list(NULL,methods))
  every_hit <- rep(FALSE,Mseq)
  mixture_hit <- rep(FALSE,Mseq)
  bnd <- function(b,t) b[min(t,Tlooks)]
  for(t in seq_len(looks)) {
    batch <- matrix(0,Mseq,K)
    for(i in seq_len(nbatch)) {
      batch <- batch+sweep(matrix(rnorm(Mseq*K),Mseq,K)%*%L,2,mu,"+")
      n_obs <- (t-1)*nbatch+i
      every_hit <- every_hit | rowMeans(exp(.2*(cumS+batch)-n_obs*.2^2/2))>=1/alpha
      mixture_hit <- mixture_hit | stream_mixture(cumS+batch,n_obs)>=1/alpha
    }
    cumS <- cumS+batch; n <- t*nbatch; Z <- cumS/sqrt(n)
    zs <- rowSums(Z)/sqrt(sum(S)); zmax <- apply(Z,1,max)
    crossed <- cbind(zs>=bnd(gs_bounds$stouffer_pocock,t),zs>=bnd(gs_bounds$stouffer_obf,t),
                     zmax>=bnd(gs_bounds$bonf_pocock,t),zmax>=bnd(gs_bounds$bonf_obf,t),
                     rowMeans(exp(.2*cumS-n*.2^2/2))>=1/alpha,every_hit,
                     stream_mixture(cumS,n)>=1/alpha,mixture_hit)
    first[crossed & first==looks+1L] <- t
  }
  dep <- if(rho==0) "indep" else "equicorr"
  planned <- summary_rows(first<=Tlooks,scenario,dep); planned$horizon <- "planned"
  extended <- summary_rows(first<=looks,scenario,dep); extended$horizon <- "extended"
  rbind(planned,extended)
}
set.seed(202610)
extension <- rbind(seq_extend("null",.5,rep(0,K)),seq_extend("sparse",.5,c(.15,0,0,0)),
                   seq_extend("dense",.5,rep(.07,K)),seq_extend("null",0,rep(0,K)))
write.csv(extension,file.path(out,"extension_results.csv"),row.names=FALSE)

# Benchmark with Wald's approximate SPRT boundaries, nominal alpha=.05 and beta=.2.
# Achieved error probabilities are recorded; these boundaries do not exactly match them.
sprt_run <- function(delta,truth,M=50000L,beta=.2) {
  nNP <- ceiling(((qnorm(1-alpha)+qnorm(.8))/delta)^2)
  A <- log((1-beta)/alpha); B <- log(beta/(1-alpha))
  S <- numeric(M); N <- rep(NA_real_,M); rej <- rep(FALSE,M); alive <- rep(TRUE,M); n <- 0L
  while(any(alive)) {
    n <- n+1L; S[alive] <- S[alive]+rnorm(sum(alive),truth)
    llr <- delta*S-n*delta^2/2
    up <- alive & llr>=A; down <- alive & llr<=B
    rej[up] <- TRUE; N[up|down] <- n; alive <- alive & !(up|down)
  }
  data.frame(delta=delta,truth=if(truth==0) "null" else "alternative",n_fixed=nNP,
             mean_n_ratio=mean(N)/nNP,mcse_ratio=sd(N)/sqrt(M)/nNP,
             rejections=sum(rej),reject=mean(rej),mcse_reject=sd(as.numeric(rej))/sqrt(M),
             upper_boundary=exp(A),lower_boundary=exp(B),replicates=M)
}
set.seed(202611)
sprt <- do.call(rbind,lapply(c(.2,.1),function(d) rbind(sprt_run(d,d),sprt_run(d,0))))
write.csv(sprt,file.path(out,"sprt_results.csv"),row.names=FALSE)

# Supplementary fixed-sample design with K=20 components: does the ordering among
# procedures valid under arbitrary dependence persist beyond K=4?
K20 <- 20L; M20 <- 50000L
k20_decisions <- function(Z,S) {
  P <- pnorm(Z,lower.tail=FALSE); k <- ncol(Z)
  C <- rowMeans(tan((.5-pmin(pmax(P,1e-15),1-1e-15))*pi))
  cbind(Bonferroni=apply(P,1,min)<=alpha/k,
        Stouffer_COR=pnorm(rowSums(Z)/sqrt(sum(S)),lower.tail=FALSE)<=alpha,
        Cauchy=.5-atan(C)/pi<=alpha,
        E_uniform=rowMeans(uniform_cal(P))>=1/alpha,
        E_step_alphaK=rowMeans((P<=alpha/k)/(alpha/k))>=1/alpha,
        E_LR_alpha=rowMeans(exp(sqrt(2*log(1/alpha))*Z-log(1/alpha)))>=1/alpha)
}
design20 <- list(list("null",rep(0,K20)),list("one",c(3.5,rep(0,K20-1))),
                 list("five",c(rep(2,5),rep(0,K20-5))),list("dense",rep(1,K20)))
set.seed(202608)
fixed20 <- do.call(rbind,lapply(c(0,.5),function(rho) do.call(rbind,lapply(design20,function(d) {
  S <- (1-rho)*diag(K20)+rho
  rej <- k20_decisions(mvtnorm::rmvnorm(M20,mean=d[[2]],sigma=S),S)
  stopifnot(identical(rej[,"Bonferroni"],rej[,"E_step_alphaK"]))
  summary_rows(rej,d[[1]],if(rho==0) "indep" else "equicorr")
}))))
write.csv(fixed20,file.path(out,"fixed_results_K20.csv"),row.names=FALSE)

# Optional continuation across studies: study k yields z_k ~ N(mu,1), independent.
# Studies are added one at a time, up to ten, and the analysis stops at the first rejection.
# Classical comparators include Bonferroni-adjusted Fisher at ten looks and
# Stouffer with exact Pocock and O'Brien-Fleming boundaries for ten studies.
# Fisher is also re-tested naively or used once at ten studies;
# products of per-study e-values stop when the running product reaches 1/alpha.
cross_run <- function(mu,M=100000L,Kmax=10L,lam=1) {
  Z <- matrix(rnorm(M*Kmax,mu),M,Kmax); P <- pnorm(Z,lower.tail=FALSE)
  cumsum_cols <- function(x) { for(k in seq_len(ncol(x))[-1]) x[,k] <- x[,k-1]+x[,k]; x }
  fisher <- cumsum_cols(-2*log(P))
  pf <- sapply(seq_len(Kmax),function(k) pchisq(fisher[,k],2*k,lower.tail=FALSE))
  stopifnot(Kmax==Tlooks)
  zs <- sweep(cumsum_cols(Z),2,sqrt(seq_len(Kmax)),"/")
  hits <- list(Fisher_repeated=pf<=alpha,
               Fisher_fixed10=cbind(matrix(FALSE,M,Kmax-1),pf[,Kmax]<=alpha),
               Fisher_Bonf_looks=pf<=alpha/Kmax,
               Stouffer_Pocock=sweep(zs,2,gs_bounds$stouffer_pocock,">="),
               Stouffer_OBF=sweep(zs,2,gs_bounds$stouffer_obf,">="),
               Product_fU=cumsum_cols(log(uniform_cal(P)))>=log(1/alpha),
               Product_half=cumsum_cols(log(.5)-.5*log(P))>=log(1/alpha),
               Product_LR=cumsum_cols(lam*Z-lam^2/2)>=log(1/alpha),
               Product_tuned=cumsum_cols(log(alpha)+sqrt(2*log(1/alpha))*Z)>=log(1/alpha))
  do.call(rbind,lapply(names(hits),function(m) {
    h <- hits[[m]]; first <- rep(NA_integer_,M)
    for(k in seq_len(Kmax)) first[is.na(first) & h[,k]] <- k
    rej <- !is.na(first); studies <- if(m=="Fisher_fixed10") rep(Kmax,M) else ifelse(rej,first,Kmax)
    data.frame(mu=mu,method=m,rejections=sum(rej),replicates=M,rate=mean(rej),
               mcse=sqrt(mean(rej)*(1-mean(rej))/M),mean_studies=mean(studies),
               mcse_studies=sd(studies)/sqrt(M),sum_studies=sum(studies),
               sum_squared_studies=sum(studies^2))
  }))
}
set.seed(202613)
continuation <- do.call(rbind,lapply(c(0,.5,1),cross_run))
write.csv(continuation,file.path(out,"continuation_results.csv"),row.names=FALSE)

# Expected log e-value per study (growth rate) under Z ~ N(mu,1): deterministic integrals.
growth <- function(logf,mu) integrate(function(z) logf(z)*dnorm(z-mu),mu-12,mu+12,rel.tol=1e-10)$value
lam_star <- sqrt(2*log(1/alpha))
growth_rates <- do.call(rbind,lapply(c(.5,1,lam_star/2,2),function(mu) data.frame(
  mu=mu,z_test_power=pnorm(qnorm(1-alpha)-mu,lower.tail=FALSE),
  continuous_mixture=growth(function(z) log(uniform_cal(pnorm(z,lower.tail=FALSE))),mu),
  half_power=growth(function(z) log(.5)-.5*pnorm(z,lower.tail=FALSE,log.p=TRUE),mu),
  threshold_tuned=log(alpha)+lam_star*mu, oracle_LR=mu^2/2)))
stopifnot(abs(growth_rates$threshold_tuned[3])<1e-12)
write.csv(growth_rates,file.path(out,"growth_rates.csv"),row.names=FALSE)

capture.output(sessionInfo(),file=file.path(out,"session_info.txt"))
capture.output(list(alpha=alpha,mu_horizon=mu_tuned,gs_bounds=gs_bounds,fixed=fixed,
                    sequential=sequential,examples=worked,stopping_design=stopping_design,
                    stopping=stopping,fixed20=fixed20,
                    extension=extension,sprt=sprt,
                    continuation=continuation,growth_rates=growth_rates),
               file=file.path(out,"simulation_summary.txt"))
cat("All numerical checks passed. Results written to",out,"\n")
