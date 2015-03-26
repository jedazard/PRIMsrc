##########################################################################################################################################
# PRIMsrc
##########################################################################################################################################

##########################################################################################################################################
# SURVIVAL INTERNAL SUBROUTINES
# (never to be called by end-user)
##########################################################################################################################################

##########################################################################################################################################
#################
# Usage         :
################
#                   cv.box.rep(x, times, status,
#                              B, K, arg,
#                              cvtype,
#                              probval, timeval,
#                              varsign, initcutpts,
#                              parallel, seed)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.box.rep <- function(x, times, status,
                       B, K, arg,
                       cvtype,
                       probval, timeval,
                       varsign, initcutpts,
                       parallel, seed) {

  CV.maxsteps <- numeric(B)
  CV.nsteps.lhr <- numeric(B)
  CV.nsteps.lrt <- numeric(B)
  CV.nsteps.cer <- numeric(B)
  CV.trace <- vector(mode="list", length=B)
  CV.boxind <- vector(mode="list", length=B)
  CV.boxcut <- vector(mode="list", length=B)
  CV.support <- vector(mode="list", length=B)
  CV.lhr <- vector(mode="list", length=B)
  CV.lrt <- vector(mode="list", length=B)
  CV.cer <- vector(mode="list", length=B)
  CV.time.bar <- vector(mode="list", length=B)
  CV.prob.bar <- vector(mode="list", length=B)
  CV.max.time.bar <- vector(mode="list", length=B)
  CV.min.prob.bar <- vector(mode="list", length=B)
  b <- 1
  k <- 0
  success <- TRUE
  while (b <= B) {
    cat("replicate : ", b, "\n", sep="")
    if (!parallel) {
      set.seed(seed[b])
      cat("seed : ", seed[b], "\n", sep="")
    }
    if (cvtype == "averaged") {
      CVBOX <- cv.ave.box(x=x, times=times, status=status,
                          K=K, arg=arg,
                          probval=probval, timeval=timeval,
                          varsign=varsign, initcutpts=initcutpts,
                          seed=seed[b])
    } else if (cvtype == "combined") {
      CVBOX <- cv.comb.box(x=x, times=times, status=status,
                           K=K, arg=arg,
                           probval=probval, timeval=timeval,
                           varsign=varsign, initcutpts=initcutpts,
                           seed=seed[b])
    } else if (cvtype == "none") {
      CVBOX <- cv.comb.box(x=x, times=times, status=status,
                           K=1, arg=arg,
                           probval=probval, timeval=timeval,
                           varsign=varsign, initcutpts=initcutpts,
                           seed=seed[b])
    } else {
      stop("Invalid CV type option \n")
    }
    if (!CVBOX$drop) {
      CV.maxsteps[b] <- CVBOX$cvfit$cv.maxsteps
      CV.nsteps.lhr[b] <- CVBOX$cvfit$cv.nsteps.lhr
      CV.nsteps.lrt[b] <- CVBOX$cvfit$cv.nsteps.lrt
      CV.nsteps.cer[b] <- CVBOX$cvfit$cv.nsteps.cer
      CV.trace[[b]] <- CVBOX$cvfit$cv.trace
      CV.boxind[[b]] <- CVBOX$cvfit$cv.boxind
      CV.boxcut[[b]] <- CVBOX$cvfit$cv.boxcut
      CV.support[[b]] <- CVBOX$cvfit$cv.stats$cv.support
      CV.lhr[[b]] <- CVBOX$cvfit$cv.stats$cv.lhr
      CV.lrt[[b]] <- CVBOX$cvfit$cv.stats$cv.lrt
      CV.cer[[b]] <- CVBOX$cvfit$cv.stats$cv.cer
      CV.time.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.time.bar
      CV.prob.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.prob.bar
      CV.max.time.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.max.time.bar
      CV.min.prob.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.min.prob.bar
      b <- b + 1
      k <- 0
    } else {
      cat("Could not find one step in at least one of the folds within replicate #", b ,". Retrying replicate with new seed\n", sep="")
      seed[b] <- seed[b] + B
      k <- k + 1
      if (k == B) {
        cat("Could not complete requested replications after ", B ," successive trials. Exiting replications\n", sep="")
        b <- B
        success <- FALSE
      }
    }
  }

  return(list("cv.maxsteps"=CV.maxsteps,
              "cv.nsteps.lhr"=CV.nsteps.lhr,
              "cv.nsteps.lrt"=CV.nsteps.lrt,
              "cv.nsteps.cer"=CV.nsteps.cer,
              "cv.trace"=CV.trace,
              "cv.boxind"=CV.boxind,
              "cv.boxcut"=CV.boxcut,
              "cv.support"=CV.support,
              "cv.lhr"=CV.lhr,
              "cv.lrt"=CV.lrt,
              "cv.cer"=CV.cer,
              "cv.time.bar"=CV.time.bar,
              "cv.prob.bar"=CV.prob.bar,
              "cv.max.time.bar"=CV.max.time.bar,
              "cv.min.prob.bar"=CV.min.prob.bar,
              "success"=success,
              "seed"=seed))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    cv.pval (x, times, status,
#                             cvtype,
#                             varsign, initcutpts,
#                             A, K, arg, obs.chisq,
#                             parallel, conf)
#
################
# Description   :
################
#                   Parallel computation of the cross-validated LRT p-value at a given peeling step
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.pval <- function(x, times, status,
                    cvtype,
                    varsign, initcutpts,
                    A, K, arg, obs.chisq,
                    parallel, conf) {

  if (!parallel) {
    null.chisq <- cv.null(x=x, times=times, status=status,
                          cvtype=cvtype,
                          varsign=varsign, initcutpts=initcutpts,
                          A=A, K=K, arg=arg)
  } else {
    if (conf$type == "SOCK") {
      cl <- makeCluster(spec=conf$names,
                        type=conf$type,
                        homogeneous=conf$homo,
                        outfile=conf$outfile,
                        verbose=conf$verbose)
    } else {
      cl <- makeCluster(spec=conf$cpus,
                        type=conf$type,
                        homogeneous=conf$homo,
                        outfile=conf$outfile,
                        verbose=conf$verbose)
    }
    clusterSetRNGStream(cl=cl, iseed=NULL)
    null.cl <- clusterCall(cl=cl, fun=cv.null,
                           x=x, times=times, status=status,
                           cvtype=cvtype,
                           varsign=varsign, initcutpts=initcutpts,
                           A=ceiling(A/conf$cpus), K=K, arg=arg)
    stopCluster(cl)
    null.chisq <- cbindlist(null.cl)
  }
  cvl <- nrow(null.chisq)
  pval <- numeric(cvl)
  for (l in 1:cvl) {
    pval[l] <- mean((null.chisq[l,] >= obs.chisq[l]), na.rm=TRUE)
  }

  return(pval)
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                   cv.null (x, times, status,
#                            cvtype,
#                            varsign, initcutpts,
#                            A, K, arg)
#
################
# Description   :
################
#                   Internal subroutine for computation of the cross-validated LRT null distribution at a given peeling step
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.null <- function(x, times, status,
                    cvtype,
                    varsign, initcutpts,
                    A, K, arg) {

  n <- nrow(x)
  null.chisq <- vector(mode="list", length=A)
  a <- 1
  while (a <= A) {
    cat("Permutation sample: ", a, "\n")
    perm.ind <- sample(x = 1:n, size = n, replace = FALSE, prob = NULL)
    perm.times <- times[perm.ind]
    perm.status <- status[perm.ind]
    if (cvtype == "averaged") {
      obj <- tryCatch({cv.ave.box(x=x, times=perm.times, status=perm.status, varsign=varsign, initcutpts=initcutpts, K=K, arg=arg, probval=NULL, timeval=NULL, seed=NULL)}, error=function(w){NULL})
      if (is.list(obj)) {
        null.chisq[[a]] <- obj$cvfit$cv.stats$cv.lrt
        a <- a + 1
      } else {
        cat("Permutation sample dropped... \n")
      }
    } else if (cvtype == "combined") {
      obj <- tryCatch({cv.comb.box(x=x, times=perm.times, status=perm.status, varsign=varsign, initcutpts=initcutpts, K=K, arg=arg, probval=NULL, timeval=NULL, seed=NULL)}, error=function(w){NULL})
      if (is.list(obj)) {
        null.chisq[[a]] <- obj$cvfit$cv.stats$cv.lrt
        a <- a + 1
      } else {
        cat("Permutation sample dropped... \n")
      }
    } else if (cvtype == "none") {
      obj <- tryCatch({cv.comb.box(x=x, times=perm.times, status=perm.status, varsign=varsign, initcutpts=initcutpts, K=1, arg=arg, probval=NULL, timeval=NULL, seed=NULL)}, error=function(w){NULL})
      if (is.list(obj)) {
        null.chisq[[a]] <- obj$cvfit$cv.stats$cv.lrt
        a <- a + 1
      } else {
        cat("Permutation sample dropped... \n")
      }
    } else {
      stop("Invalid CV type option \n")
    }
  }

  return(t(list2mat(null.chisq)))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                   cv.ave.box (x, times, status,
#                               probval, timeval,
#                               varsign, initcutpts,
#                               K, arg, seed)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.ave.box <- function(x, times, status,
                       probval, timeval,
                       varsign, initcutpts,
                       K, arg, seed) {

  n <- nrow(x)
  p <- ncol(x)

  fold.obj <- cv.ave.fold(x=x, times=times, status=status,
                          probval=probval, timeval=timeval,
                          varsign=varsign, initcutpts=initcutpts,
                          K=K, arg=arg, seed=seed)
  trace.list <- fold.obj$trace
  steps.list <- fold.obj$steps
  boxcut.list <- fold.obj$boxcut
  drop <- fold.obj$drop

  # Cross-validated minimum length from all folds
  CV.Lm <- min(fold.obj$nsteps)

  # Get the variable traces
  # Variable traces are first stacked and truncated in a matrix where folds are by rows and steps by columns
  CV.trace <- list2mat(list=trace.list, trunc=CV.Lm)
  CV.trace <- t(CV.trace)
  dimnames(CV.trace) <- list(paste("step", 0:(CV.Lm-1), sep=""), 1:K)

  # Truncate the cross-validated quantities from all folds to the same cross-validated length
  for (k in 1:K) {
    boxcut.list[[k]] <- boxcut.list[[k]][1:CV.Lm,]
    steps.list[[k]] <- steps.list[[k]][1:CV.Lm]
  }

  # Compute the averaged box statistics for each step from all the folds
  # Each entry or row signifies a step
  CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), colnames(x)))
  CV.support <- rep(NA, CV.Lm)
  names(CV.support) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.lhr <- rep(NA, CV.Lm)
  names(CV.lhr) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.lrt <- rep(NA, CV.Lm)
  names(CV.lrt) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.cer <- rep(NA, CV.Lm)
  names(CV.cer) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.time.bar <- rep(NA, CV.Lm)
  names(CV.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.prob.bar <- rep(NA, CV.Lm)
  names(CV.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.max.time.bar <- rep(NA, CV.Lm)
  names(CV.max.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.min.prob.bar <- rep(NA, CV.Lm)
  names(CV.min.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  for (l in 1:CV.Lm) {
    summincut <- matrix(NA, K, p)
    sumtime <- rep(NA, K)
    sumprob <- rep(NA, K)
    summaxtime <- rep(NA, K)
    summinprob <- rep(NA, K)
    sumlhr <- rep(NA, K)
    sumlrt <- rep(NA, K)
    sumci <- rep(NA, K)
    sumsupport <- rep(NA, K)
    for (k in 1:K) {
      outbounds <- steps.list[[k]][[l]]
      if (!is.null(outbounds)) {
        summincut[k,] <- boxcut.list[[k]][l,]
        sumlhr[k] <- steps.list[[k]][[l]][[1]]
        sumlrt[k] <- steps.list[[k]][[l]][[2]]
        sumci[k] <- steps.list[[k]][[l]][[3]]
        sumsupport[k] <- steps.list[[k]][[l]][[4]]
        sumtime[k] <- steps.list[[k]][[l]][[5]]
        sumprob[k] <- steps.list[[k]][[l]][[6]]
        summaxtime[k] <- steps.list[[k]][[l]][[7]]
        summinprob[k] <- steps.list[[k]][[l]][[8]]
      }
    }
    CV.boxcut[l, ] <- colMeans(summincut, na.rm=TRUE)
    CV.lhr[l] <- mean(sumlhr, na.rm=TRUE)
    CV.lrt[l] <- mean(sumlrt, na.rm=TRUE)
    CV.cer[l] <- mean(sumci, na.rm=TRUE)
    CV.support[l] <- mean(sumsupport, na.rm=TRUE)
    CV.time.bar[l] <- mean(sumtime, na.rm=TRUE)
    CV.prob.bar[l] <- mean(sumprob, na.rm=TRUE)
    CV.max.time.bar[l] <- mean(summaxtime, na.rm=TRUE)
    CV.min.prob.bar[l] <- mean(summinprob, na.rm=TRUE)
  }

  # Box peeling rules for each step
  CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), colnames(x))))
  for (j in 1:p) {
    if (varsign[j] > 0) {
      ss <- ">="
    } else {
      ss <- "<="
    }
    CV.rules[, j] <- paste(colnames(x)[j], ss, format(x=CV.boxcut[, j], digits=3, nsmall=3), sep="")
  }

  # Get the box membership indicator vector of all observations for each step from all the folds
  # Based on the corresponding averaged box over the folds
  CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
  for (l in 1:CV.Lm) {
    boxcut <- CV.boxcut[l, ] * varsign
    x.cut <- t(t(x) * varsign)
    x.ind <- t(t(x.cut) >= boxcut)
    CV.boxind[l,] <- (rowMeans(x.ind) == 1)  # Set as TRUE which observations are inside the box boudaries for all axes directions
  }
  rownames(CV.boxind) <- paste("step", 0:(CV.Lm-1), sep="")
  colnames(CV.boxind) <- rownames(x)

  # Applying the cross-validation criterion to the profiles
  # Cross-validated optimal length from all folds
  # By maximization of the LHR (between in and out box test samples)
  if (all(is.na(CV.lhr))) {
    CV.L.lhr <- NA
  } else {
    CV.L.lhr <- which.max(CV.lhr)
  }
  # By maximization of the LRT (between in and out box test samples)
  if (all(is.na(CV.lrt))) {
    CV.L.lrt <- NA
  } else {
    CV.L.lrt <- which.max(CV.lrt)
  }
  # By minimization of the CER (between predicted and observed inbox test samples survival times)
  if (all(is.na(CV.cer))) {
    CV.L.cer <- NA
  } else {
    CV.L.cer <- which.min(CV.cer)
  }

  # Box statistics for each step
  CV.stats <-  data.frame("cv.support"=CV.support,
                          "cv.lhr"=CV.lhr,
                          "cv.lrt"=CV.lrt,
                          "cv.cer"=CV.cer,
                          "cv.time.bar"=CV.time.bar,
                          "cv.prob.bar"=CV.prob.bar,
                          "cv.max.time.bar"=CV.max.time.bar,
                          "cv.min.prob.bar"=CV.min.prob.bar)
  rownames(CV.stats) <- paste("step", 0:(CV.Lm-1), sep="")

  # Create the return object 'CV.fit'
  CV.fit <- list("cv.nsteps.lhr"=CV.L.lhr,
                 "cv.nsteps.lrt"=CV.L.lrt,
                 "cv.nsteps.cer"=CV.L.cer,
                 "cv.maxsteps"=CV.Lm,
                 "cv.boxcut"=CV.boxcut,
                 "cv.rules"=CV.rules,
                 "cv.stats"=CV.stats,
                 "cv.trace"=CV.trace,
                 "cv.boxind"=CV.boxind)

  return(list("x"=x, "times"=times, "status"=status,
              "cvfit"=CV.fit, "drop"=drop, "seed"=seed))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                   cv.comb.box (x, times, status,
#                                probval, timeval,
#                                varsign, initcutpts,
#                                K, arg, seed)
#
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.comb.box <- function(x, times, status,
                        probval, timeval,
                        varsign, initcutpts,
                        K, arg, seed) {
  n <- nrow(x)
  p <- ncol(x)

  fold.obj <- cv.comb.fold(x=x, times=times, status=status,
                           varsign=varsign, initcutpts=initcutpts,
                           K=K, arg=arg, seed=seed)
  ord <- fold.obj$key
  times.list <- fold.obj$cvtimes
  status.list <- fold.obj$cvstatus
  trace.list <- fold.obj$trace
  boxind.list <- fold.obj$boxind
  boxcut.list <- fold.obj$boxcut

  # Cross-validated minimum length from all folds
  CV.Lm <- min(fold.obj$nsteps)

  # Concatenates the observations of test times and status from all folds
  # Re-ordered by initial order of observations
  CV.times <- unlist(times.list)[ord]
  CV.status <- unlist(status.list)[ord]

  # Get the variable traces
  # Variable traces are first stacked and truncated in a matrix where folds are by rows and steps by columns
  CV.trace <- list2mat(list=trace.list, trunc=CV.Lm)
  CV.trace <- t(CV.trace)
  dimnames(CV.trace) <- list(paste("step", 0:(CV.Lm-1), sep=""), 1:K)

  # Get the test box membership indicator vector of all observations for each step from all the folds
  # Based on the combined membership indicator vectors over the folds
  # Re-ordered by initial order of observations
  CV.boxind <- cbindlist(boxind.list, trunc=CV.Lm)[,ord]
  rownames(CV.boxind) <- paste("step", 0:(CV.Lm-1), sep="")
  colnames(CV.boxind) <- rownames(x)

  # Get the combined boxcut (truncated to the same cross-validated length) for each step from all the folds
  # using the circumscribing box to the conmbined test set in-box samples over all the folds
  CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), colnames(x)))
  tmparray <- list2array(list=boxcut.list, trunc=CV.Lm)
  for (l in 1:CV.Lm) {
    for (j in 1:p) {
        if (varsign[j] > 0) {
          CV.boxcut[l,j] <- min(x[CV.boxind[l,],j])
        } else {
          CV.boxcut[l,j] <- max(x[CV.boxind[l,],j])
        }
    }
  }

  # Box peeling rules for each step
  CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), colnames(x))))
  for (j in 1:p) {
    if (varsign[j] > 0) {
      ss <- ">="
    } else {
      ss <- "<="
    }
    CV.rules[, j] <- paste(colnames(x)[j], ss, format(x=CV.boxcut[, j], digits=3, nsmall=3), sep="")
  }

  # Compute the combined test box statistics from all folds for all steps, each entry or row signifies a step
  CV.support <- rep(NA, CV.Lm)
  names(CV.support) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.lhr <- rep(NA, CV.Lm)
  names(CV.lhr) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.lrt <- rep(NA, CV.Lm)
  names(CV.lrt) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.cer <- rep(NA, CV.Lm)
  names(CV.cer) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.time.bar <- rep(NA, CV.Lm)
  names(CV.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.prob.bar <- rep(NA, CV.Lm)
  names(CV.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.max.time.bar <- rep(NA, CV.Lm)
  names(CV.max.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.min.prob.bar <- rep(NA, CV.Lm)
  names(CV.min.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  timemat <- matrix(NA, nrow=CV.Lm, ncol=n)
  probmat <- matrix(NA, nrow=CV.Lm, ncol=n)
  ind.rem <- numeric(0)
  for (l in 1:CV.Lm) {
    boxind <- CV.boxind[l,]
    boxind1 <- 1*boxind
    if ((l == 1) && (sum(boxind, na.rm=TRUE) != 0)) {
      surv.fit <- survfit(Surv(CV.times[boxind], CV.status[boxind]) ~ 1)
      timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
      probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
      CV.lhr[l] <- 0
      CV.lrt[l] <- 0
      CV.cer[l] <- 1
      CV.support[l] <- 1
    } else if ((sum(boxind, na.rm=TRUE) != length(boxind[!is.na(boxind)])) && (sum(boxind, na.rm=TRUE) != 0)) {
      surv.fit <- survfit(Surv(CV.times[boxind], CV.status[boxind]) ~ 1)
      timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
      probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
      surv.formula <- (Surv(CV.times, CV.status) ~ 1 + boxind1)
      coxobj <- coxph(surv.formula, singular.ok=TRUE, iter.max=1)
      CV.lhr[l] <- coxobj$coef
      CV.lrt[l] <- survdiff(surv.formula, rho=0)$chisq
      predobj <- predict(object=coxobj, type="lp", reference="sample")
      CV.cer[l] <- rcorr.cens(x=predobj, S=Surv(CV.times, CV.status))['C Index']
      CV.support[l] <- mean(boxind, na.rm=TRUE)
    } else {
      timemat[l, ] <- NA
      probmat[l, ] <- NA
      CV.lhr[l] <- 0
      CV.lrt[l] <- 0
      CV.cer[l] <- 1
      CV.support[l] <- NA
      ind.rem <- c(ind.rem, l)
    }
  }
  if (length(ind.rem) != CV.Lm) {
    drop <- FALSE
    endobj <- endpoints (ind=ind.rem, timemat=timemat, probmat=probmat, timeval=timeval, probval=probval)
    time.bar <- endobj$time.bar
    prob.bar <- endobj$prob.bar
    max.time.bar <- endobj$max.time.bar
    min.prob.bar <- endobj$min.prob.bar
    for (l in 1:CV.Lm) {
      if (!(l %in% ind.rem)) {
        CV.time.bar[l] <- time.bar[l]
        CV.prob.bar[l] <- prob.bar[l]
        CV.max.time.bar[l] <-  max.time.bar[l]
        CV.min.prob.bar[l] <- min.prob.bar[l]
      } else {
        CV.time.bar[l] <- NA
        CV.prob.bar[l] <- NA
        CV.max.time.bar[l] <- NA
        CV.min.prob.bar[l] <- NA
      }
    }
  } else {
    cat("Dropped !\n", sep="")
    drop <- TRUE
    CV.time.bar <- rep(NA, CV.Lm)
    CV.prob.bar <- rep(NA, CV.Lm)
    CV.max.time.bar <- rep(NA, CV.Lm)
    CV.min.prob.bar <- rep(NA, CV.Lm)
  }

  # Applying the cross-validation criterion to the profiles
  # Cross-validated optimal length from all folds
  # By maximization of the LHR (between in and out box test samples)
  if (all(is.na(CV.lhr))) {
    CV.L.lhr <- NA
  } else {
    CV.L.lhr <- which.max(CV.lhr)
  }
  # By maximization of the LRT (between in and out box test samples)
  if (all(is.na(CV.lrt))) {
    CV.L.lrt <- NA
  } else {
    CV.L.lrt <- which.max(CV.lrt)
  }
  # By minimization of the CER (between predicted and observed inbox test samples survival times)
  if (all(is.na(CV.cer))) {
    CV.L.cer <- NA
  } else {
    CV.L.cer <- which.min(CV.cer)
  }

  # Box statistics for each step
  CV.stats <-  data.frame("cv.support"=CV.support,
                          "cv.lhr"=CV.lhr,
                          "cv.lrt"=CV.lrt,
                          "cv.cer"=CV.cer,
                          "cv.time.bar"=CV.time.bar,
                          "cv.prob.bar"=CV.prob.bar,
                          "cv.max.time.bar"=CV.max.time.bar,
                          "cv.min.prob.bar"=CV.min.prob.bar)
  rownames(CV.stats) <- paste("step", 0:(CV.Lm-1), sep="")

  # Create the return object 'CV.fit'
  CV.fit <- list("cv.nsteps.lhr"=CV.L.lhr,
                 "cv.nsteps.lrt"=CV.L.lrt,
                 "cv.nsteps.cer"=CV.L.cer,
                 "cv.maxsteps"=CV.Lm,
                 "cv.boxcut"=CV.boxcut,
                 "cv.rules"=CV.rules,
                 "cv.stats"=CV.stats,
                 "cv.trace"=CV.trace,
                 "cv.boxind"=CV.boxind)

  return(list("x"=x, "times"=times, "status"=status,
              "cvfit"=CV.fit, "drop"=drop, "seed"=seed))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    cv.ave.fold (x, times, status,
#                                 probval, timeval,
#                                 varsign, initcutpts,
#                                 K, arg, seed)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.ave.fold <- function(x, times, status,
                        probval, timeval,
                        varsign, initcutpts,
                        K, arg, seed) {

  drop <- FALSE
  folds <- cv.folds(n=nrow(x), K=K, seed=seed)
  steps <- vector(mode="list", length=K)
  boxcut <- vector(mode="list", length=K)
  trace <- vector(mode="list", length=K)
  nsteps <- numeric(K)

  for (k in 1:K) {
    cat("Overall fold : ", k, "\n", sep="")
    # Initialize training and test data
    if (K == 1) {
      traindata <- testdata <- x[folds$subsets[(folds$which == k)], , drop=FALSE]
      traintime <- testtime <- times[folds$subsets[(folds$which == k)]]
      trainstatus <- teststatus <- status[folds$subsets[(folds$which == k)]]
    } else {
      traindata <- x[folds$subsets[(folds$which != k)], , drop=FALSE]
      traintime <- times[folds$subsets[(folds$which != k)]]
      trainstatus <- status[folds$subsets[(folds$which != k)]]
      testdata <- x[folds$subsets[(folds$which == k)], , drop=FALSE]
      testtime <- times[folds$subsets[(folds$which == k)]]
      teststatus <- status[folds$subsets[(folds$which == k)]]
    }
    peelobj <- cv.ave.peel(traindata=traindata, trainstatus=trainstatus, traintime=traintime,
                           testdata=testdata, teststatus=teststatus, testtime=testtime,
                           probval=probval, timeval=timeval,
                           varsign=varsign, initcutpts=initcutpts, K=K, arg=arg, seed=seed)

    # Store the test set data/results from each fold (Note: observations of times and status from each fold are ordered in the list)
    steps[[k]] <- peelobj$steps
    nsteps[[k]] <- peelobj$nsteps
    boxcut[[k]] <- peelobj$boxcut
    trace[[k]] <- peelobj$trace
    drop <- (drop || peelobj$drop)
  }

  return(list("nsteps"=nsteps, "steps"=steps, "boxcut"=boxcut, "trace"=trace, "drop"=drop))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    cv.comb.fold (x, times, status,
#                                  varsign, initcutpts,
#                                  K, arg, seed)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.comb.fold <- function(x, times, status,
                         varsign, initcutpts,
                         K, arg, seed) {

  folds <- cv.folds(n=nrow(x), K=K, seed=seed)
  cvtimes <- vector(mode="list", length=K)
  cvstatus <- vector(mode="list", length=K)
  boxind <- vector(mode="list", length=K)
  boxcut <- vector(mode="list", length=K)
  trace <- vector(mode="list", length=K)
  nsteps <- numeric(K)

  for (k in 1:K) {
    cat("Overall fold : ", k, "\n", sep="")
    if (K == 1) {
      traindata <- testdata <- x[folds$subsets[(folds$which == k)], , drop=FALSE]
      traintime <- testtime <- times[folds$subsets[(folds$which == k)]]
      trainstatus <- teststatus <- status[folds$subsets[(folds$which == k)]]
    } else {
      traindata <- x[folds$subsets[(folds$which != k)], , drop=FALSE]
      traintime <- times[folds$subsets[(folds$which != k)]]
      trainstatus <- status[folds$subsets[(folds$which != k)]]
      testdata <- x[folds$subsets[(folds$which == k)], , drop=FALSE]
      testtime <- times[folds$subsets[(folds$which == k)]]
      teststatus <- status[folds$subsets[(folds$which == k)]]
    }
    peelobj <- cv.comb.peel(traindata=traindata, trainstatus=trainstatus, traintime=traintime,
                            testdata=testdata, teststatus=teststatus, testtime=testtime,
                            varsign=varsign, initcutpts=initcutpts, K=K, arg=arg, seed=seed)

    # Store the test set data/results from each fold
    # Note: the order of observations of times and status from each fold is kept in the list
    nsteps[k] <- peelobj$nsteps
    cvtimes[[k]] <- testtime
    cvstatus[[k]] <- teststatus
    boxind[[k]] <- peelobj$testindmat
    boxcut[[k]] <- peelobj$boxcut
    trace[[k]] <- peelobj$trace
  }

  return(list("nsteps"=nsteps, "key"=folds$key,
              "cvtimes"=cvtimes, "cvstatus"=cvstatus,
              "boxind"=boxind, "boxcut"=boxcut, "trace"=trace))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    cv.ave.peel (traindata, trainstatus, traintime,
#                                 testdata, teststatus, testtime,
#                                 probval, timeval,
#                                 varsign, initcutpts,
#                                 K, arg, seed)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.ave.peel <- function(traindata, trainstatus, traintime,
                        testdata, teststatus, testtime,
                        probval, timeval,
                        varsign, initcutpts,
                        K, arg, seed) {

  # Training the model
  peelobj <- peel.box(traindata=traindata, traintime=traintime, trainstatus=trainstatus,
                      varsign=varsign, initcutpts=initcutpts,
                      arg=arg, seed=seed)
  nsteps <- peelobj$nsteps

  # Compute the box statistics for all steps, each entry or row signifies a step
  steps.list <- vector(mode="list", length=nsteps)
  timemat <- matrix(NA, nrow=nsteps, ncol=nrow(testdata))
  probmat <- matrix(NA, nrow=nsteps, ncol=nrow(testdata))
  ind.rem <- numeric(0)
  for (l in 1:nsteps) {
    # Extract the rule and sign as one vector
    boxcut <- peelobj$boxcut[l, ] * varsign
    test.cut <- t(t(testdata) * varsign)
    test.ind <- t(t(test.cut) >= boxcut)
    # Set as TRUE which observations are TRUE for all covariates
    test.ind <- (rowMeans(test.ind) == 1)
    test.ind1 <- 1*test.ind
    if ((l == 1) && (sum(test.ind, na.rm=TRUE) != 0)) {
      lhr <- 0
      lrt <- 0
      cer <- 1
      support <- 1
      steps.list[[l]] <- c(lhr, lrt, cer, support)
      names(steps.list[[l]]) <- NULL
      surv.fit <- survfit(Surv(testtime[test.ind], teststatus[test.ind]) ~ 1)
      timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
      probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
    } else if ((sum(test.ind, na.rm=TRUE) != length(test.ind[!is.na(test.ind)])) && (sum(test.ind, na.rm=TRUE) != 0)) {
      surv.formula <- (Surv(testtime, teststatus) ~ 1 + test.ind1)
      coxobj <- coxph(surv.formula, singular.ok=TRUE, iter.max=1)
      lhr <- coxobj$coef
      lrt <- survdiff(surv.formula, rho=0)$chisq
      predobj <- predict(object=coxobj, type="lp", reference="sample")
      cer <- rcorr.cens(x=predobj, S=Surv(testtime, teststatus))['C Index']
      support <- mean(test.ind)
      steps.list[[l]] <- c(lhr, lrt, cer, support)
      names(steps.list[[l]]) <- NULL
      surv.fit <- survfit(Surv(testtime[test.ind], teststatus[test.ind]) ~ 1)
      timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
      probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
    } else {
      lhr <- 0
      lrt <- 0
      cer <- 1
      support <- NA
      steps.list[[l]] <- c(lhr, lrt, cer, support)
      names(steps.list[[l]]) <- NULL
      timemat[l, ] <- NA
      probmat[l, ] <- NA
      ind.rem <- c(ind.rem, l)
    }
  }
  if (length(ind.rem) != nsteps) {
    drop <- FALSE
    endobj <- endpoints (ind=ind.rem, timemat=timemat, probmat=probmat, timeval=timeval, probval=probval)
    time.bar <- endobj$time.bar
    prob.bar <- endobj$prob.bar
    max.time.bar <- endobj$max.time.bar
    min.prob.bar <- endobj$min.prob.bar
    for (l in 1:nsteps) {
      if (!(l %in% ind.rem)) {
        steps.list[[l]] <- c(steps.list[[l]], time.bar[l], prob.bar[l], max.time.bar[l], min.prob.bar[l])
      } else {
        steps.list[[l]] <- c(steps.list[[l]], NA, NA, NA, NA)
      }
    }
  } else {
    cat("Dropped !\n", sep="")
    drop <- TRUE
    for (l in 1:nsteps) {
      steps.list[[l]] <- rep(NA, 8)
    }
  }

  return(list("steps"=steps.list, "nsteps"=nsteps, "boxcut"=peelobj$boxcut, "trace"=peelobj$trace, "drop"=drop))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    cv.comb.peel (traindata, trainstatus, traintime,
#                                  testdata, teststatus, testtime,
#                                  varsign, initcutpts,
#                                  K, arg, seed)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.comb.peel <- function(traindata, trainstatus, traintime,
                         testdata, teststatus, testtime,
                         varsign, initcutpts,
                         K, arg, seed) {
  # Training the model
  peelobj <- peel.box(traindata=traindata, traintime=traintime, trainstatus=trainstatus,
                      varsign=varsign, initcutpts=initcutpts,
                      arg=arg, seed=seed)
  nsteps <- peelobj$nsteps

  # Create matrix that will represent the index of the test data that is within the box (i.e. for each step)
  test.ind.mat <- matrix(NA, nrow=nsteps, ncol=nrow(testdata))
  for (l in 1:nsteps) {
    # Extract the rule and sign as one vector
    boxcut <- peelobj$boxcut[l, ] * varsign
    test.cut <- t(t(testdata) * varsign)
    test.ind <- t(t(test.cut) >= boxcut)
    # Set as TRUE which observations are TRUE for all covariates
    test.ind.mat[l, ] <- (rowMeans(test.ind) == 1)
  }

  return(list("testindmat"=test.ind.mat, "nsteps"=nsteps, "boxcut"=peelobj$boxcut, "trace"=peelobj$trace))
}
##########################################################################################################################################





##########################################################################################################################################
#################
# Usage         :
################
#                    peel.box (traindata, traintime, trainstatus,
#                              varsign, initcutpts,
#                              arg, seed)
#
################
# Description   :
################
#                    Internal subroutine for fitting a Survival Bump Hunting (SBH) model. Search the first box of the recursive coverage 
#                    (outer) loop of our Patient Recursive Survival Peeling (PRSP) algorithm for the case of a survival (possibly censored) response.
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

peel.box <- function(traindata, traintime, trainstatus,
                     varsign, initcutpts,
                     arg, seed) {

  if (!is.null(seed))
    set.seed(seed)

  alpha <- NULL
  beta <- NULL
  minn <- NULL
  L <- NULL
  peelcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=arg, split=",")) ))

  digits <- getOption("digits")

  # Ensures that the training data is a numeric matrix
  traindata <- as.matrix(traindata)
  mode(traindata) <- "numeric"

  # Constants
  n <- nrow(traindata)                                   # Number of samples
  p <- ncol(traindata)                                   # Number of initially pre-selected covariates
  beta <- max(minn/n, beta)                              # Minimal box support thresholded to 10 points
  ncut <- ceiling(log(beta) / log(1 - (1/n)))            # Maximal number of peeling steps

  # Initializations of variable trace and box boundaries
  vartrace <- numeric(ncut)
  boxcut <- matrix(data=NA, nrow=ncut, ncol=p)

  # Initializations
  boxcutpts <- initcutpts
  boxmass <- 1
  boxes <- matrix(data=FALSE, nrow=n, ncol=p)            # Initial logical matrix of box membership indicator by dimension
  sel <- rep(TRUE, n)                                    # Initial logical vector of in-box samples
  xsel <- traindata                                      # Initial selection of samples from training data
  varpeel <- (apply(traindata, 2, "var") > 10^(-digits)) # Initial selection of covariates for peeling
  continue <- TRUE
  if (!(is.null(L))) {
    switch <- 1
  } else {
    L <- 1
    switch <- 0
  }
  l <- 0
  lhrlj <- matrix(NA, ncut, p)
  lrtlj <- matrix(NA, ncut, p)

  while ((boxmass >= beta) & (l*switch < L) & (continue)) {
    l <- l + 1
    xsign <- t(t(xsel) * varsign)

    # Potential cutpts by dimension
    cutpts.sign <- updatecut(x=xsign, fract=alpha)
    cutpts <- cutpts.sign * varsign

    # Update box membership indicator by dimension
    boxes <- as.matrix(t((t(traindata) * varsign) >= as.vector(cutpts.sign)) & sel)

    vmd <- rep(NA, p)
    for (j in 1:p) {
      boxes1j <- 1 * boxes[,j]
      if ((sum(boxes1j) != length(boxes1j)) && (sum(boxes1j) != 0)) {
        # Rate of increase of LHR (between in and out box)
        if (peelcriterion == "hr") {
          lhrlj[l,j] <- coxph(Surv(traintime, trainstatus) ~ 1 + boxes1j, singular.ok=TRUE, iter.max=1)$coef
          if (l == 1) {
            vmd[j] <- (lhrlj[l,j] - 0) / (1 - mean(boxes1j))
          } else {
            vmd[j] <- (lhrlj[l,j] - lhrlj[l-1,j]) / (boxmass - mean(boxes1j))
          }
        # Rate of increase of LRT (between in and out box)
        } else if (peelcriterion == "lr") {
          lrtlj[l,j] <- survdiff(Surv(traintime, trainstatus) ~ 1 + boxes1j, rho=0)$chisq
          if (l == 1) {
            vmd[j] <- (lrtlj[l,j] - 0) / (1 - mean(boxes1j))
          } else {
            vmd[j] <- (lrtlj[l,j] - lrtlj[l-1,j]) / (boxmass - mean(boxes1j))
          }
        } else {
          stop("Invalid peeling criterion \n")
        }
      } else {
        varpeel[j] <- FALSE
      }
    }

    # If the previous attempted peeling succeeded
    if (sum(varpeel) > 0 && (!is.empty(vmd[(!is.nan(vmd)) & (!is.infinite(vmd)) & (!is.na(vmd))]))) {
      # Maximizing the rate of increase of LHR or LRT (peeling criterion).
      # Only one variable (the first one in rank) is selected in case of ties
      varj <- which(vmd == max(vmd[(!is.nan(vmd)) & (!is.infinite(vmd))], na.rm=TRUE))[1]
      # Updating
      sel <- boxes[, varj, drop=TRUE]
      boxmass <- mean(1 * sel)
      xsel <- traindata[sel, ,drop=FALSE]
      varpeel <- (apply(xsel, 2, "var") > 10^(-digits))
      boxcutpts[varj] <- cutpts[varj]
      # Saving trained box quantities of interest for the current peeling step
      boxcut[l, ] <- boxcutpts
      vartrace[l] <- varj
    # Else exiting the loop and decrementing the peeling step number since the last attempted failed in that case
    } else {
      continue <- FALSE
      l <- l - 1
    }
  }

  if (l == 0) {
    # Taking the first step box covering all the data
    boxcut <- as.matrix(initcutpts)
    vartrace <- 0
  } else if (l >= 1) {
    # Prepending the first step box covering all the data
    boxcut <- rbind(initcutpts, boxcut[1:l, , drop=FALSE])
    vartrace <- c(0, vartrace[1:l])
  }

  rownames(boxcut) <- paste("step", 0:l, sep="")
  colnames(boxcut) <- colnames(traindata)
  names(vartrace) <- paste("step", 0:l, sep="")

  # Returning the final results, considering the starting point as step #0
  return(list("nsteps"=l+1,
              "boxcut"=boxcut,
              "trace"=vartrace))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    cv.folds (n, K, seed=NULL)
#
################
# Description   :
################
#                    Split n observations into K groups to be used for K-fold cross-validation.
#                    K should thereby be chosen such that all groups are of approximately equal size.
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.folds <- function (n, K, seed=NULL) {

  if (!is.null(seed))
    set.seed(seed)

  n <- round(rep(n, length.out = 1))
  if (!isTRUE(n > 0))
    stop("'n' must be positive")
  K <- round(rep(K, length.out = 1))
  if (!isTRUE((K >= 1) && K <= n))
    stop(paste("'K' outside allowable range {1,...,", n, "} \n", sep=""))
  if (K == 1) {
    subs <- seq_len(n)
  } else if (K == n) {
    subs <- seq_len(n)
  } else {
    subs <- sample(n)
  }
  w <- rep(seq_len(K), length.out=n)
  ord <- numeric(0)
  for (k in 1:K) {
    ord <- c(ord, subs[(w == k)])
  }
  key <- pmatch(x=1:n, table=ord)
  folds <- list(n=n, K=K, subsets=subs, which=w, key=key, seed=seed)

  return(folds)
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    endpoints (ind, timemat, probmat, timeval, probval)
#
################
# Description   :
################
#                    Return the maximum time value and the corresponding minimum survival probability
#                    for every box of the trajectory (box peeling sequence). Also return the time value
#                    and/or the corresponding survival probability depending on what is specified by the user.
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

endpoints <- function(ind, timemat, probmat, timeval, probval) {

  N <- nrow(timemat) #N <- nrow(probmat)
  L <- N-length(ind)
  if (!(is.empty(ind))) {
    timemat <- timemat[-ind, , drop=FALSE]
    probmat <- probmat[-ind, , drop=FALSE]
  }
  min.prob.bar <- apply(probmat, 1, min, na.rm=TRUE)
  max.time.bar <- apply(timemat, 1, max, na.rm=TRUE)
  if (is.null(probval) && is.null(timeval)) {
    prob.bar <- rep(NA, L)
    time.bar <- rep(NA, L)
  } else if (!is.null(probval)) {
    prob.bar <- rep(probval, L)
    ind.probmat <- (probmat <= probval)
    ind.probmat[is.na(ind.probmat)] <- FALSE
    time.bar <- numeric(L)
    for (l in 1:L) {
      if (probval >= min.prob.bar[l]) {
        time.bar[l] <- min(timemat[l,which(ind.probmat[l,,drop=TRUE])])
      } else {
        time.bar[l] <- NA
      }
    }
  } else if (!is.null(timeval)) {
    time.bar <- rep(timeval, L)
    ind.timemat <- (timemat >= timeval)
    ind.timemat[is.na(ind.timemat)] <- FALSE
    prob.bar <- numeric(L)
    for (l in 1:L) {
      if (timeval <= max.time.bar[l]) {
        prob.bar[l] <- max(probmat[l,which(ind.timemat[l,,drop=TRUE])])
      } else {
        prob.bar[l] <- NA
      }
    }
  }

  return(list("time.bar"=time.bar,
              "prob.bar"=prob.bar,
              "max.time.bar"=max.time.bar,
              "min.prob.bar"=min.prob.bar))

}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    updatecut (x, fract)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

updatecut <- function(x, fract) {

  p <- dim(x)[2]
  cutpts <- apply(x, 2, "quantile", type=7, probs=fract)

  for (j in 1:p) {
    xunique <- sort(unique(x[, j]))
    if (length(xunique) == 1) {
      cutpts[j] <- min(xunique)
    } else {
      if (isTRUE(all.equal(as.single(cutpts[j]), as.single(min(xunique))))) {
        cutpts[j] <- xunique[2]
      }
    }
  }

  return(cutpts)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   lapply.array (X, trunc=NULL, sub=NULL, fill=NA, MARGIN=1:2, FUN, ...)
#
################
#
################
# Description   :
################
#                   Compute FUN element-wise on entries of a list of matrices (even of different row or column numbers)
#
#
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

lapply.array <- function (X, trunc=NULL, sub=NULL, fill=NA, MARGIN=1:2, FUN, ...) {
  x <- list2array(list=X, sub=sub, trunc=trunc, fill=fill)
  return(apply(X=x, MARGIN=MARGIN, FUN=FUN, ...))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   lapply.mat (X, trunc=NULL, sub=NULL, fill=NA, MARGIN=2, FUN, ...)
#
################
#
################
# Description   :
################
#                   Compute FUN element-wise on entries of a list of vectors (even of different lengths)
#
#
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

lapply.mat <- function (X, trunc=NULL, sub=NULL, fill=NA, MARGIN=2, FUN, ...) {
  x <- list2mat(list=X, sub=sub, trunc=trunc, fill=fill)
  return(apply(X=x, MARGIN=MARGIN, FUN=FUN, ...))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   list2array (list, trunc=NULL, sub=NULL, fill=NA)
#
################
#
################
# Description   :
################
#                   Internal function to bind a list of matrices
#                   (even of different row or column numbers)
#                   by third dimension of a single 3D array
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

list2array <- function (list, trunc=NULL, sub=NULL, fill=NA) {
  if (!is.empty(list)) {
    if (is.null(sub)) {
      my.list <- list
    } else {
      L <- length(list)
      my.list <- vector(mode="list", length=L)
      for (i in 1:L) {
        my.list[[i]] <- list[[i]][[sub]]
      }
    }
    max.row <- max(sapply(my.list, nrow))
    max.col <- max(sapply(my.list, ncol))
    corrected.list <- lapply(my.list, function(x) {rbind(x, matrix(data=NA, nrow=max.row - nrow(x), ncol=ncol(x)))})
    corrected.list <- lapply(corrected.list, function(x) {cbind(x, matrix(data=fill, nrow=nrow(x), ncol=max.col - ncol(x)))})
    my.array <- array(data=fill, dim=c(nrow(corrected.list[[1]]), ncol(corrected.list[[1]]), length(corrected.list)))
    for(i in 1:length(corrected.list)) {
      my.array[,,i] <- corrected.list[[i]]
    }
    if (!is.null(trunc)) {
      if (trunc == "min") {
        trunc <- min(sapply(my.list, length))
      }
      my.array <- my.array[1:trunc,,,drop=FALSE]
    }
  } else {
    my.array <- array(data=fill, dim=c(0,0,0))
  }
  return(my.array)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   list2mat (list, trunc=NULL, sub=NULL, fill=NA)
#
################
#
################
# Description   :
################
#                   Internal function to bind a list of vectors
#                   (even of different length) by rows into a single matrix
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

list2mat <- function (list, trunc=NULL, sub=NULL, fill=NA) {
  if (!is.empty(list)) {
    if (is.null(sub)) {
      my.list <- list
    } else {
      L <- length(list)
      my.list <- vector(mode="list", length=L)
      for (i in 1:L) {
        my.list[[i]] <- list[[i]][[sub]]
      }
    }
    max.len <- max(sapply(my.list, length))
    corrected.list <- lapply(my.list, function(x) {c(x, rep(fill, max.len - length(x)))})
    my.mat <- do.call(rbind, corrected.list)
    if (!is.null(trunc)) {
      if (trunc == "min") {
        trunc <- min(sapply(my.list, length))
      }
      my.mat <- my.mat[, 1:trunc, drop=FALSE]
    }
  } else {
    my.mat <- matrix(data=fill, nrow=0, ncol=0)
  }
  return(my.mat)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   cbindlist (list, trunc)
#
################
#
################
# Description   :
################
#                   Internal function to bind a list of matrices by columns
#                   (even of different number of rows) into a single matrix
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cbindlist <- function(list, trunc) {
  if (!is.empty(list)) {
    max.row <- max(sapply(list, nrow))
    corrected.list <- lapply(list, function(x) {rbind(x, matrix(data=NA, nrow=max.row - nrow(x), ncol=ncol(x)))})
    my.mat <- corrected.list[[1]]
    lcl <- length(corrected.list)
    if (lcl > 1) {
      for(i in 2:lcl){
        my.mat <- cbind(my.mat, corrected.list[[i]])
      }
    }
    if (missing(trunc)) {
      trunc <- max.row
    }
    my.mat <- my.mat[1:trunc,,drop=FALSE]
  } else {
    my.mat <- matrix(data=NA, nrow=0, ncol=0)
  }
  return(my.mat)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   is.empty(x)
################
#
################
# Description   :
################
#                   Internal function to represent the empty array, matrix, or vector of zero dimension or length.
#                   Often returned by expressions and functions whose value is undefined.
#                   It returns a logical: TRUE if its argument is empty and FALSE otherwise.
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

is.empty <- function(x) {
  if (is.vector(x)) {
    if((length(x) == 0) || (x == "")) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else if (is.matrix(x) || is.data.frame(x)) {
    return( ((nrow(x) == 0) || (ncol(x) == 0)) )
  } else {
    return( ((length(x) == 0) || (x == "")) )
  }
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   myround (x, digits = 0)
################
#
################
# Description   :
################
#                   Go to the nearest digit rounding
#                   Note that for rounding off a .5, the "go to the even digit" standard is NOT used here.
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

myround <- function (x, digits = 0) {
    upround <- function (x, digits = 0) {
        return(ceiling(x*10^(digits))/10^digits)
    }
    dnround <- function (x, digits = 0) {
        return(floor(x*10^digits)/10^digits)
    }
    i <- (x - trunc(x) >= 0.5)
    x[i] <- upround(x[i], digits=digits)
    x[!i] <- dnround(x[!i], digits=digits)
    return(x)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   .onAttach (libname, pkgname)
################
#
################
# Description   :
################
#                   Startup initializations
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

.onAttach <- function(libname, pkgname) {
    SSver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), fields="Version")
    packageStartupMessage(paste(pkgname, SSver))
    packageStartupMessage("Type PRIMsrc.news() to see new features, changes, and bug fixes\n")
}
##########################################################################################################################################
