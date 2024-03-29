#===============================================================================================================================#
# PRIMsrc
#===============================================================================================================================#

#===============================================================================================================================#
# SURVIVAL INTERNAL SUBROUTINES
# (never to be called by end-user)
#===============================================================================================================================#

#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                   cv.presel(X,
#                             y,
#                             delta,
#                             B,
#                             K,
#                             vs,
#                             vstype,
#                             vsarg,
#                             vscons,
#                             cv,
#                             groups,
#                             onese,
#                             parallel.vs,
#                             parallel.rep,
#                             conf,
#                             verbose,
#                             seed)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

cv.presel <- function(X,
                      y,
                      delta,
                      B,
                      K,
                      vs,
                      vstype,
                      vsarg,
                      vscons,
                      cv,
                      groups,
                      onese,
                      parallel.vs,
                      parallel.rep,
                      conf,
                      verbose,
                      seed) {
  
  # Treatment of variable screening
  if (vs) {
    
    cat("Screening of informative covariates ... \n", sep="")
    
    # Parsing and evaluating all 'vsarg' parameters
    if (vstype == "prsp") {
      alpha <- NULL
      beta <- NULL
      msize <- NULL
      peelcriterion <- NULL
      cvcriterion <- NULL
    } else if (vstype == "pcqr") {
      tau <- NULL
      alpha <- NULL
      nalpha <- NULL
      nlambda <- NULL
    } else if (vstype == "ppl") {
      alpha <- NULL
      nalpha <- NULL
      nlambda <- NULL
    } else if (vstype == "spca") {
      n.thres <- NULL
      n.pcs <- NULL
      n.var <- NULL
    } else {
      stop("Invalid variable screening type option. Exiting ... \n\n")
    }
    eval(parse(text = unlist(strsplit(x = vsarg, split = ","))))
    
    replic <- 1:B
    
    if (!parallel.rep) {
      
      if (!is.null(seed)) {
        seed <- (0:(B-1)) + seed
      }
      
      CV.type.presel.obj <- cv.type.presel(X=X,
                                           y=y,
                                           delta=delta,
                                           replic=replic,
                                           K=K,
                                           vstype=vstype,
                                           vsarg=vsarg,
                                           vscons=vscons,
                                           cv=cv,
                                           groups=groups,
                                           onese=onese,
                                           parallel.vs=parallel.vs,
                                           parallel.rep=parallel.rep,
                                           conf=conf,
                                           verbose=verbose,
                                           seed=seed)
      
    } else {
      
      # Parallel back-end registration
      if (conf$type == "SOCKET") {
        cluster <- parallel::makeCluster(spec=conf$spec,
                                         type="PSOCK",
                                         homogeneous=conf$homo,
                                         outfile=conf$outfile,
                                         verbose=conf$verbose)
      } else if (conf$type == "MPI") {
        cluster <- parallel::makeCluster(spec=conf$spec,
                                         type="MPI",
                                         homogeneous=conf$homo,
                                         outfile=conf$outfile,
                                         verbose=conf$verbose)
      } else {
        stop("Unrecognized cluster type: you must specify a \"SOCKET\" or \"MPI\" cluster type. Exiting ... \n\n")
      }
      
      # Parallelization
      parallel::clusterSetRNGStream(cl=cluster, iseed=seed)
      obj.cl <- parallel::clusterApply(cl=cluster, x=replic, fun=cv.type.presel,
                                       X=X,
                                       y=y,
                                       delta=delta,
                                       K=K,
                                       vstype=vstype,
                                       vsarg=vsarg,
                                       vscons=vscons,
                                       cv=cv,
                                       groups=groups,
                                       onese=onese,
                                       parallel.vs=parallel.vs,
                                       parallel.rep=parallel.rep,
                                       conf=conf,
                                       verbose=verbose,
                                       seed=NULL)
      
      # Stopping the cluster and cleaning all MPI states
      parallel::stopCluster(cl=cluster)
      
      # Collect the results from all replicates
      CV.type.presel.obj <- list("varsel"=vector(mode="list", length=B),
                                 "varsign"=vector(mode="list", length=B),
                                 "boxstat.mean"=vector(mode="list", length=B),
                                 "success"=logical(B))
      for (b in 1:B) {
        CV.type.presel.obj$varsel[[b]] <- obj.cl[[b]]$varsel
        CV.type.presel.obj$varsign[[b]] <- obj.cl[[b]]$varsign
        CV.type.presel.obj$boxstat.mean[[b]] <- obj.cl[[b]]$boxstat.mean
        CV.type.presel.obj$success[b] <- obj.cl[[b]]$success
      }
      CV.type.presel.obj$vstype <- obj.cl[[1]]$vstype
      CV.type.presel.obj$vsarg <- obj.cl[[1]]$vsarg
      CV.type.presel.obj$vscons <- obj.cl[[1]]$vscons
      
    }
    
    vstype <- CV.type.presel.obj$vstype
    vsarg <- CV.type.presel.obj$vsarg
    vscons <- CV.type.presel.obj$vscons
    varsel.list <- CV.type.presel.obj$varsel
    varsign.list <- CV.type.presel.obj$varsign
    boxstat.profile.list <- CV.type.presel.obj$boxstat.mean
    success <-  any(CV.type.presel.obj$success)
    
    # Parsing 'vsarg' to retrieve 'cvcriterion', 'M', and 'msize'
    if (vstype == "prsp") {
      cvcriterion <- vsarg$cvcriterion
    }
    msize <- vsarg$msize
    M <- length(msize)
    
    # Get the fitted model
    if (!success) {
      # Failed to select any covariates
      varsel <- NA
      varsign <- NA
      if (vstype == "prsp") {
        boxstat.profile <- matrix(data=NA, nrow=B, ncol=M)
        colnames(boxstat.profile) <- paste("msize=", msize, sep="")
        boxstat.mean <- rep(NA, M)
        boxstat.se <- rep(NA, M)
        names(boxstat.mean) <- paste("msize=", msize, sep="")
        names(boxstat.se) <- paste("msize=", msize, sep="")
      } else {
        boxstat.profile <- NA
        boxstat.mean <- NA
        boxstat.se <- NA
      }
      boxstat.opt <- NA
      boxstat.1se <- NA
      success <- FALSE
    } else {
      # Get the averaged CV profiles of screening criterion for each PRSP model from all replicates
      if (vstype == "prsp") {
        # Get the averaged CV profiles
        boxstat.profile <- list2mat(list=boxstat.profile.list, fill=NA, coltrunc=M)
        colnames(boxstat.profile) <- paste("msize=", msize, sep="")
        boxstat.mean <- rep(NA, M)
        boxstat.se <- rep(NA, M)
        names(boxstat.mean) <- paste("msize=", msize, sep="")
        names(boxstat.se) <- paste("msize=", msize, sep="")
        if (cvcriterion == "lhr") {
          for (m in 1:M) {
            sumlhr <- rep(x=NA, times=B)
            for (b in 1:B) {
              sumlhr[b] <- boxstat.profile.list[[b]][[m]]
            }
            boxstat.mean[m] <- mean(sumlhr, na.rm=TRUE)
            boxstat.se[m] <- sd(sumlhr, na.rm=TRUE)
          }
          if (all(is.na(boxstat.mean)) || is.empty(boxstat.mean)) {
            boxstat.opt <- NA
            boxstat.1se <- NA
          } else {
            boxstat.opt <- which.max(boxstat.mean)
            w <- boxstat.mean >= boxstat.mean[boxstat.opt]-boxstat.se[boxstat.opt]
            if (all(is.na(w)) || is.empty(w)) {
              boxstat.1se <- NA
            } else {
              boxstat.1se <- min(which(w))
            }
          }
        } else if (cvcriterion == "lrt") {
          for (m in 1:M) {
            sumlrt <- rep(x=NA, times=B)
            for (b in 1:B) {
              sumlrt[b] <- boxstat.profile.list[[b]][[m]]
            }
            boxstat.mean[m] <- mean(sumlrt, na.rm=TRUE)
            boxstat.se[m] <- sd(sumlrt, na.rm=TRUE)
          }
          if (all(is.na(boxstat.mean)) || is.empty(boxstat.mean)) {
            boxstat.opt <- NA
            boxstat.1se <- NA
          } else {
            boxstat.opt <- which.max(boxstat.mean)
            w <- boxstat.mean >= boxstat.mean[boxstat.opt]-boxstat.se[boxstat.opt]
            if (all(is.na(w)) || is.empty(w)) {
              boxstat.1se <- NA
            } else {
              boxstat.1se <- min(which(w))
            }
          }
        } else if (cvcriterion == "cer") {
          for (m in 1:M) {
            sumcer <- rep(x=NA, times=B)
            for (b in 1:B) {
              sumcer[b] <- boxstat.profile.list[[b]][[m]]
            }
            boxstat.mean[m] <- mean(sumcer, na.rm=TRUE)
            boxstat.se[m] <- sd(sumcer, na.rm=TRUE)
          }
          if (all(is.na(boxstat.mean)) || is.empty(boxstat.mean)) {
            boxstat.opt <- NA
            boxstat.1se <- NA
          } else {
            boxstat.opt <- which.min(boxstat.mean)
            w <- boxstat.mean <= boxstat.mean[boxstat.opt]+boxstat.se[boxstat.opt]
            if (all(is.na(w)) || is.empty(w)) {
              boxstat.1se <- NA
            } else {
              boxstat.1se <- min(which(w))
            }
          }
        } else {
          stop("Invalid CV criterion option. Exiting ... \n\n")
        }
      } else {
        boxstat.profile <- NA
        boxstat.mean <- NA
        boxstat.se <- NA
        boxstat.opt <- NA
        boxstat.1se <- NA
      }
      # Get the selected covariates with their signs from all replicates for the selected model
      sel.l <- unlist(varsel.list)
      sign.l <- unlist(varsign.list)
      na <- (is.na(sel.l))
      nna <- length(which(na))
      sel.l <- sel.l[!na]
      sign.l <- sign.l[!na]
      if ((is.empty(sel.l)) || (is.empty(sign.l)) || (all(sign.l == 0)))  {
        varsel <- NA
        varsign <- NA
      } else {
        vartab <- table(names(sel.l), useNA="no")
        nvotes <- B - nna
        w <- names(which(vartab >= ceiling(nvotes*vscons)))
        if (all((is.na(w))) || (is.empty(w))) {
          varsel <- NA
          varsign <- NA
        } else{
          varsel <- pmatch(x=w, table=colnames(X), nomatch = NA_integer_, duplicates.ok = FALSE)
          names(varsel) <- w
          signtab <- table(names(sign.l), sign.l, useNA="no")
          if (length(unique(sign.l)) == 1) {
            signfreq <- rep(x=unique(sign.l), times=nrow(signtab))
            names(signfreq) <- rownames(signtab)
            varsign <- signfreq[w]
          } else if (all(sign.l != 1)) {
            signfreq <- apply(signtab[,c("-1","0"),drop=FALSE], 1, which.max)
            signfreq[signfreq==1] <- -1
            signfreq[signfreq==2] <- 0
            varsign <- signfreq[w]
          } else if (all(sign.l != -1)) {
            signfreq <- apply(signtab[,c("0","1"),drop=FALSE], 1, which.max)
            signfreq[signfreq==1] <- 0
            signfreq[signfreq==2] <- 1
            varsign <- signfreq[w]
          } else {
            signfreq <- apply(signtab[,c("-1","1"),drop=FALSE], 1, which.max)
            signfreq[signfreq==1] <- -1
            signfreq[signfreq==2] <- 1
            varsign <- signfreq[w]
          }
          ord <- order(as.numeric(varsel))
          varsel <- varsel[ord]
          varsign <- varsign[ord]
        }
      }
      if (all(is.na(varsel))) {
        varsel <- NA
        varsign <- NA
        if (vstype == "prsp") {
          boxstat.profile <- matrix(data=NA, nrow=B, ncol=M)
          colnames(boxstat.profile) <- paste("msize=", msize, sep="")
          boxstat.mean <- rep(NA, M)
          boxstat.se <- rep(NA, M)
          names(boxstat.mean) <- paste("msize=", msize, sep="")
          names(boxstat.se) <- paste("msize=", msize, sep="")
        } else {
          boxstat.profile <- NA
          boxstat.mean <- NA
          boxstat.se <- NA
        }
        boxstat.opt <- NA
        boxstat.1se <- NA
        success <- FALSE
      } else {
        success <- TRUE
      }
    }
    
  } else {
    
    cat("No screening of covariates. \n")
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    n <- nrow(X)
    p <- ncol(X)
    
    # Maximum Model Size
    msize <- p
    M <- length(msize)
    
    # Returning argument `vsarg` of variable selection parameters
    vsarg <- list(NA)
    
    # Determination of directions of peeling of all covariates
    if (cv) {
      cv.fit <- glmnet::cv.glmnet(x=X,
                                  y=survival::Surv(time=y, event=delta),
                                  alpha=0,
                                  nlambda=100,
                                  nfolds=pmax(3,K),
                                  family="cox",
                                  maxit=1e+05)
      fit <- glmnet::glmnet(x=X,
                            y=survival::Surv(time=y, event=delta),
                            alpha=0,
                            family="cox",
                            maxit=1e+05)
      if (onese) {
        cv.coef <- as.numeric(coef(fit, s=cv.fit$lambda.1se))
      } else {
        cv.coef <- as.numeric(coef(fit, s=cv.fit$lambda.min))
      }
    } else {
      fit <- glmnet::glmnet(x=X,
                            y=survival::Surv(time=y, event=delta),
                            alpha=0,
                            family="cox",
                            maxit=1e+05)
      cv.coef <- as.numeric(coef(fit, s=fit$lambda[50]))
    }
    w <- which(!(is.na(cv.coef)) & (cv.coef != 0))
    if (is.empty(w)) {
      varsel <- NA
      varsign <- NA
      success <- FALSE
    } else {
      varsel <- (1:ncol(X))[w]
      varsign <- sign(cv.coef[w])
      names(varsel) <- colnames(X)[w]
      names(varsign) <- colnames(X)[w]
      success <- TRUE      
    }
    boxstat.profile <- matrix(data=NA, nrow=B, ncol=M)
    boxstat.mean <- rep(NA, M)
    boxstat.se <- rep(NA, M)
    boxstat.opt <- NA
    boxstat.1se <- NA    
  }
  
  return(list("vstype"=vstype,
              "vsarg"=vsarg,
              "vscons"=vscons,
              "varsel"=varsel,
              "varsign"=varsign,
              "boxstat.profile"=boxstat.profile,
              "boxstat.mean"=boxstat.mean,
              "boxstat.se"=boxstat.se,
              "boxstat.opt"=boxstat.opt,
              "boxstat.1se"=boxstat.1se,
              "success"=success,
              "seed"=seed))
  
}
#===============================================================================================================================#





#===============================================================================================================================#
#
#===============================================================================================================================#

cv.type.presel <- function(X,
                           y,
                           delta,
                           replic,
                           K,
                           vstype,
                           vsarg,
                           vscons,
                           cv,
                           groups,
                           onese,
                           parallel.vs,
                           parallel.rep,
                           conf,
                           verbose,
                           seed) {
  
  if (!parallel.rep) {
    
    B <- length(replic)
    success <- logical(B)
    varsel <- vector(mode="list", length=B)
    varsign <- vector(mode="list", length=B)
    boxstat.mean <- vector(mode="list", length=B)
    b <- 1
    while (b <= B) {
      cat("Replicate : ", b, "\n", sep="")
      cat("seed : ", seed[b], "\n", sep="")
      if (vstype == "prsp") {
        CVSEL <- cv.prsp(X=X,
                         y=y,
                         delta=delta,
                         K=K,
                         vsarg=vsarg,
                         vscons=vscons,
                         cv=cv,
                         groups=groups,
                         onese=onese,
                         parallel=parallel.vs,
                         conf=conf,
                         verbose=verbose,
                         seed=seed[b])
        boxstat.mean[[b]] <- CVSEL$boxstat.mean
      } else if (vstype == "pcqr") {
        CVSEL <- cv.pcqr(X=X,
                         y=y,
                         delta=delta,
                         K=K,
                         vsarg=vsarg,
                         vscons=vscons,
                         cv=cv,
                         onese=onese,
                         parallel=parallel.vs,
                         conf=conf,
                         verbose=verbose,
                         seed=seed[b])
        boxstat.mean[[b]] <- NULL
      } else if (vstype == "ppl") {
        CVSEL <- cv.ppl(X=X,
                        y=y,
                        delta=delta,
                        K=K,
                        vsarg=vsarg,
                        vscons=vscons,
                        cv=cv,
                        onese=onese,
                        parallel=parallel.vs,
                        conf=conf,
                        verbose=verbose,
                        seed=seed[b])
        boxstat.mean[[b]] <- NULL
      } else if (vstype == "spca") {
        CVSEL <- cv.spca(X=X,
                         y=y,
                         delta=delta,
                         K=K,
                         vsarg=vsarg,
                         vscons=vscons,
                         cv=cv,
                         parallel=parallel.vs,
                         conf=conf,
                         verbose=verbose,
                         seed=seed[b])
        boxstat.mean[[b]] <- NULL
      } else {
        stop("Invalid variable screening type option. Exiting ... \n\n")
      }
      varsel[[b]] <- CVSEL$varsel
      varsign[[b]] <- CVSEL$varsign
      success[b] <- CVSEL$success
      b <- b + 1
    }
    
  } else {
    
    if (vstype == "prsp") {
      CVSEL <- cv.prsp(X=X,
                       y=y,
                       delta=delta,
                       K=K,
                       vsarg=vsarg,
                       vscons=vscons,
                       cv=cv,
                       groups=groups,
                       onese=onese,
                       parallel=parallel.vs,
                       conf=conf,
                       verbose=verbose,
                       seed=seed)
      boxstat.mean <- CVSEL$boxstat.mean
    } else if (vstype == "pcqr") {
      CVSEL <- cv.pcqr(X=X,
                       y=y,
                       delta=delta,
                       K=K,
                       vsarg=vsarg,
                       vscons=vscons,
                       cv=cv,
                       onese=onese,
                       parallel=parallel.vs,
                       conf=conf,
                       verbose=verbose,
                       seed=seed)
      boxstat.mean <- NULL
    } else if (vstype == "ppl") {
      CVSEL <- cv.ppl(X=X,
                      y=y,
                      delta=delta,
                      K=K,
                      vsarg=vsarg,
                      vscons=vscons,
                      cv=cv,
                      onese=onese,
                      parallel=parallel.vs,
                      conf=conf,
                      verbose=verbose,
                      seed=seed)
      boxstat.mean <- NULL
    } else if (vstype == "spca") {
      CVSEL <- cv.spca(X=X,
                       y=y,
                       delta=delta,
                       K=K,
                       vsarg=vsarg,
                       vscons=vscons,
                       cv=cv,
                       parallel=parallel.vs,
                       conf=conf,
                       verbose=verbose,
                       seed=seed)
      boxstat.mean <- NULL
    } else {
      stop("Invalid variable screening type option. Exiting ... \n\n")
    }
    varsel <- CVSEL$varsel
    varsign <- CVSEL$varsign
    success <- CVSEL$success
    
  }
  
  return(list("vstype"=vstype,
              "vsarg"=CVSEL$vsarg,
              "vscons"=vscons,
              "varsel"=varsel,
              "varsign"=varsign,
              "boxstat.mean"=boxstat.mean,
              "success"=success,
              "seed"=seed))
}
#===============================================================================================================================#




#===============================================================================================================================#
# Variables screening by univariate PRSP screening (PRSP)
# CV over number of ordered statistics used to compute variable frequencies
#===============================================================================================================================#

cv.prsp <- function(X,
                    y,
                    delta,
                    K,
                    vsarg,
                    vscons,
                    cv,
                    groups,
                    onese,
                    parallel,
                    conf,
                    verbose,
                    seed) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Parsing and evaluating 'vsarg' argument to evaluate all parameters
  alpha <- NULL
  beta <- NULL
  msize <- NULL
  peelcriterion <- NULL
  cvcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=vsarg, split=",")) ))
  
  if (is.null(msize)) {
    cv <- TRUE
  } else {
    cv <- FALSE
  }
  
  # Maximal possible number of peeling steps
  Lmax <- ceiling(log(beta) / log(1 - alpha))
  
  # Maximum Model Size or Vector of Model Sizes
  if (!is.null(msize)) {
    Smax <- max(1,floor(p))
    if (msize < 1) {
      cat("Warning: Parameter `msize` was less than the minimal allowed value and was reset to ", 1, ".\n", sep="")
      msize <- max(1,floor(msize))
    }
    if (msize > Smax) {
      cat("Warning: Parameter `msize` was greater than the maximal allowed vale and was reset to ", Smax, ".\n", sep="")
      msize <- min(Smax, msize)
    }
  } else {
    Smax <- max(1,floor(p/5))
    msize <- unique(ceiling(seq(from=max(1,floor(Smax/100)), to=Smax, length=floor(100*Smax/p))))
    cat("Model sizes to be explored by cross-validation: ", msize, "\n", sep=" ")
  }
  M <- length(msize)
  
  # Creating argument `arg` of `cv.prsp.tune` for variable selection parameters
  arg <- paste("alpha=", alpha,
               ",beta=", beta,
               ",L=", Lmax,
               ",peelcriterion=\"", peelcriterion, "\"",
               ",cvcriterion=\"", cvcriterion, "\"", sep="")
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Parallel back-end registration
  if (parallel) {
    if (conf$type == "SOCKET") {
      cluster <- parallel::makeCluster(spec=conf$spec,
                                       type="PSOCK",
                                       homogeneous=conf$homo,
                                       outfile=file.path(getwd(), "output_SOCK_inner.txt", fsep=.Platform$file.sep),
                                       verbose=conf$verbose)
    } else if (conf$type == "MPI") {
      cluster <- parallel::makeCluster(spec=conf$spec,
                                       type="MPI",
                                       homogeneous=conf$homo,
                                       outfile=file.path(getwd(), "output_MPI_inner.txt", fsep=.Platform$file.sep),
                                       verbose=conf$verbose)
    }
  }
  
  folds <- cv.folds(y=delta, K=K, seed=seed)
  varsel.list <- vector(mode="list", length=K)
  varsign.list <- vector(mode="list", length=K)
  boxstat.list <- vector(mode="list", length=K)
  k <- 1
  while (k <= K) {
    if ((verbose) && (cv)) cat("Fold : ", k, "\n", sep="")
    # Initialize training and test data
    if (K == 1) {
      traindata <- testdata <- X[folds$perm[(folds$which == k)], , drop=FALSE]
      traintime <- testtime <- y[folds$perm[(folds$which == k)]]
      trainstatus <- teststatus <- delta[folds$perm[(folds$which == k)]]
      traingroups <- testgroups <- groups[folds$perm[(folds$which == k)]]
    } else {
      traindata <- X[folds$perm[(folds$which != k)], , drop=FALSE]
      traintime <- y[folds$perm[(folds$which != k)]]
      trainstatus <- delta[folds$perm[(folds$which != k)]]
      traingroups <- groups[folds$perm[(folds$which != k)]]
      testdata <- X[folds$perm[(folds$which == k)], , drop=FALSE]
      testtime <- y[folds$perm[(folds$which == k)]]
      teststatus <- delta[folds$perm[(folds$which == k)]]
      testgroups <- groups[folds$perm[(folds$which == k)]]
    }
    # Training the model for each length of ranked statistics
    if (!parallel) {
      cv.obj <- cv.prsp.tune(traindata=traindata,
                             traintime=traintime,
                             trainstatus=trainstatus,
                             traingroups=traingroups,
                             testdata=testdata,
                             testtime=testtime,
                             teststatus=teststatus,
                             testgroups=testgroups,
                             arg=arg,
                             msize=msize,
                             parallel=parallel,
                             verbose=verbose)
    } else {
      parallel::clusterSetRNGStream(cl=cluster, iseed=seed)
      obj.cl <- parallel::clusterApply(cl=cluster, x=msize, fun=cv.prsp.tune,
                                       traindata=traindata,
                                       traintime=traintime,
                                       trainstatus=trainstatus,
                                       traingroups=traingroups,
                                       testdata=testdata,
                                       testtime=testtime,
                                       teststatus=teststatus,
                                       testgroups=testgroups,
                                       arg=arg,
                                       parallel=parallel,
                                       verbose=verbose)
      cv.obj <- list("varsel"=vector(mode="list", length=M),
                     "varsign"=vector(mode="list", length=M),
                     "boxstat"=vector(mode="list", length=M),
                     "success"=TRUE)
      for (m in 1:M) {
        cv.obj$varsel[[m]] <- obj.cl[[m]]$varsel
        cv.obj$varsign[[m]] <- obj.cl[[m]]$varsign
        cv.obj$boxstat[[m]] <- obj.cl[[m]]$boxstat
        cv.obj$success <- (cv.obj$success & obj.cl[[m]]$success)
      }
    }
    if (cv.obj$success) {
      # Store the test set results for each model from each fold
      boxstat.list[[k]] <- cv.obj$boxstat
      varsel.list[[k]] <- cv.obj$varsel
      varsign.list[[k]] <- cv.obj$varsign
    } else {
      boxstat.list[[k]] <- as.list(rep(NA, M))
      varsel.list[[k]] <- as.list(rep(NA, M))
      varsign.list[[k]] <- as.list(rep(NA, M))
    }
    k <- k + 1
  }
  
  # Stopping the cluster and cleaning all MPI states
  if (parallel) {
    parallel::stopCluster(cl=cluster)
  }
  
  # Get the fitted model
  if (all(is.na(varsel.list))) {
    # Failed to select any covariates
    varsel <- NA
    varsign <- NA
    boxstat.mean <- rep(NA, M)
    boxstat.opt <- NA
    boxstat.1se <- NA
    msize <- NA
    success <- FALSE
  } else {
    # Get the averaged CV profiles of screening criterion for each model from all folds
    boxstat.mean <- rep(NA, M)
    boxstat.se <- rep(NA, M)
    if (cvcriterion == "lhr") {
      for (m in 1:M) {
        sumlhr <- rep(x=NA, times=K)
        for (k in 1:K) {
          sumlhr[k] <- boxstat.list[[k]][[m]]
        }
        boxstat.mean[m] <- mean(sumlhr, na.rm=TRUE)
        boxstat.se[m] <- sd(sumlhr, na.rm=TRUE) / sqrt(n/K)
      }
      if (all(is.na(boxstat.mean)) || is.empty(boxstat.mean)) {
        boxstat.opt <- NA
        boxstat.1se <- NA
      } else {
        boxstat.opt <- which.max(boxstat.mean)
        w <- boxstat.mean >= boxstat.mean[boxstat.opt]-boxstat.se[boxstat.opt]
        if (all(is.na(w)) || is.empty(w)) {
          boxstat.1se <- NA
        } else {
          boxstat.1se <- min(which(w))
        }
      }
    } else if (cvcriterion == "lrt") {
      for (m in 1:M) {
        sumlrt <- rep(x=NA, times=K)
        for (k in 1:K) {
          sumlrt[k] <- boxstat.list[[k]][[m]]
        }
        boxstat.mean[m] <- mean(sumlrt, na.rm=TRUE)
        boxstat.se[m] <- sd(sumlrt, na.rm=TRUE) / sqrt(n/K)
      }
      if (all(is.na(boxstat.mean)) || is.empty(boxstat.mean)) {
        boxstat.opt <- NA
        boxstat.1se <- NA
      } else {
        boxstat.opt <- which.max(boxstat.mean)
        w <- boxstat.mean >= boxstat.mean[boxstat.opt]-boxstat.se[boxstat.opt]
        if (all(is.na(w)) || is.empty(w)) {
          boxstat.1se <- NA
        } else {
          boxstat.1se <- min(which(w))
        }
      }
    } else if (cvcriterion == "cer") {
      for (m in 1:M) {
        sumcer <- rep(x=NA, times=K)
        for (k in 1:K) {
          sumcer[k] <- boxstat.list[[k]][[m]]
        }
        boxstat.mean[m] <- mean(sumcer, na.rm=TRUE)
        boxstat.se[m] <- sd(sumcer, na.rm=TRUE) / sqrt(n/K)
      }
      if (all(is.na(boxstat.mean)) || is.empty(boxstat.mean)) {
        boxstat.opt <- NA
        boxstat.1se <- NA
      } else {
        boxstat.opt <- which.min(boxstat.mean)
        w <- boxstat.mean <= boxstat.mean[boxstat.opt]+boxstat.se[boxstat.opt]
        if (all(is.na(w)) || is.empty(w)) {
          boxstat.1se <- NA
        } else {
          boxstat.1se <- min(which(w))
        }
      }
    } else {
      stop("Invalid CV criterion option. Exiting ... \n\n")
    }
    # Get the selected covariates with their signs from all folds for each model
    varsel <- vector(mode="list", length=M)
    varsign <- vector(mode="list", length=M)
    for (m in 1:M) {
      sel.m <- numeric(0)
      sign.m <- numeric(0)
      for (k in 1:K) {
        sel.m <- c(sel.m, varsel.list[[k]][[m]])
        sign.m <- c(sign.m, varsign.list[[k]][[m]])
      }
      na <- (is.na(sel.m))
      nna <- sum(na)
      sel.m <- sel.m[!na]
      sign.m <- sign.m[!na]
      if ((is.empty(sel.m)) || (is.empty(sign.m)) || (all(sign.m == 0)))  {
        varsel[[m]] <- NA
        varsign[[m]] <- NA
      } else {
        vartab <- table(names(sel.m), useNA="no")
        nvotes <- K - nna
        w <- names(which(vartab >= ceiling(nvotes*vscons)))
        if ((all(is.na(w))) || (is.empty(w))) {
          varsel[[m]] <- NA
          varsign[[m]] <- NA
        } else{
          varsel[[m]] <- pmatch(x=w, table=colnames(X), nomatch = NA_integer_, duplicates.ok = FALSE)
          names(varsel[[m]]) <- w
          signtab <- table(names(sign.m), sign.m, useNA="no")
          if (length(unique(sign.m)) == 1) {
            signfreq <- rep(x=unique(sign.m), times=nrow(signtab))
            names(signfreq) <- rownames(signtab)
            varsign[[m]] <- signfreq[w]
          } else if (all(sign.m != 1)) {
            signfreq <- apply(signtab[,c("-1","0"),drop=FALSE], 1, which.max)
            signfreq[signfreq==1] <- -1
            signfreq[signfreq==2] <- 0
            varsign[[m]] <- signfreq[w]
          } else if (all(sign.m != -1)) {
            signfreq <- apply(signtab[,c("0","1"),drop=FALSE], 1, which.max)
            signfreq[signfreq==1] <- 0
            signfreq[signfreq==2] <- 1
            varsign[[m]] <- signfreq[w]
          } else {
            signfreq <- apply(signtab[,c("-1","1"),drop=FALSE], 1, which.max)
            signfreq[signfreq==1] <- -1
            signfreq[signfreq==2] <- 1
            varsign[[m]] <- signfreq[w]
          }
          ord <- order(as.numeric(varsel[[m]]))
          varsel[[m]] <- varsel[[m]][ord]
          varsign[[m]] <- varsign[[m]][ord]
        }
      }
    }
    if (all(is.na(varsel))) {
      varsel <- NA
      varsign <- NA
      boxstat.mean <- rep(NA, M)
      boxstat.opt <- NA
      boxstat.1se <- NA
      success <- FALSE
      msize <- NA
    } else {
      if (onese) {
        varsel <- varsel[[boxstat.1se]]
        varsign <- varsign[[boxstat.1se]]
      } else {
        varsel <- varsel[[boxstat.opt]]
        varsign <- varsign[[boxstat.opt]]
      }
      success <- TRUE
    }
  }
  
  # Returning updated argument `vsarg` of variable selection parameters
  vsarg <- list("alpha"=alpha,
                "beta"=beta,
                "peelcriterion"=peelcriterion,
                "cvcriterion"=cvcriterion,
                "msize"=msize)
  
  return(list("vsarg"=vsarg,
              "vscons"=vscons,
              "varsel"=varsel,
              "varsign"=varsign,
              "boxstat.mean"=boxstat.mean,
              "boxstat.opt"=boxstat.opt,
              "boxstat.1se"=boxstat.1se,
              "success"=success,
              "seed"=seed))
}
#===============================================================================================================================#




#===============================================================================================================================#
#
#===============================================================================================================================#

cv.prsp.tune <- function(traindata,
                         traintime,
                         trainstatus,
                         traingroups,
                         testdata,
                         testtime,
                         teststatus,
                         testgroups,
                         arg,
                         msize,
                         parallel,
                         verbose) {
  
  # Parsing and evaluating 'arg' parameters to evaluate 'cvcriterion'
  alpha <- NULL
  beta <- NULL
  L <- NULL
  peelcriterion <- NULL
  cvcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=arg, split=",")) ))
  
  if (!parallel) {
    
    M <- length(msize)
    varsel <- vector(mode="list", length=M)
    varsign <- vector(mode="list", length=M)
    boxstat <- vector(mode="list", length=M)
    peelobj <- cv.prsp.univ(traindata=traindata,
                            traintime=traintime,
                            trainstatus=trainstatus,
                            traingroups=traingroups,
                            arg=arg)
    m <- 1
    while (m <= M) {
      if (verbose) cat("Model Size: ", msize[m], "\n")
      if (is.empty(peelobj$varsel)) {
        success <- FALSE
        varsel <- NULL
        varsign <- NULL
        boxstat <- NULL
        m <- M + 1
      } else {
        success <- TRUE
        # Retrieving screened variables and box from the univariate screening
        varsel[[m]] <- peelobj$varsel[1:msize[m]]
        varsign[[m]] <- peelobj$varsign[1:msize[m]]
        # Extract the rule and sign as one vector
        varcut <- peelobj$varcut[1:msize[m]]
        test.cut <- t(t(testdata[,varsel[[m]], drop=FALSE]) * varsign[[m]])
        test.ind <- t(t(test.cut) >= varcut * varsign[[m]])
        # Select which test observations are `TRUE` for all screened covariates
        pred <- (rowMeans(test.ind) == 1)
        test.ind1 <- 1*pred
        if ((sum(pred, na.rm=TRUE) != length(pred[!is.na(pred)])) && (sum(pred, na.rm=TRUE) != 0)) {
          surv.formula <- (survival::Surv(testtime, teststatus) ~ 1 + test.ind1)
          if (cvcriterion == "lhr") {
            coxobj <- survival::coxph(surv.formula, na.action=na.exclude, iter.max=20, timefix=TRUE, method="efron")
            boxstat[[m]] <- coxobj$coef
          } else if (cvcriterion == "lrt") {
            boxstat[[m]] <- survival::survdiff(surv.formula, rho=0)$chisq
          } else if (cvcriterion == "cer") {
            coxobj <- survival::coxph(surv.formula, na.action=na.exclude, iter.max=20, timefix=TRUE, method="efron")
            predobj <- predict(object=coxobj, type="lp", reference="sample", na.action=na.exclude)
            boxstat[[m]] <- Hmisc::rcorr.cens(x=predobj, S=survival::Surv(testtime, teststatus))['C Index']
          } else {
            stop("Invalid CV criterion option. Exiting ... \n\n")
          }
        } else {
          if (cvcriterion == "lhr") {
            boxstat[[m]] <- 0
          } else if (cvcriterion == "lrt") {
            boxstat[[m]] <- 0
          } else if (cvcriterion == "cer") {
            boxstat[[m]] <- 1
          } else {
            stop("Invalid CV criterion option. Exiting ... \n\n")
          }
        }
        m <- m + 1
      }
    }
    
  } else {
    
    peelobj <- cv.prsp.univ(traindata=traindata,
                            traintime=traintime,
                            trainstatus=trainstatus,
                            traingroups=traingroups,
                            arg=arg)
    if (is.empty(peelobj$varsel)) {
      success <- FALSE
      varsel <- NULL
      varsign <- NULL
      boxstat <- NULL
    } else {
      # Retrieving screened variables and box from the univariate screening
      success <- TRUE
      varsel <- peelobj$varsel[1:msize]
      varsign <- peelobj$varsign[1:msize]
      varcut <- peelobj$varcut[1:msize]
      # Extract the rule and sign as one vector
      test.cut <- t(t(testdata[,varsel, drop=FALSE]) * varsign)
      test.ind <- t(t(test.cut) >= varcut * varsign)
      # Set as TRUE which observations are TRUE for all covariates
      pred <- (rowMeans(test.ind) == 1)
      test.ind1 <- 1*pred
      if ((sum(pred, na.rm=TRUE) != length(pred[!is.na(pred)])) && (sum(pred, na.rm=TRUE) != 0)) {
        surv.formula <- (survival::Surv(testtime, teststatus) ~ 1 + test.ind1)
        if (cvcriterion == "lhr") {
          coxobj <- survival::coxph(surv.formula, na.action=na.exclude, iter.max=20, timefix=TRUE, method="efron")
          boxstat <- coxobj$coef
        } else if (cvcriterion == "lrt") {
          boxstat <- survival::survdiff(surv.formula, rho=0)$chisq
        } else if (cvcriterion == "cer") {
          coxobj <- survival::coxph(surv.formula, na.action=na.exclude, iter.max=20, timefix=TRUE, method="efron")
          predobj <- predict(object=coxobj, type="lp", reference="sample", na.action=na.exclude)
          boxstat <- Hmisc::rcorr.cens(x=predobj, S=survival::Surv(testtime, teststatus))['C Index']
        } else {
          stop("Invalid CV criterion option. Exiting ... \n\n")
        }
      } else {
        if (cvcriterion == "lhr") {
          boxstat <- 0
        } else if (cvcriterion == "lrt") {
          boxstat <- 0
        } else if (cvcriterion == "cer") {
          boxstat <- 1
        } else {
          stop("Invalid CV criterion option. Exiting ... \n\n")
        }
      }
    }
  }
  
  return(list("varsel"=varsel,
              "varsign"=varsign,
              "boxstat"=boxstat,
              "success"=success))
}
#===============================================================================================================================#



#===============================================================================================================================#
#
#===============================================================================================================================#

cv.prsp.univ <- function(traindata,
                         traintime,
                         trainstatus,
                         traingroups,
                         arg) {
  
  # Parsing and evaluating 'arg' to evaluate all parameters
  alpha <- NULL
  beta <- NULL
  L <- NULL
  peelcriterion <- NULL
  cvcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=arg, split=",")) ))
  
  # Creating argument `arg` of `cv.comb.peel` for PRSP parameters
  arg <- paste("alpha=", alpha,
               ",beta=", beta,
               ",L=", L,
               ",peelcriterion=\"", peelcriterion, "\"", sep="")
  
  # Univariate screening
  p <- ncol(traindata)
  varstat <- numeric(p)
  varsel <- 1:p
  varsign <- numeric(p)
  varcut <- numeric(p)
  names(varsel) <- colnames(traindata)
  names(varsign) <- colnames(traindata)
  names(varcut) <- colnames(traindata)
  for (j in 1:p) {
    prsp.fit <- prsp(traindata=traindata[, j, drop=FALSE],
                     traintime=traintime,
                     trainstatus=trainstatus,
                     varsign=NULL,
                     initcutpts=NULL,
                     arg=arg,
                     traingroups=traingroups)
    cat("covariate: ", j, "\t (maxstep: ", prsp.fit$maxsteps, ")\n", sep="")
    varcut[j] <- prsp.fit$boxcut[1,1,drop=TRUE]
    varsign[j] <- prsp.fit$varsign
    Lm <- prsp.fit$maxsteps
    l <- 1
    boxstat <- numeric(Lm)
    while (l <= Lm) {
      # Extract the rule and sign as one vector
      boxcut <- prsp.fit$boxcut[l,1,drop=TRUE]
      train.cut <- t(t(traindata[,j,drop=FALSE]) * varsign[j])
      train.ind <- t(t(train.cut) >= boxcut * varsign[j])
      # Select which test observations are `TRUE`
      pred <- (rowMeans(train.ind) == 1)
      train.ind1 <- 1*pred
      if ((sum(pred, na.rm=TRUE) != length(pred[!is.na(pred)])) && (sum(pred, na.rm=TRUE) != 0)) {
        surv.formula <- (survival::Surv(traintime, trainstatus) ~ 1 + train.ind1)
        if (cvcriterion == "lhr") {
          coxobj <- survival::coxph(surv.formula, na.action=na.exclude, iter.max=20, timefix=TRUE, method="efron")
          boxstat[l] <- coxobj$coef
        } else if (cvcriterion == "lrt") {
          boxstat[l] <- survival::survdiff(surv.formula, rho=0)$chisq
        } else if (cvcriterion == "cer") {
          coxobj <- survival::coxph(surv.formula, na.action=na.exclude, iter.max=20, timefix=TRUE, method="efron")
          predobj <- predict(object=coxobj, type="lp", reference="sample", na.action=na.exclude)
          boxstat[l] <- Hmisc::rcorr.cens(x=predobj, S=survival::Surv(traintime, trainstatus))['C Index']
        } else {
          stop("Invalid CV criterion option. Exiting ... \n\n")
        }
      } else {
        if (cvcriterion == "lhr") {
          boxstat[l] <- 0
        } else if (cvcriterion == "lrt") {
          boxstat[l] <- 0
        } else if (cvcriterion == "cer") {
          boxstat[l] <- 1
        } else {
          stop("Invalid CV criterion option. Exiting ... \n\n")
        }
      }
      l <- l + 1
    }
    if (cvcriterion == "lhr") {
      if (all(is.na(boxstat)) || is.empty(boxstat)) {
        varstat[j] <- NA
      } else {
        boxstat.opt <- which.max(boxstat)
        varstat[j] <- boxstat[boxstat.opt]
      }
    } else if (cvcriterion == "lrt") {
      if (all(is.na(boxstat)) || is.empty(boxstat)) {
        varstat[j] <- NA
      } else {
        boxstat.opt <- which.max(boxstat)
        varstat[j] <- boxstat[boxstat.opt]
      }
    } else if (cvcriterion == "cer") {
      if (all(is.na(boxstat)) || is.empty(boxstat)) {
        varstat[j] <- NA
      } else {
        boxstat.opt <- which.min(boxstat)
        varstat[j] <- boxstat[boxstat.opt]
      }
    } else {
      stop("Invalid CV criterion option. Exiting ... \n\n")
    }
  }
  
  # Order by variable statistics
  if (cvcriterion == "lhr") {
    ord <- order(varstat, decreasing=TRUE, na.last=TRUE)
  } else if (cvcriterion == "lrt") {
    ord <- order(varstat, decreasing=TRUE, na.last=TRUE)
  } else if (cvcriterion == "cer") {
    ord <- order(varstat, decreasing=FALSE, na.last=TRUE)
  } else {
    stop("Invalid CV criterion option. Exiting ... \n\n")
  }
  varsel <- varsel[ord]
  varsign <- varsign[ord]
  varcut <- varcut[ord]
  
  return(list("varsel"=varsel,
              "varsign"=varsign,
              "varcut"=varcut))
}
#===============================================================================================================================#



#===============================================================================================================================#
# Variables screening by Penalized Censored Quantile Regression model (PCQR)
# CV of variable selection using full cross-validation of both parameters alpha and lambda
#===============================================================================================================================#

cv.pcqr <- function(X,
                    y,
                    delta,
                    K,
                    vsarg,
                    vscons,
                    cv,
                    onese,
                    parallel,
                    conf,
                    verbose,
                    seed) {
  
  # Parsing and evaluating 'vsarg' argument to evaluate all parameters
  tau <- NULL
  alpha <- NULL
  nalpha <- NULL
  nlambda <- NULL
  eval(parse( text=unlist(strsplit(x=vsarg, split=",")) ))
  
  M <- 2
  msize <- numeric(M)
  
  if (is.null(alpha)) {
    enalpha <- seq(from=0, to=1, length.out=nalpha)
  } else {
    nalpha <- 1
    enalpha <- alpha
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Parallel back-end registration
  if (parallel) {
    if (conf$type == "SOCKET") {
     cluster <- parallel::makeCluster(spec=conf$spec,
                                       type="PSOCK",
                                       homogeneous=conf$homo,
                                       outfile=file.path(getwd(), "output_SOCK_inner.txt", fsep=.Platform$file.sep),
                                       verbose=conf$verbose)
    } else if (conf$type == "MPI") {
      cluster <- parallel::makeCluster(spec=conf$spec,
                                       type="MPI",
                                       homogeneous=conf$homo,
                                       outfile=file.path(getwd(), "output_MPI_inner.txt", fsep=.Platform$file.sep),
                                       verbose=conf$verbose)
    }
  }
  
  folds <- cv.folds(y=delta, K=K, seed=seed)
  varsel.list <- vector(mode="list", length=K)
  varsign.list <- vector(mode="list", length=K)
  k <- 1
  while (k <= K) {
    if ((verbose) && (cv)) cat("Fold : ", k, "\n", sep="")
    # Initialize training and test data
    if (K == 1) {
      traindata <- testdata <- X[folds$perm[(folds$which == k)], , drop=FALSE]
      traintime <- testtime <- y[folds$perm[(folds$which == k)]]
      trainstatus <- teststatus <- delta[folds$perm[(folds$which == k)]]
    } else {
      traindata <- X[folds$perm[(folds$which != k)], , drop=FALSE]
      traintime <- y[folds$perm[(folds$which != k)]]
      trainstatus <- delta[folds$perm[(folds$which != k)]]
      testdata <- X[folds$perm[(folds$which == k)], , drop=FALSE]
      testtime <- y[folds$perm[(folds$which == k)]]
      teststatus <- delta[folds$perm[(folds$which == k)]]
    }
    enlambda <- vector(mode="list", length=nalpha)
    cv.errmu <- vector(mode="list", length=nalpha)
    cv.errse <- vector(mode="list", length=nalpha)
    for (i in 1:nalpha) {
      cv.fit <- cv.pcqreg(X=traindata,
                          y=traintime,
                          delta=trainstatus,
                          method="cquantile",
                          gamma=IQR(traintime)/10,
                          tau=tau,
                          nfolds=pmax(3,K),
                          type.loss="deviance",
                          alpha=enalpha[i],
                          nlambda=nlambda,
                          lambda.min=ifelse(nrow(X)>ncol(X), 0.001, 0.05),
                          preprocess="standardize",
                          screen="ASR",
                          max.iter=1e4,
                          eps=1e-7,
                          dfmax=ncol(X)+1,
                          penalty.factor=rep(x=1, times=ncol(X)),
                          parallel=parallel,
                          cluster=cluster,
                          verbose=verbose,
                          seed=seed)
      cv.errmu[[i]] <- cv.fit$cve
      cv.errse[[i]] <- cv.fit$cvse
      enlambda[[i]] <- cv.fit$lambda
    }
    cv.errmu <- list2mat(list=cv.errmu, coltrunc="max", fill=NA)
    cv.errse <- list2mat(list=cv.errse, coltrunc="max", fill=NA)
    cv.errmu.index <- as.matrix(which(x=(cv.errmu == as.numeric(cv.errmu)[which.min(cv.errmu)]), arr.ind=TRUE, useNames=FALSE))
    w <- which.min(cv.errmu.index[, 1, drop=TRUE])
    index.min <- as.numeric(cv.errmu.index[w,,drop=TRUE])
    ww <- as.matrix(which(x=cv.errmu <= cv.errmu[index.min[1],index.min[2]] + cv.errse[index.min[1],index.min[2]], arr.ind=TRUE, useNames=TRUE))
    index.1se <- c(min(ww[,1]), min(ww[ww[,1]==min(ww[,1]),2]))
    if (is.empty(index.min) || is.empty(index.1se)) {
      varsel.list[[k]] <- list(NA, NA)
      varsign.list[[k]] <- list(NA, NA)
    } else {
      enalpha.min <- enalpha[index.min[1]]
      enalpha.1se <- enalpha[index.1se[1]]
      enlambda.min <- enlambda[[index.min[1]]][index.min[2]]
      enlambda.1se <- enlambda[[index.1se[1]]][index.1se[2]]
      fit <- pcqreg(X=traindata,
                    y=traintime,
                    delta=trainstatus,
                    method="cquantile",
                    gamma=IQR(traintime)/10,
                    tau=tau,
                    alpha=ifelse(test=onese, yes=enalpha.1se, no=enalpha.min),
                    nlambda=nlambda,
                    lambda.min=ifelse(nrow(X)>ncol(X), 0.001, 0.05),
                    lambda=c(enlambda.1se, enlambda.min),
                    preprocess="standardize",
                    screen="ASR",
                    max.iter=1e4,
                    eps=1e-7,
                    dfmax=ncol(X)+1,
                    penalty.factor=rep(x=1, times=ncol(X)),
                    verbose=verbose)
      w <- apply(X=fit$beta, MARGIN=2, FUN=function(x) {sum(!(is.na(x)) & (x != 0))})
      if (all(w == 0)) {
        varsel.list[[k]] <- list(NA, NA)
        varsign.list[[k]] <- list(NA, NA)
      } else {
        if (length(fit$lambda) <= 2) {
          cv.coef <- -coef.pcqreg(fit, lambda=c(enlambda.1se, enlambda.min), exact=TRUE)[-1,]
          cv.coef <- list(cv.coef[,1], cv.coef[,2])
          varsel <- lapply(cv.coef, FUN=function(x) {which(!(is.na(x)) & (x != 0))})
          if (all(is.na(varsel))) {
            varsel.list[[k]] <- list(NA, NA)
            varsign.list[[k]] <- list(NA, NA)
          } else {
            varsel.list[[k]] <- varsel
            varsign.list[[k]] <- lapply(cv.coef, function(x) sign(x))
          }
        } else {
          cv.coef <- -coef.pcqreg(fit, lambda=c(enlambda.1se, enlambda.min), exact=FALSE)[-1,]
          cv.coef <- list(cv.coef[,1], cv.coef[,2])
          varsel <- lapply(cv.coef, FUN=function(x) {which(!(is.na(x)) & (x != 0))})
          if (all(is.na(varsel))) {
            varsel.list[[k]] <- list(NA, NA)
            varsign.list[[k]] <- list(NA, NA)
          } else {
            varsel.list[[k]] <- varsel
            varsign.list[[k]] <- lapply(cv.coef, function(x) sign(x))
          }
        }
      }
    }
    k <- k + 1
  }
  
  # Stopping the cluster and cleaning all MPI states
  if (parallel) {
    parallel::stopCluster(cl=cluster)
  }
  
  # Get the fitted model
  if (all(is.na(varsel.list))) {
    # Failed to select any covariates
    varsel <- NA
    varsign <- NA
    enalpha.min <- NA
    enalpha.1se <- NA
    enlambda.min <- NA
    enlambda.1se <- NA
    msize <- rep(NA, M)
    success <- FALSE
  } else {
    # Get the selected covariates with their signs from all folds for each optimal lambda value
    varsel <- vector(mode="list", length=M)
    varsign <- vector(mode="list", length=M)
    for (l in 1:M) {
      sel.l <- numeric(0)
      sign.l <- numeric(0)
      for (k in 1:K) {
        sel.l <- c(sel.l, varsel.list[[k]][[l]])
        sign.l <- c(sign.l, varsign.list[[k]][[l]])
      }
      na <- (is.na(sel.l))
      nna <- length(which(na))
      sel.l <- sel.l[!na]
      sign.l <- sign.l[!na]
      if ((is.empty(sel.l)) || (is.empty(sign.l)) || (all(sign.l == 0)))  {
        varsel[[l]] <- NA
        varsign[[l]] <- NA
      } else {
        vartab <- table(names(sel.l), useNA="no")
        nvotes <- K - nna
        w <- names(which(vartab >= ceiling(nvotes*vscons)))
        if ((all(is.na(w))) || (is.empty(w))) {
          varsel[[l]] <- NA
          varsign[[l]] <- NA
        } else {
          varsel[[l]] <- pmatch(x=w, table=colnames(X), nomatch = NA_integer_, duplicates.ok = FALSE)
          names(varsel[[l]]) <- w
          signtab <- table(names(sign.l), sign.l, useNA="no")
          if (length(unique(sign.l)) == 1) {
            signfreq <- rep(x=unique(sign.l), times=nrow(signtab))
            names(signfreq) <- rownames(signtab)
            varsign[[l]] <- signfreq[w]
          } else if (all(sign.l != 1)) {
            signfreq <- apply(signtab[,c("-1","0"),drop=FALSE], 1, which.max)
            signfreq[signfreq==1] <- -1
            signfreq[signfreq==2] <- 0
            varsign[[l]] <- signfreq[w]
          } else if (all(sign.l != -1)) {
            signfreq <- apply(signtab[,c("0","1"),drop=FALSE], 1, which.max)
            signfreq[signfreq==1] <- 0
            signfreq[signfreq==2] <- 1
            varsign[[l]] <- signfreq[w]
          } else {
            signfreq <- apply(signtab[,c("-1","1"),drop=FALSE], 1, which.max)
            signfreq[signfreq==1] <- -1
            signfreq[signfreq==2] <- 1
            varsign[[l]] <- signfreq[w]
          }
          ord <- order(as.numeric(varsel[[l]]))
          varsel[[l]] <- varsel[[l]][ord]
          varsign[[l]] <- varsign[[l]][ord]
        }
      }
    }
    if (all(is.na(varsel))) {
      varsel <- NA
      varsign <- NA
      enalpha.min <- NA
      enalpha.1se <- NA
      enlambda.min <- NA
      enlambda.1se <- NA
      msize <- rep(NA, M)
      success <- FALSE
    } else {
      if (onese) {
        varsel <- varsel[[1]]
        varsign <- varsign[[1]]
        msize[1] <- length(varsel)
      } else {
        varsel <- varsel[[2]]
        varsign <- varsign[[2]]
        msize[2] <- length(varsel)
      }
      success <- TRUE
    }
  }
  
  # Returning updated argument `vsarg` of variable selection parameters
  vsarg <- list("tau"=tau,
                "alpha"=alpha,
                "nalpha"=nalpha,
                "nlambda"=nlambda,
                "msize"=msize)
  
  return(list("vsarg"=vsarg,
              "vscons"=vscons,
              "varsel"=varsel,
              "varsign"=varsign,
              "enalpha.min"=enalpha.min,
              "enalpha.1se"=enalpha.1se,
              "enlambda.min"=enlambda.min,
              "enlambda.1se"=enlambda.1se,
              "success"=success,
              "seed"=seed))
}
#===============================================================================================================================#



#===============================================================================================================================#
#
#===============================================================================================================================#

cv.pcqreg <- function(X,
                      y,
                      delta,
                      method,
                      gamma,
                      tau,
                      nfolds,
                      type.loss,
                      alpha,
                      nlambda,
                      lambda.min,
                      preprocess,
                      screen,
                      max.iter,
                      eps,
                      dfmax,
                      penalty.factor,
                      parallel,
                      cluster,
                      verbose,
                      seed) {
  
  cvf <- function(i, XX, y, delta, folds, cv.args, loss.args) {
    
    cv.args$X <- XX[folds$perm[(folds$which != i)], , drop=FALSE]
    cv.args$y <- y[folds$perm[(folds$which != i)]]
    cv.args$delta <- delta[folds$perm[(folds$which != i)]]
    X2 <- XX[folds$perm[(folds$which == i)], , drop=FALSE]
    y2 <- y[folds$perm[(folds$which == i)]]
    fit.i <- do.call("pcqreg", cv.args) # ensures the cross validation uses the same gamma for huber loss
    yhat <- matrix(data=predict.pcqreg(object=fit.i, X=X2, type="response", exact=FALSE),
                   nrow=length(y2))
    
    return(list(pe=loss.pcqreg(y=y2, yhat=yhat, args=loss.args),
                nl=length(fit.i$lambda)))
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n <- length(y)
  
  fit <- pcqreg(X=X,
                y=y,
                delta=delta,
                method=method,
                gamma=gamma,
                tau=tau,
                alpha=alpha,
                nlambda=nlambda,
                lambda.min=lambda.min,
                preprocess=preprocess,
                screen=screen,
                max.iter=max.iter,
                eps=eps,
                dfmax=dfmax,
                penalty.factor=penalty.factor,
                verbose=verbose)
  
  cv.args <- list("method"=method,
                  "gamma"=gamma,
                  "tau"=tau,
                  "alpha"=alpha,
                  "lambda"=fit$lambda,
                  "nlambda"=nlambda,
                  "lambda.min"=lambda.min,
                  "preprocess"=preprocess,
                  "screen"=screen,
                  "max.iter"=max.iter,
                  "eps"=eps,
                  "dfmax"=dfmax,
                  "penalty.factor"=penalty.factor,
                  "verbose"=verbose)
  
  loss.args <- list("method"=method,
                    "gamma"=gamma,
                    "tau"=tau,
                    "type.loss"=type.loss)
  
  if (is.null(delta)) {
    folds <- cv.folds(y=y, K=nfolds, seed=seed)
  } else {
    folds <- cv.folds(y=delta, K=nfolds, seed=seed)
  }
  
  E <- matrix(NA, nrow=n, ncol=length(cv.args$lambda))
  if (parallel) {
    parallel::clusterSetRNGStream(cl=cluster, iseed=seed)
    fold.results <- parallel::parLapply(cl=cluster,
                                        X=1:nfolds,
                                        fun=cvf,
                                        XX=X,
                                        y=y,
                                        delta=delta,
                                        folds=folds,
                                        cv.args=cv.args,
                                        loss.args=loss.args)
    for (i in 1:nfolds) {
      fit.i <- fold.results[[i]]
      E[folds$perm[(folds$which == i)], 1:fit.i$nl] <- fit.i$pe
    }
  } else {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    for (i in 1:nfolds) {
      fit.i <- cvf(i=i,
                   XX=X,
                   y=y,
                   delta=delta,
                   folds=folds,
                   cv.args=cv.args,
                   loss.args=loss.args)
      E[folds$perm[(folds$which == i)], 1:fit.i$nl] <- fit.i$pe
    }
  }
  
  ## Eliminate saturated lambda values
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[,ind]
  lambda <- cv.args$lambda[ind]
  
  ## Results
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, sd) / sqrt(n/nfolds)
  index.min <- which.min(cve)
  
  # adjust the selection using 1-SD method
  index.1se <- min(which(cve < cve[index.min]+cvse[index.min]))
  
  # Output
  return(list(cve=cve,
              cvse=cvse,
              type.loss=type.loss,
              lambda=lambda,
              fit=fit,
              lambda.1se=lambda[index.1se],
              lambda.min=lambda[index.min]))
}
#===============================================================================================================================#




#===============================================================================================================================#
#
#===============================================================================================================================#

loss.pcqreg <- function(y,
                        yhat,
                        args) {
  
  hloss <- function(r, gamma) {
    rr <- abs(r)
    ifelse(rr <= gamma, rr^2/(2*gamma), rr-gamma/2)
  }
  
  qloss <- function(r, tau) {
    ifelse(r <= 0, (tau-1)*r, tau*r)
  }
  
  r <- y - yhat
  type.loss <- args$type.loss
  if (type.loss == "deviance") {
    method <- args$method
    if (method == "huber") {
      gamma <- args$gamma
      val <- hloss(r, gamma)
    } else if ((method == "quantile") || (method == "cquantile")) {
      tau <- args$tau
      val <- qloss(r, tau)
    } else {
      val <- r^2
    }
  } else if (type.loss == "mse") {
    val <- r^2
  } else {
    val <- abs(r)
  }
  return(val)
}
#===============================================================================================================================#




#===============================================================================================================================#
#
#===============================================================================================================================#

pcqreg <- function (X,
                    y,
                    delta,
                    method,
                    gamma,
                    tau,
                    alpha,
                    nlambda,
                    lambda.min,
                    lambda,
                    preprocess,
                    screen,
                    max.iter,
                    eps,
                    dfmax,
                    penalty.factor,
                    verbose) {
  
  # Error checking
  if (missing(lambda) && nlambda < 2)
    stop("nlambda should be at least 2. Exiting ... \n\n")
  if (alpha < 0 || alpha > 1)
    stop("alpha should be between 0 and 1. Exiting ... \n\n")
  if (method == "huber" && !missing(gamma) && gamma <= 0)
    stop("gamma should be positive for Huber loss. Exiting ... \n\n")
  if ((method == "quantile" || method == "cquantile")  && (tau < 0 || tau > 1))
    stop("tau should be between 0 and 1 for quantile loss. Exiting ... \n\n")
  if (length(penalty.factor) != ncol(X))
    stop("the length of penalty.factor should equal the number of columns of X. Exiting ... \n\n")
  
  # Flag for user-supplied lambda
  if (missing(lambda)) {
    lambda <- double(nlambda)
    user <- 0
  } else {
    nlambda <- length(lambda)
    user <- 1
  }
  
  # Saving function call for output
  call <- match.call()
  
  # Include a column for intercept
  n <- nrow(X)
  XX <- cbind(rep(1,n), X)
  penalty.factor <- c(0, penalty.factor) # no penalty for intercept term
  p <- ncol(XX)
  
  shift <- 0
  if (method == "huber") {
    shift <- if(gamma > sd(y))
      mean(y)
    else
      median(y)
  } else if (method == "ls") {
    shift <- mean(y)
  } else if (method == "quantile") {
    shift <- quantile(y, tau)
  } else if (method == "cquantile") {
    obj <- tryCatch({coef(object=quantreg::crq(formula=survival::Surv(time=y, event=delta, type="right") ~ 1, taus=tau, method="Portnoy"),taus=tau)},
                    error=function(x){NULL})
    if ((is.na(obj)) || (is.empty(obj))) {
      shift <- quantile(y, tau)
    } else {
      shift <- coef(object=quantreg::crq(formula=survival::Surv(time=y, event=delta, type="right") ~ 1, taus=tau, method="Portnoy"),taus=tau)
    }
  }
  yy <- y - shift
  
  # Flags for preprocessing and screening
  ppflag = switch(preprocess, standardize = 1L, rescale = 2L)
  scrflag = switch(screen, ASR = 1L, SR = 2L, none = 0L)
  
  # Fitting
  if (alpha > 0) {
    if (method == "huber") {
      fit <- .C("C_huber", double(p*nlambda), integer(nlambda), as.double(lambda), integer(1), integer(1), as.double(XX), as.double(yy),
                as.double(penalty.factor), as.double(gamma), as.double(alpha), as.double(eps), as.double(lambda.min), as.integer(nlambda),
                as.integer(n), as.integer(p), as.integer(ppflag), as.integer(scrflag), 1L, as.integer(dfmax), as.integer(max.iter), as.integer(user), as.integer(verbose))
    } else if ((method == "quantile") || (method == "cquantile")) {
      fit <- .C("C_quantile", double(p*nlambda), integer(nlambda), as.double(lambda), integer(1), integer(1), as.double(XX), as.double(yy),
                as.double(penalty.factor), as.double(tau), as.double(alpha), as.double(eps), as.double(lambda.min), as.integer(nlambda),
                as.integer(n), as.integer(p), as.integer(ppflag), as.integer(scrflag), 1L, as.integer(dfmax), as.integer(max.iter), as.integer(user), as.integer(verbose))
    } else {
      fit <- .C("C_squared", double(p*nlambda), integer(nlambda), as.double(lambda), integer(1), integer(1), as.double(XX), as.double(yy),
                as.double(penalty.factor), as.double(alpha), as.double(eps), as.double(lambda.min), as.integer(nlambda),
                as.integer(n), as.integer(p), as.integer(ppflag), as.integer(scrflag), 1L, as.integer(dfmax), as.integer(max.iter), as.integer(user), as.integer(verbose))
    }
    beta <- matrix(fit[[1]],nrow = p)
    iter <- fit[[2]]
    lambda <- fit[[3]]
    saturated <- fit[[4]]
    nv <- fit[[5]]
    # Eliminate saturated lambda values
    ind <- !is.na(iter)
    beta <- beta[, ind]
    iter <- iter[ind]
    lambda <- lambda[ind]
  } else {
    if (method == "huber") {
      fit <- .C("C_huber_l2", double(p*nlambda), integer(nlambda), as.double(lambda), as.double(XX), as.double(yy),
                as.double(penalty.factor), as.double(gamma), as.double(eps), as.double(lambda.min), as.integer(nlambda),
                as.integer(n), as.integer(p), as.integer(ppflag), 1L, as.integer(max.iter), as.integer(user), as.integer(verbose))
    } else if ((method == "quantile") || (method == "cquantile")) {
      fit <- .C("C_quantile_l2", double(p*nlambda), integer(nlambda), as.double(lambda), as.double(XX), as.double(yy),
                as.double(penalty.factor), as.double(tau), as.double(eps), as.double(lambda.min), as.integer(nlambda),
                as.integer(n), as.integer(p), as.integer(ppflag), 1L, as.integer(max.iter), as.integer(user), as.integer(verbose))
    } else {
      fit <- .C("C_squared_l2", double(p*nlambda), integer(nlambda), as.double(lambda), as.double(XX), as.double(yy),
                as.double(penalty.factor), as.double(eps), as.double(lambda.min), as.integer(nlambda), as.integer(n),
                as.integer(p), as.integer(ppflag), 1L, as.integer(max.iter), as.integer(user), as.integer(verbose))
    }
    beta <- matrix(fit[[1]],nrow = p)
    iter <- fit[[2]]
    lambda <- fit[[3]]
    saturated <- 0
    nv <- 0
  }
  
  # Intercept
  beta[1,] <- beta[1,] + shift
  
  # Names
  vnames <- colnames(X)
  if (is.null(vnames)) vnames=paste0("V", seq(p-1))
  vnames <- c("(Intercept)", vnames)
  dimnames(beta) <- list(vnames, paste0("L", 1:length(lambda)))
  
  # Output
  return(list(call=call,
              beta=beta,
              iter=iter,
              saturated=saturated,
              lambda=lambda,
              alpha=alpha,
              gamma=if (method == "huber") gamma else NULL,
              tau=if ((method == "quantile") || (method == "cquantile")) tau else NULL,
              penalty.factor=penalty.factor[-1],
              method=method,
              nv=nv))
}
#===============================================================================================================================#




#===============================================================================================================================#
#
#===============================================================================================================================#

predict.pcqreg <- function(object,
                           X,
                           lambda,
                           type=c("response","coefficients","nvars"),
                           exact) {
  
  type <- match.arg(type)
  
  if ((missing(X)) && (type == "response"))
    stop("Need to supply 'X'")
  
  beta <- coef.pcqreg(object=object, lambda=lambda, exact=exact)
  if (type == "coefficients")
    return(beta)
  
  if (is.matrix(beta)) {
    b0 <- beta[1,]
    b <- beta[-1,]
  } else {
    b0 <- beta[1]
    b <- beta[-1]
  }
  if (type == "nvars") {
    if (is.matrix(b))
      return(apply(b!=0, 2, sum))
    else
      return(sum(b!=0))
  }
  if (type == "response")
    return(sweep(X %*% b, 2, b0, "+"))
  
  return(NULL)
}
#===============================================================================================================================#




#===============================================================================================================================#
#
#===============================================================================================================================#

coef.pcqreg <- function(object,
                        lambda,
                        exact) {
  
  if (missing(lambda)) {
    beta <- object$beta
  } else if (exact) {
    # augment the lambda sequence with the new values, and refit the model
    ls <- object$lambda
    ind <- match(lambda, ls, 0)
    if (any(ind == 0)) {
      ls <- unique(rev(sort(c(lambda,ls))))
      object <- update(object, lambda=ls)
      ind <- match(lambda, ls)
    }
    beta <- object$beta[, ind]
  } else {
    # use linear interpolation to estimate coefficients for supplied lambda
    ls <- object$lambda
    lambda[lambda>max(ls)] <- max(ls)
    lambda[lambda<min(ls)] <- min(ls)
    ind <- approx(x=ls, y=seq(ls), xout=lambda)$y
    left <- floor(ind)
    right <- ceiling(ind)
    weight <- ind %% 1
    beta <- (1-weight)*object$beta[,left] + weight*object$beta[,right]
    if (length(lambda) > 1) colnames(beta) <- round(lambda, 4)
  }
  return(beta)
}
#===============================================================================================================================#




#===============================================================================================================================#
# Variables screening by Penalized Partial-Likelihood (PPL)
# CV of variable selection using full cross-validation of both parameters alpha and lambda
#===============================================================================================================================#

cv.ppl <- function(X,
                   y,
                   delta,
                   K,
                   vsarg,
                   vscons,
                   cv,
                   onese,
                   parallel,
                   conf,
                   verbose,
                   seed) {
  
  # Parsing and evaluating 'vsarg' argument to evaluate all parameters
  alpha <- NULL
  nalpha <- NULL
  nlambda <- NULL
  eval(parse( text=unlist(strsplit(x=vsarg, split=",")) ))
  
  M <- 2
  msize <- numeric(M)
  
  if (is.null(alpha)) {
    enalpha <- seq(from=0, to=1, length.out=nalpha)
  } else {
    nalpha <- 1
    enalpha <- alpha
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Parallel back-end registration
  if (parallel) {
    if (conf$type == "SOCKET") {
      cluster <- parallel::makeCluster(spec=conf$spec,
                                       type="PSOCK",
                                       homogeneous=conf$homo,
                                       outfile=file.path(getwd(), "output_SOCK_inner.txt", fsep=.Platform$file.sep),
                                       verbose=conf$verbose)
    } else if (conf$type == "MPI") {
      cluster <- parallel::makeCluster(spec=conf$spec,
                                       type="MPI",
                                       homogeneous=conf$homo,
                                       outfile=file.path(getwd(), "output_MPI_inner.txt", fsep=.Platform$file.sep),
                                       verbose=conf$verbose)
    }
  }
  
  folds <- cv.folds(y=delta, K=K, seed=seed)
  varsel.list <- vector(mode="list", length=K)
  varsign.list <- vector(mode="list", length=K)
  k <- 1
  while (k <= K) {
    if ((verbose) && (cv)) cat("Fold : ", k, "\n", sep="")
    # Initialize training and test data
    if (K == 1) {
      traindata <- testdata <- X[folds$perm[(folds$which == k)], , drop=FALSE]
      traintime <- testtime <- y[folds$perm[(folds$which == k)]]
      trainstatus <- teststatus <- delta[folds$perm[(folds$which == k)]]
    } else {
      traindata <- X[folds$perm[(folds$which != k)], , drop=FALSE]
      traintime <- y[folds$perm[(folds$which != k)]]
      trainstatus <- delta[folds$perm[(folds$which != k)]]
      testdata <- X[folds$perm[(folds$which == k)], , drop=FALSE]
      testtime <- y[folds$perm[(folds$which == k)]]
      teststatus <- delta[folds$perm[(folds$which == k)]]
    }
    enlambda <- vector(mode="list", length=nalpha)
    cv.errmu <- vector(mode="list", length=nalpha)
    cv.errse <- vector(mode="list", length=nalpha)
    
    for (i in 1:nalpha) {
      cv.fit <- glmnet::cv.glmnet(x=traindata,
                                  y=survival::Surv(time=traintime, event=trainstatus),
                                  alpha=enalpha[i],
                                  nlambda=nlambda,
                                  nfolds=pmax(3,K),
                                  family="cox",
                                  parallel=parallel,
                                  maxit=1e5)
      cv.errmu[[i]] <- cv.fit$cvm
      cv.errse[[i]] <- cv.fit$cvsd
      enlambda[[i]] <- cv.fit$lambda
    }
    cv.errmu <- list2mat(list=cv.errmu, coltrunc="max", fill=NA)
    cv.errse <- list2mat(list=cv.errse, coltrunc="max", fill=NA)
    cv.errmu.index <- as.matrix(which(x=(cv.errmu == as.numeric(cv.errmu)[which.min(cv.errmu)]), arr.ind=TRUE, useNames=FALSE))
    w <- which.min(cv.errmu.index[, 1, drop=TRUE])
    index.min <- as.numeric(cv.errmu.index[w,,drop=TRUE])
    ww <- as.matrix(which(x=cv.errmu <= cv.errmu[index.min[1],index.min[2]] + cv.errse[index.min[1],index.min[2]], arr.ind=TRUE, useNames=TRUE))
    index.1se <- c(min(ww[,1]), min(ww[ww[,1]==min(ww[,1]),2]))
    if (is.empty(index.min) || is.empty(index.1se)) {
      varsel.list[[k]] <- list(NA, NA)
      varsign.list[[k]] <- list(NA, NA)
    } else {
      enalpha.min <- enalpha[index.min[1]]
      enalpha.1se <- enalpha[index.1se[1]]
      enlambda.min <- enlambda[[index.min[1]]][index.min[2]]
      enlambda.1se <- enlambda[[index.1se[1]]][index.1se[2]]
      fit <- glmnet::glmnet(x=traindata,
                            y=survival::Surv(time=traintime, event=trainstatus),
                            alpha=ifelse(test=onese, yes=enalpha.1se, no=enalpha.min),
                            family="cox",
                            maxit=1e5)
      w <- apply(X=fit$beta, MARGIN=2, FUN=function(x) {sum(!(is.na(x)) & (x != 0))})
      if (all(w == 0)) {
        varsel.list[[k]] <- list(NA, NA)
        varsign.list[[k]] <- list(NA, NA)
      } else {
        cv.coef <- coef(object=fit, s=c(enlambda.1se, enlambda.min))
        cv.coef <- list(cv.coef[,1], cv.coef[,2])
        varsel <- lapply(cv.coef, FUN=function(x) {which(!(is.na(x)) & (x != 0))})
        if (all(is.na(varsel))) {
          varsel.list[[k]] <- list(NA, NA)
          varsign.list[[k]] <- list(NA, NA)
        } else {
          varsel.list[[k]] <- varsel
          varsign.list[[k]] <- lapply(cv.coef, function(x) sign(x))
        }
      }
    }
    k <- k + 1
  }
  
  # Stopping the cluster and cleaning all MPI states
  if (parallel) {
    parallel::stopCluster(cl=cluster)
  }
  
  # Get the fitted model
  if (all(is.na(varsel.list))) {
    # Failed to select any covariates
    varsel <- NA
    varsign <- NA
    enalpha.min <- NA
    enalpha.1se <- NA
    enlambda.min <- NA
    enlambda.1se <- NA
    msize <- rep(NA, M)
    success <- FALSE
  } else {
    # Get the selected covariates with their signs from all folds for each optimal lambda value
    varsel <- vector(mode="list", length=M)
    varsign <- vector(mode="list", length=M)
    for (l in 1:M) {
      sel.l <- numeric(0)
      sign.l <- numeric(0)
      for (k in 1:K) {
        sel.l <- c(sel.l, varsel.list[[k]][[l]])
        sign.l <- c(sign.l, varsign.list[[k]][[l]])
      }
      na <- (is.na(sel.l))
      nna <- length(which(na))
      sel.l <- sel.l[!na]
      sign.l <- sign.l[!na]
      if ((is.empty(sel.l)) || (is.empty(sign.l)) || (all(sign.l == 0)))  {
        varsel[[l]] <- NA
        varsign[[l]] <- NA
      } else {
        vartab <- table(names(sel.l), useNA="no")
        nvotes <- K - nna
        w <- names(which(vartab >= ceiling(nvotes*vscons)))
        if ((all(is.na(w))) || (is.empty(w))) {
          varsel[[l]] <- NA
          varsign[[l]] <- NA
        } else {
          varsel[[l]] <- pmatch(x=w, table=colnames(X), nomatch = NA_integer_, duplicates.ok = FALSE)
          names(varsel[[l]]) <- w
          signtab <- table(names(sign.l), sign.l, useNA="no")
          if (length(unique(sign.l)) == 1) {
            signfreq <- rep(x=unique(sign.l), times=nrow(signtab))
            names(signfreq) <- rownames(signtab)
            varsign[[l]] <- signfreq[w]
          } else if (all(sign.l != 1)) {
            signfreq <- apply(signtab[,c("-1","0"),drop=FALSE], 1, which.max)
            signfreq[signfreq==1] <- -1
            signfreq[signfreq==2] <- 0
            varsign[[l]] <- signfreq[w]
          } else if (all(sign.l != -1)) {
            signfreq <- apply(signtab[,c("0","1"),drop=FALSE], 1, which.max)
            signfreq[signfreq==1] <- 0
            signfreq[signfreq==2] <- 1
            varsign[[l]] <- signfreq[w]
          } else {
            signfreq <- apply(signtab[,c("-1","1"),drop=FALSE], 1, which.max)
            signfreq[signfreq==1] <- -1
            signfreq[signfreq==2] <- 1
            varsign[[l]] <- signfreq[w]
          }
          ord <- order(as.numeric(varsel[[l]]))
          varsel[[l]] <- varsel[[l]][ord]
          varsign[[l]] <- varsign[[l]][ord]
        }
      }
    }
    if (all(is.na(varsel))) {
      varsel <- NA
      varsign <- NA
      enalpha.min <- NA
      enalpha.1se <- NA
      enlambda.min <- NA
      enlambda.1se <- NA
      msize <- rep(NA, M)
      success <- FALSE
    } else {
      if (onese) {
        varsel <- varsel[[1]]
        varsign <- varsign[[1]]
        msize[1] <- length(varsel)
      } else {
        varsel <- varsel[[2]]
        varsign <- varsign[[2]]
        msize[2] <- length(varsel)
      }
      success <- TRUE
    }
  }
  
  # Returning updated argument `vsarg` of variable selection parameters
  vsarg <- list("alpha"=alpha,
                "nalpha"=nalpha,
                "nlambda"=nlambda,
                "msize"=msize)
  
  return(list("vsarg"=vsarg,
              "vscons"=vscons,
              "varsel"=varsel,
              "varsign"=varsign,
              "enalpha.min"=enalpha.min,
              "enalpha.1se"=enalpha.1se,
              "enlambda.min"=enlambda.min,
              "enlambda.1se"=enlambda.1se,
              "success"=success,
              "seed"=seed))
}
#===============================================================================================================================#




#===============================================================================================================================#
# Variables screening by SPCA (SPCA)
# CV over number of SPCs and ordered CPH coefficients
#===============================================================================================================================#

cv.spca <- function(X,
                    y,
                    delta,
                    K,
                    vsarg,
                    vscons,
                    cv,
                    parallel,
                    conf,
                    verbose,
                    seed) {
  
  # Parsing and evaluating 'vsarg' argument to evaluate all parameters
  n.thres <- NULL
  n.pcs <- NULL
  n.var <- NULL
  eval(parse( text=unlist(strsplit(x=vsarg, split=",")) ))
  
  M <- 1
  msize <- numeric(M)
  
  if (ncol(X) <= n.var) {
    n.var <- pmin(ncol(X), n.var) - 1
    cat("Warning: Parameter `n.var` was greater or equal than the dimensionality of covariates and was reset to ", ncol(X) - 1, ".\n", sep="")
  }
  
  if (ncol(X) <= n.pcs) {
    n.pcs <- pmin(ncol(X), n.pcs) - 1
    cat("Warning: Parameter `n.pcs` was greater or equal than the dimensionality of covariates and was reset to ", ncol(X) - 1, ".\n", sep="")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Parallel back-end registration
  if (parallel) {
    if (conf$type == "SOCKET") {
      cluster <- parallel::makeCluster(spec=conf$spec,
                                       type="PSOCK",
                                       homogeneous=conf$homo,
                                       outfile=file.path(getwd(), "output_SOCK_inner.txt", fsep=.Platform$file.sep),
                                       verbose=conf$verbose)
    } else if (conf$type == "MPI") {
      cluster <- parallel::makeCluster(spec=conf$spec,
                                       type="MPI",
                                       homogeneous=conf$homo,
                                       outfile=file.path(getwd(), "output_MPI_inner.txt", fsep=.Platform$file.sep),
                                       verbose=conf$verbose)
    }
  }
  
  folds <- cv.folds(y=delta, K=K, seed=seed)
  varsel.list <- vector(mode="list", length=K)
  varsign.list <- vector(mode="list", length=K)
  k <- 1
  while (k <= K) {
    if ((verbose) && (cv)) cat("Fold : ", k, "\n", sep="")
    # Initialize training and test data
    if (K == 1) {
      traindata <- testdata <- X[folds$perm[(folds$which == k)], , drop=FALSE]
      traintime <- testtime <- y[folds$perm[(folds$which == k)]]
      trainstatus <- teststatus <- delta[folds$perm[(folds$which == k)]]
    } else {
      traindata <- X[folds$perm[(folds$which != k)], , drop=FALSE]
      traintime <- y[folds$perm[(folds$which != k)]]
      trainstatus <- delta[folds$perm[(folds$which != k)]]
      testdata <- X[folds$perm[(folds$which == k)], , drop=FALSE]
      testtime <- y[folds$perm[(folds$which == k)]]
      teststatus <- delta[folds$perm[(folds$which == k)]]
    }
    
    # Structure the data
    pca.train.data <- list(x=t(traindata), y=traintime, censoring.status=trainstatus)
    pca.test.data  <- list(x=t(testdata),  y=testtime,  censoring.status=teststatus)
    
    # Compute trained Wald scores for each variable
    ws.train <- superpc::superpc.train(data=pca.train.data, type="survival", s0.perc=NULL)
    lower <- quantile(x=abs(ws.train$feature.scores), probs = 0)
    upper <- quantile(x=abs(ws.train$feature.scores), probs = 1 - (n.var/ncol(traindata)))
    var.thres <- seq(from=lower, to=upper, length=n.thres)
    
    if (cv) {
      # Finding the optimal trained PCR model by internal CV of the Likelihood Ratio Statistic (LRS)
      # Return LRS from full CV w.r.t. threshold values and #PCs
      superpccv <- superpc::superpc.cv(fit=ws.train,
                                       data=pca.train.data,
                                       n.threshold=n.thres,
                                       n.fold=K,
                                       n.components=n.pcs,
                                       min.features=n.var,
                                       max.features=ncol(traindata))
      lrs <- superpccv$scor
      max.mat <- which(lrs == max(lrs, na.rm=TRUE), arr.ind=TRUE)
      max.mat <- max.mat[1,,drop=FALSE]
      max.npcs <- max.mat[1]
      max.thr <- max.mat[2]
    } else {
      max.npcs <- n.pcs[length(n.pcs)]
      max.thr <- n.thres[length(n.thres)]
    }
    
    # Trained model used to predict the survival times on the test set
    superpcfit <- superpc::superpc.predict(object=ws.train,
                                           data=pca.train.data,
                                           newdata=pca.test.data,
                                           prediction.type="continuous",
                                           threshold=var.thres[max.thr],
                                           n.components=max.npcs)
    
    # List of selected variables (exceeding threshold) used in each fold
    varsel <- which(superpcfit$which.features)
    if (is.empty(varsel)) {
      varsel.list[[k]] <- NA
      varsign.list[[k]] <- NA
    } else {
      varsel.list[[k]] <- varsel
      # Feature selection for supervised principal components
      # Forms reduced model to approximate the supervised principal component predictor.
      ft <- superpc::superpc.predict.red(fit=ws.train,
                                         data=pca.train.data,
                                         data.test=pca.test.data,
                                         threshold=var.thres[max.thr],
                                         n.components=max.npcs,
                                         n.shrinkage=n.thres,
                                         prediction.type="continuous")
      # Importance scores of selected features for PC#1
      scores <- ft$import[varsel,1]
      # Cox scores for selected features in order of decreasing absolute value of importance score for PC#1
      lt <- superpc::superpc.listfeatures(data=pca.train.data,
                                          train.obj=ws.train,
                                          fit.red=ft,
                                          component.number=1)
      ordscor <- order(abs(scores), decreasing = TRUE)
      rawscor <- as.numeric(lt[,"Raw-score"])
      names(rawscor) <- names(scores[ordscor])
      varsign.list[[k]] <- sign(rawscor)
      # Determining the separating threshold of median survival time and formation of the two prognostic groups
      pred.pcr <- superpcfit$v.pred.1df
      pred.ind <- (pred.pcr <= median(pred.pcr))
    }
    k <- k + 1
  }
  
  # Stopping the cluster and cleaning all MPI states
  if (parallel) {
    parallel::stopCluster(cl=cluster)
  }
  
  # Get the fitted model
  if (all(is.na(varsel.list))) {
    # Failed to select any covariates
    varsel <- NA
    varsign <- NA
    msize <- NA
    success <- FALSE
  } else {
    # Get the selected covariates with their signs from all folds
    sel.l <- unlist(varsel.list)
    sign.l <- unlist(varsign.list)
    na <- (is.na(sel.l))
    nna <- length(which(na))
    sel.l <- sel.l[!na]
    sign.l <- sign.l[!na]
    if ((is.empty(sel.l)) || (is.empty(sign.l)) || (all(sign.l == 0)))  {
      varsel <- NA
      varsign <- NA
    } else {
      vartab <- table(names(sel.l), useNA="no")
      nvotes <- K - nna
      w <- names(which(vartab >= ceiling(nvotes*vscons)))
      if ((all(is.na(w))) || (is.empty(w))) {
        varsel <- NA
        varsign <- NA
      } else{
        varsel <- pmatch(x=w, table=colnames(X), nomatch = NA_integer_, duplicates.ok = FALSE)
        names(varsel) <- w
        signtab <- table(names(sign.l), sign.l, useNA="no")
        if (length(unique(sign.l)) == 1) {
          signfreq <- rep(x=unique(sign.l), times=nrow(signtab))
          names(signfreq) <- rownames(signtab)
          varsign <- signfreq[w]
        } else if (all(sign.l != 1)) {
          signfreq <- apply(signtab[,c("-1","0"),drop=FALSE], 1, which.max)
          signfreq[signfreq==1] <- -1
          signfreq[signfreq==2] <- 0
          varsign <- signfreq[w]
        } else if (all(sign.l != -1)) {
          signfreq <- apply(signtab[,c("0","1"),drop=FALSE], 1, which.max)
          signfreq[signfreq==1] <- 0
          signfreq[signfreq==2] <- 1
          varsign <- signfreq[w]
        } else {
          signfreq <- apply(signtab[,c("-1","1"),drop=FALSE], 1, which.max)
          signfreq[signfreq==1] <- -1
          signfreq[signfreq==2] <- 1
          varsign <- signfreq[w]
        }
        ord <- order(as.numeric(varsel))
        varsel <- varsel[ord]
        varsign <- varsign[ord]
      }
    }
    if (all(is.na(varsel))) {
      varsel <- NA
      varsign <- NA
      msize <- NA
      success <- FALSE
    } else {
      msize <- length(varsel)
      success <- TRUE
    }
  }
  
  # Returning updated argument `vsarg` of variable selection parameters
  vsarg <- list("n.thres"=n.thres,
                "n.pcs"=n.pcs,
                "n.var"=n.var,
                "msize"=msize)
  
  return(list("vsarg"=vsarg,
              "vscons"=vscons,
              "varsel"=varsel,
              "varsign"=varsign,
              "success"=success,
              "seed"=seed))
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                   cv.box(X,
#                          y,
#                          delta,
#                          B,
#                          K,
#                          cv,
#                          cvtype,
#                          cvarg,
#                          varsign,
#                          initcutpts,
#                          groups,
#                          decimals,
#                          probval,
#                          timeval,
#                          parallel,
#                          conf,
#                          verbose,
#                          seed)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

cv.box <- function(X,
                   y,
                   delta,
                   B,
                   K,
                   cv,
                   cvtype,
                   cvarg,
                   varsign,
                   initcutpts,
                   groups,
                   decimals,
                   probval,
                   timeval,
                   parallel,
                   conf,
                   verbose,
                   seed) {
  
  seed <- seed[1]
  replic <- 1:B
  
  if (!parallel) {
    
    if (!is.null(seed)) {
      seed <- (0:(B-1)) + seed
    }
    
    CV.type.box.obj <- cv.type.box(X=X,
                                   y=y,
                                   delta=delta,
                                   replic=replic,
                                   K=K,
                                   cvarg=cvarg,
                                   cv=cv,
                                   cvtype=cvtype,
                                   varsign=varsign,
                                   initcutpts=initcutpts,
                                   groups=groups,
                                   decimals=decimals,
                                   probval=probval,
                                   timeval=timeval,
                                   parallel=parallel,
                                   verbose=verbose,
                                   seed=seed)
    
  } else {
    
    # Parallel back-end registration
    if (conf$type == "SOCKET") {
      cluster <- parallel::makeCluster(spec=conf$spec,
                                       type="PSOCK",
                                       homogeneous=conf$homo,
                                       outfile=conf$outfile,
                                       verbose=conf$verbose)
    } else if (conf$type == "MPI") {
      cluster <- parallel::makeCluster(spec=conf$spec,
                                       type="MPI",
                                       homogeneous=conf$homo,
                                       outfile=conf$outfile,
                                       verbose=conf$verbose)
    } else {
      stop("Unrecognized cluster type: you must specify a \"SOCKET\" or \"MPI\" cluster type. Exiting ... \n\n")
    }
    
    # Parallelization
    parallel::clusterSetRNGStream(cl=cluster, iseed=seed)
    obj.cl <- parallel::clusterApply(cl=cluster, x=replic, fun=cv.type.box,
                                     X=X,
                                     y=y,
                                     delta=delta,
                                     K=K,
                                     cvarg=cvarg,
                                     cv=cv,
                                     cvtype=cvtype,
                                     varsign=varsign,
                                     initcutpts=initcutpts,
                                     groups=groups,
                                     decimals=decimals,
                                     probval=probval,
                                     timeval=timeval,
                                     parallel=parallel,
                                     verbose=verbose,
                                     seed=NULL)
    
    # Stopping the cluster and cleaning all MPI states
    parallel::stopCluster(cl=cluster)
    
    # Collect the results for each step from all the replicates
    CV.type.box.obj <- list("cv.maxsteps"=numeric(0),
                            "cv.trace"=vector(mode="list", length=B),
                            "cv.boxind"=vector(mode="list", length=B),
                            "cv.boxcut"=vector(mode="list", length=B),
                            "cv.support"=vector(mode="list", length=B),
                            "cv.size"=vector(mode="list", length=B),
                            "cv.bmp.lhr"=vector(mode="list", length=B),
                            "cv.bmp.lrt"=vector(mode="list", length=B),
                            "cv.bmp.cer"=vector(mode="list", length=B),
                            "cv.grp.lhr"=vector(mode="list", length=B),
                            "cv.grp.lrt"=vector(mode="list", length=B),
                            "cv.grp.cer"=vector(mode="list", length=B),
                            "cv.time.bar"=vector(mode="list", length=B),
                            "cv.prob.bar"=vector(mode="list", length=B),
                            "cv.max.time.bar"=vector(mode="list", length=B),
                            "cv.min.prob.bar"=vector(mode="list", length=B),
                            "success"=logical(B))
    for (b in 1:B) {
      CV.type.box.obj$cv.maxsteps <- c(CV.type.box.obj$cv.maxsteps, obj.cl[[b]]$cv.maxsteps)
      CV.type.box.obj$cv.trace[[b]] <- obj.cl[[b]]$cv.trace
      CV.type.box.obj$cv.boxind[[b]] <- obj.cl[[b]]$cv.boxind
      CV.type.box.obj$cv.boxcut[[b]] <- obj.cl[[b]]$cv.boxcut
      CV.type.box.obj$cv.support[[b]] <- obj.cl[[b]]$cv.support
      CV.type.box.obj$cv.size[[b]] <- obj.cl[[b]]$cv.size
      CV.type.box.obj$cv.bmp.lhr[[b]] <- obj.cl[[b]]$cv.bmp.lhr
      CV.type.box.obj$cv.bmp.lrt[[b]] <- obj.cl[[b]]$cv.bmp.lrt
      CV.type.box.obj$cv.bmp.cer[[b]] <- obj.cl[[b]]$cv.bmp.cer
      CV.type.box.obj$cv.grp.lhr[[b]] <- obj.cl[[b]]$cv.grp.lhr
      CV.type.box.obj$cv.grp.lrt[[b]] <- obj.cl[[b]]$cv.grp.lrt
      CV.type.box.obj$cv.grp.cer[[b]] <- obj.cl[[b]]$cv.grp.cer
      CV.type.box.obj$cv.time.bar[[b]] <- obj.cl[[b]]$cv.time.bar
      CV.type.box.obj$cv.prob.bar[[b]] <- obj.cl[[b]]$cv.prob.bar
      CV.type.box.obj$cv.max.time.bar[[b]] <- obj.cl[[b]]$cv.max.time.bar
      CV.type.box.obj$cv.min.prob.bar[[b]] <- obj.cl[[b]]$cv.min.prob.bar
      CV.type.box.obj$success[b] <- obj.cl[[b]]$success
    }
    
  }
  
  return(list("cv.maxsteps"=CV.type.box.obj$cv.maxsteps,
              "cv.trace"=CV.type.box.obj$cv.trace,
              "cv.boxind"=CV.type.box.obj$cv.boxind,
              "cv.boxcut"=CV.type.box.obj$cv.boxcut,
              "cv.support"=CV.type.box.obj$cv.support,
              "cv.size"=CV.type.box.obj$cv.size,
              "cv.bmp.lhr"=CV.type.box.obj$cv.bmp.lhr,
              "cv.bmp.lrt"=CV.type.box.obj$cv.bmp.lrt,
              "cv.bmp.cer"=CV.type.box.obj$cv.bmp.cer,
              "cv.grp.lhr"=CV.type.box.obj$cv.grp.lhr,
              "cv.grp.lrt"=CV.type.box.obj$cv.grp.lrt,
              "cv.grp.cer"=CV.type.box.obj$cv.grp.cer,
              "cv.time.bar"=CV.type.box.obj$cv.time.bar,
              "cv.prob.bar"=CV.type.box.obj$cv.prob.bar,
              "cv.max.time.bar"=CV.type.box.obj$cv.max.time.bar,
              "cv.min.prob.bar"=CV.type.box.obj$cv.min.prob.bar,
              "success"=all(CV.type.box.obj$success),
              "seed"=seed))
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                   cv.type.box(X,
#                               y,
#                               delta,
#                               replic,
#                               K,
#                               cv,
#                               cvtype,
#                               cvarg,
#                               varsign,
#                               initcutpts,
#                               decimals,
#                               groups,
#                               probval,
#                               timeval,
#                               parallel,
#                               verbose,
#                               seed)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

cv.type.box <- function(X,
                        y,
                        delta,
                        replic,
                        K,
                        cv,
                        cvtype,
                        cvarg,
                        varsign,
                        initcutpts,
                        groups,
                        decimals,
                        probval,
                        timeval,
                        parallel,
                        verbose,
                        seed) {
  
  if (!parallel) {
    
    B <- length(replic)
    success <- logical(B)
    CV.maxsteps <- numeric(B)
    CV.boxind <- vector(mode="list", length=B)
    CV.boxcut <- vector(mode="list", length=B)
    CV.support <- vector(mode="list", length=B)
    CV.size<- vector(mode="list", length=B)
    CV.bmp.lhr <- vector(mode="list", length=B)
    CV.bmp.lrt <- vector(mode="list", length=B)
    CV.bmp.cer <- vector(mode="list", length=B)
    CV.grp.lhr <- vector(mode="list", length=B)
    CV.grp.lrt <- vector(mode="list", length=B)
    CV.grp.cer <- vector(mode="list", length=B)
    CV.trace <- vector(mode="list", length=B)
    CV.time.bar <- vector(mode="list", length=B)
    CV.prob.bar <- vector(mode="list", length=B)
    CV.max.time.bar <- vector(mode="list", length=B)
    CV.min.prob.bar <- vector(mode="list", length=B)
    
    b <- 1
    while (b <= B) {
      if (verbose) {
        cat("Replicate : ", b, "\n", sep="")
        cat("seed : ", seed[b], "\n", sep="")
      }
      if (cvtype == "averaged") {
        CVBOX <- cv.ave.box(X=X,
                            y=y,
                            delta=delta,
                            K=K,
                            cv=cv,
                            cvarg=cvarg,
                            groups=groups,
                            decimals=decimals,
                            probval=probval,
                            timeval=timeval,
                            varsign=varsign,
                            initcutpts=initcutpts,
                            verbose=verbose,
                            seed=seed[b])
      } else if (cvtype == "combined") {
        CVBOX <- cv.comb.box(X=X,
                             y=y,
                             delta=delta,
                             K=K,
                             cv=cv,
                             cvarg=cvarg,
                             groups=groups,
                             decimals=decimals,
                             probval=probval,
                             timeval=timeval,
                             varsign=varsign,
                             initcutpts=initcutpts,
                             verbose=verbose,
                             seed=seed[b])
      } else {
        stop("Invalid CV type option. Exiting ... \n\n")
      }
      success[b] <- CVBOX$success
      CV.maxsteps[b] <- CVBOX$cvfit$cv.maxsteps
      CV.trace[[b]] <- CVBOX$cvfit$cv.trace
      CV.boxind[[b]] <- CVBOX$cvfit$cv.boxind
      CV.boxcut[[b]] <- CVBOX$cvfit$cv.boxcut
      CV.support[[b]] <- CVBOX$cvfit$cv.stats$cv.support
      CV.size[[b]] <- CVBOX$cvfit$cv.stats$cv.size
      CV.bmp.lhr[[b]] <- CVBOX$cvfit$cv.stats$cv.bmp.lhr
      CV.bmp.lrt[[b]] <- CVBOX$cvfit$cv.stats$cv.bmp.lrt
      CV.bmp.cer[[b]] <- CVBOX$cvfit$cv.stats$cv.bmp.cer
      CV.grp.lhr[[b]] <- CVBOX$cvfit$cv.stats$cv.grp.lhr
      CV.grp.lrt[[b]] <- CVBOX$cvfit$cv.stats$cv.grp.lrt
      CV.grp.cer[[b]] <- CVBOX$cvfit$cv.stats$cv.grp.cer
      CV.time.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.time.bar
      CV.prob.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.prob.bar
      CV.max.time.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.max.time.bar
      CV.min.prob.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.min.prob.bar
      b <- b + 1
    }
    
  } else {
    
    if (cvtype == "averaged") {
      CVBOX <- cv.ave.box(X=X,
                          y=y,
                          delta=delta,
                          K=K,
                          cv=cv,
                          cvarg=cvarg,
                          groups=groups,
                          decimals=decimals,
                          probval=probval,
                          timeval=timeval,
                          varsign=varsign,
                          initcutpts=initcutpts,
                          verbose=verbose,
                          seed=seed)
    } else if (cvtype == "combined") {
      CVBOX <- cv.comb.box(X=X,
                           y=y,
                           delta=delta,
                           K=K,
                           cv=cv,
                           cvarg=cvarg,
                           groups=groups,
                           decimals=decimals,
                           probval=probval,
                           timeval=timeval,
                           varsign=varsign,
                           initcutpts=initcutpts,
                           verbose=verbose,
                           seed=seed)
    } else {
      stop("Invalid CV type option. Exiting ... \n\n")
    }
    success <- CVBOX$success
    CV.maxsteps <- CVBOX$cvfit$cv.maxsteps
    CV.trace <- CVBOX$cvfit$cv.trace
    CV.boxind <- CVBOX$cvfit$cv.boxind
    CV.boxcut <- CVBOX$cvfit$cv.boxcut
    CV.support <- CVBOX$cvfit$cv.stats$cv.support
    CV.size <- CVBOX$cvfit$cv.stats$cv.size
    CV.bmp.lhr <- CVBOX$cvfit$cv.stats$cv.bmp.lhr
    CV.bmp.lrt <- CVBOX$cvfit$cv.stats$cv.bmp.lrt
    CV.bmp.cer <- CVBOX$cvfit$cv.stats$cv.bmp.cer
    CV.grp.lhr <- CVBOX$cvfit$cv.stats$cv.grp.lhr
    CV.grp.lrt <- CVBOX$cvfit$cv.stats$cv.grp.lrt
    CV.grp.cer <- CVBOX$cvfit$cv.stats$cv.grp.cer
    CV.time.bar <- CVBOX$cvfit$cv.stats$cv.time.bar
    CV.prob.bar <- CVBOX$cvfit$cv.stats$cv.prob.bar
    CV.max.time.bar <- CVBOX$cvfit$cv.stats$cv.max.time.bar
    CV.min.prob.bar <- CVBOX$cvfit$cv.stats$cv.min.prob.bar
    
  }
  
  return(list("cv.maxsteps"=CV.maxsteps,
              "cv.trace"=CV.trace,
              "cv.boxind"=CV.boxind,
              "cv.boxcut"=CV.boxcut,
              "cv.support"=CV.support,
              "cv.size"=CV.size,
              "cv.bmp.lhr"=CV.bmp.lhr,
              "cv.bmp.lrt"=CV.bmp.lrt,
              "cv.bmp.cer"=CV.bmp.cer,
              "cv.grp.lhr"=CV.grp.lhr,
              "cv.grp.lrt"=CV.grp.lrt,
              "cv.grp.cer"=CV.grp.cer,
              "cv.time.bar"=CV.time.bar,
              "cv.prob.bar"=CV.prob.bar,
              "cv.max.time.bar"=CV.max.time.bar,
              "cv.min.prob.bar"=CV.min.prob.bar,
              "success"=success,
              "seed"=seed))
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    cv.pval (X,
#                             y,
#                             delta,
#                             K,
#                             A,
#                             pv,
#                             cv,
#                             cvtype,
#                             cvarg,
#                             groups,
#                             decimals,
#                             varsign,
#                             initcutpts,
#                             obs.stat,
#                             parallel,
#                             conf,
#                             verbose,
#                             seed)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

cv.pval <- function(X,
                    y,
                    delta,
                    K,
                    A,
                    pv,
                    cv,
                    cvtype,
                    cvarg,
                    groups,
                    decimals,
                    varsign,
                    initcutpts,
                    obs.stat,
                    parallel,
                    conf,
                    verbose,
                    seed) {
  
  if (pv) {
    
    cat("Computation of p-values ... \n")
    
    # Parsing and evaluating 'cvarg' parameters to evaluate 'L' and 'peelcriterion'
    alpha <- NULL
    beta <- NULL
    L <- NULL
    peelcriterion <- NULL
    cvcriterion <- NULL
    eval(parse( text=unlist(strsplit(x=cvarg, split=",")) ))
    
    seed <- seed[1]
    
    # Computation of log-rank permutation p-values using the permutation distribution for the null distribution
    nperm <- 1:A
    
    if (!parallel) {
      
      if (!is.null(seed)) {
        seed <- (0:(A-1)) + seed
      }
      
      null.obj <- cv.null(X=X,
                          y=y,
                          delta=delta,
                          nperm=nperm,
                          K=K,
                          cv=cv,
                          cvtype=cvtype,
                          cvarg=cvarg,
                          groups=groups,
                          peelcriterion=peelcriterion,
                          decimals=decimals,
                          varsign=varsign,
                          initcutpts=initcutpts,
                          parallel=parallel,
                          verbose=verbose,
                          seed=seed)
      
    } else {
      
      # Parallel back-end registration
      if (conf$type == "SOCKET") {
        cluster <- parallel::makeCluster(spec=conf$spec,
                                         type="PSOCK",
                                         homogeneous=conf$homo,
                                         outfile=conf$outfile,
                                         verbose=conf$verbose)
      } else if (conf$type == "MPI") {
        cluster <- parallel::makeCluster(spec=conf$spec,
                                         type="MPI",
                                         homogeneous=conf$homo,
                                         outfile=conf$outfile,
                                         verbose=conf$verbose)
      } else {
        stop("Unrecognized cluster type: you must specify a \"SOCKET\" or \"MPI\" cluster type. Exiting ... \n\n")
      }
      
      # Parallelization
      parallel::clusterSetRNGStream(cl=cluster, iseed=seed)
      null.obj <- parallel::clusterApply(cl=cluster, x=nperm, fun=cv.null,
                                         X=X,
                                         y=y,
                                         delta=delta,
                                         K=K,
                                         cv=cv,
                                         cvtype=cvtype,
                                         cvarg=cvarg,
                                         groups=groups,
                                         peelcriterion=peelcriterion,
                                         decimals=decimals,
                                         varsign=varsign,
                                         initcutpts=initcutpts,
                                         parallel=parallel,
                                         verbose=verbose,
                                         seed=NULL)
      
      # Stopping the cluster and cleaning all MPI states
      parallel::stopCluster(cl=cluster)
      
    }
    
    # Collect the results from all replicates
    null.stat <- t(list2mat(list=null.obj, coltrunc=L, fill=NA))
    
    pval <- numeric(L)
    if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
      for (l in 1:L) {
        pval[l] <- mean((null.stat[l,] >= obs.stat[l]), na.rm=TRUE)
      }
    } else if (peelcriterion == "grp") {
      for (l in 1:L) {
        pval[l] <- mean((abs(null.stat[l,]) >= obs.stat[l]), na.rm=TRUE)
      }
    } else {
      stop("Invalid peeling criterion. Exiting ...\n\n")
    }
    pval <- round(pval, digits=floor(log(base=10, A)))
    names(pval) <- paste("step", 0:(L-1), sep="")
    return(list("pval"=pval, "seed"=seed))
    
  } else {
    
    cat("No computation of p-values. \n")
    return(NULL)
    
  }
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                   cv.null (X,
#                            y,
#                            delta,
#                            nperm,
#                            K,
#                            cv,
#                            cvtype,
#                            cvarg,
#                            groups,
#                            peelcriterion,
#                            decimals,
#                            varsign,
#                            initcutpts,
#                            parallel,
#                            verbose,
#                            seed)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

cv.null <- function(X,
                    y,
                    delta,
                    nperm,
                    K,
                    cv,
                    cvtype,
                    cvarg,
                    groups,
                    peelcriterion,
                    decimals,
                    varsign,
                    initcutpts,
                    parallel,
                    verbose,
                    seed) {
  
  if (!parallel) {
    
    A <- length(nperm)
    n <- nrow(X)
    null.stat <- vector(mode="list", length=A)
    
    a <- 1
    while (a <= A) {
      if (verbose) {
        cat("Permutation : ", a, "\n")
        cat("seed : ", seed[a], "\n", sep="")
      }
      perm.ind <- sample(x = 1:n, size = n, replace = FALSE, prob = NULL)
      perm.y <- y[perm.ind]
      perm.delta <- delta[perm.ind]
      if (cvtype == "averaged") {
        peelobj <- cv.ave.box(X=X,
                              y=perm.y,
                              delta=perm.delta,
                              K=K,
                              cv=cv,
                              cvarg=cvarg,
                              varsign=varsign,
                              initcutpts=initcutpts,
                              groups=groups,
                              decimals=decimals,
                              probval=NULL,
                              timeval=NULL,
                              verbose=verbose,
                              seed=seed[a])
        if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
          null.stat[[a]] <- peelobj$cvfit$cv.stats$cv.bmp.lrt
        } else if (peelcriterion == "grp") {
          null.stat[[a]] <- peelobj$cvfit$cv.stats$cv.grp.lhr
        } else {
          stop("Invalid peeling criterion. Exiting ...\n\n")
        }
      } else if (cvtype == "combined") {
        peelobj <- cv.comb.box(X=X,
                               y=perm.y,
                               delta=perm.delta,
                               K=K,
                               cv=cv,
                               cvarg=cvarg,
                               varsign=varsign,
                               initcutpts=initcutpts,
                               groups=groups,
                               decimals=decimals,
                               probval=NULL,
                               timeval=NULL,
                               verbose=verbose,
                               seed=seed[a])
        if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
          null.stat[[a]] <- peelobj$cvfit$cv.stats$cv.bmp.lrt
        } else if (peelcriterion == "grp") {
          null.stat[[a]] <- peelobj$cvfit$cv.stats$cv.grp.lhr
        } else {
          stop("Invalid peeling criterion. Exiting ...\n\n")
        }
      } else {
        stop("Invalid CV type option. Exiting ... \n\n")
      }
      a <- a + 1
    }
    
  } else {
    
    n <- nrow(X)
    perm.ind <- sample(x = 1:n, size = n, replace = FALSE, prob = NULL)
    perm.y <- y[perm.ind]
    perm.delta <- delta[perm.ind]
    if (cvtype == "averaged") {
      peelobj <- cv.ave.box(X=X,
                            y=perm.y,
                            delta=perm.delta,
                            K=K,
                            cv=cv,
                            cvarg=cvarg,
                            decimals=decimals,
                            varsign=varsign,
                            initcutpts=initcutpts,
                            groups=groups,
                            probval=NULL,
                            timeval=NULL,
                            verbose=verbose,
                            seed=seed)
      if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
        null.stat <- peelobj$cvfit$cv.stats$cv.bmp.lrt
      } else if (peelcriterion == "grp") {
        null.stat <- peelobj$cvfit$cv.stats$cv.grp.lhr
      } else {
        stop("Invalid peeling criterion. Exiting ...\n\n")
      }
    } else if (cvtype == "combined") {
      peelobj <- cv.comb.box(X=X,
                             y=perm.y,
                             delta=perm.delta,
                             K=K,
                             cv=cv,
                             cvarg=cvarg,
                             decimals=decimals,
                             varsign=varsign,
                             initcutpts=initcutpts,
                             groups=groups,
                             probval=NULL,
                             timeval=NULL,
                             verbose=verbose,
                             seed=seed)
      if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
        null.stat <- peelobj$cvfit$cv.stats$cv.bmp.lrt
      } else if (peelcriterion == "grp") {
        null.stat <- peelobj$cvfit$cv.stats$cv.grp.lhr
      } else {
        stop("Invalid peeling criterion. Exiting ...\n\n")
      }
    } else {
      stop("Invalid CV type option. Exiting ... \n\n")
    }
  }
  
  return(null.stat)
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                   cv.ave.box (X,
#                               y,
#                               delta,
#                               K,
#                               cv,
#                               cvarg,
#                               varsign,
#                               initcutpts,
#                               groups,
#                               decimals,
#                               probval,
#                               timeval,
#                               verbose,
#                               seed)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

cv.ave.box <- function(X,
                       y,
                       delta,
                       K,
                       cv,
                       cvarg,
                       varsign,
                       initcutpts,
                       groups,
                       decimals,
                       probval,
                       timeval,
                       verbose,
                       seed) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Parsing and evaluating 'cvarg' parameters to evaluate 'beta'
  alpha <- NULL
  beta <- NULL
  L <- NULL
  peelcriterion <- NULL
  cvcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=cvarg, split=",")) ))
  
  # Creating argument `arg` of `cv.ave.peel` for PRSP parameters
  arg <- paste("alpha=", alpha,
               ",beta=", beta,
               ",L=", L,
               ",peelcriterion=\"", peelcriterion, "\"", sep="")
  
  folds <- cv.folds(y=delta, K=K, seed=seed)
  
  boxstat.list <- vector(mode="list", length=K)
  boxcut.list <- vector(mode="list", length=K)
  trace.list <- vector(mode="list", length=K)
  maxsteps <- numeric(K)
  ind <- numeric(0)

  for (k in 1:K) {
    if ((verbose) && (cv)) cat("Fold : ", k, "\n", sep="")
    # Initialize training and test data
    if (K == 1) {
      traindata <- testdata <- X[folds$perm[(folds$which == k)], , drop=FALSE]
      traintime <- testtime <- y[folds$perm[(folds$which == k)]]
      trainstatus <- teststatus <- delta[folds$perm[(folds$which == k)]]
      traingroups <- testgroups <- groups[folds$perm[(folds$which == k)]]
    } else {
      traindata <- X[folds$perm[(folds$which != k)], , drop=FALSE]
      traintime <- y[folds$perm[(folds$which != k)]]
      trainstatus <- delta[folds$perm[(folds$which != k)]]
      traingroups <- groups[folds$perm[(folds$which != k)]]
      testdata <- X[folds$perm[(folds$which == k)], , drop=FALSE]
      testtime <- y[folds$perm[(folds$which == k)]]
      teststatus <- delta[folds$perm[(folds$which == k)]]
      testgroups <- groups[folds$perm[(folds$which == k)]]
    }
    # Store the test set results from each fold
    peelobj <- tryCatch(expr={cv.ave.peel(traindata=traindata, trainstatus=trainstatus, 
                                          traintime=traintime, traingroups=traingroups,
                                          testdata=testdata, teststatus=teststatus, 
                                          testtime=testtime, testgroups=testgroups,
                                          varsign=varsign, initcutpts=initcutpts, 
                                          arg=arg, probval=probval, timeval=timeval)},
                        error=function(x){NULL})
    if (is.list(peelobj)) {
      maxsteps[k] <- peelobj$maxsteps
      boxstat.list[[k]] <- peelobj$boxstat
      boxcut.list[[k]] <- peelobj$boxcut
      trace.list[[k]] <- peelobj$trace
    } else {
      if (verbose) cat("CV fold dropped... \n")
      ind <- c(ind, k)
      maxsteps[k] <- NA
      boxstat.list[[k]] <- NA
      boxcut.list[[k]] <- NA
      trace.list[[k]] <- NA
    }
  }
  
  if (!is.empty(ind)) {
    maxsteps <- maxsteps[-ind]
    boxstat.list <- boxstat.list[-ind]
    boxcut.list <- boxcut.list[-ind]
    trace.list <- trace.list[-ind]
  }

  # Cross-validated maximum peeling length from all folds
  CV.Lm <- min(maxsteps, na.rm=TRUE)
  
  # Compute the averaged box for each step from all the folds
  # Truncate the cross-validated quantities from all folds to the same cross-validated length
  for (k in 1:K) {
    if (!is.matrix(boxcut.list[[k]])) {
      boxcut.list[[k]] <- matrix(data=NA, nrow=CV.Lm, ncol=p)
    } else {
      boxcut.list[[k]] <- boxcut.list[[k]][1:CV.Lm,,drop=FALSE]
    }
  }  
  # Each row of the final matrix signifies a step
  CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), colnames(X)))
  for (l in 1:CV.Lm) {
    summincut <- matrix(NA, K, p)
    for (k in 1:K) {
      summincut[k,] <- boxcut.list[[k]][l,]
    }
    CV.boxcut[l, ] <- colMeans(summincut, na.rm=TRUE)
  }
  rownames(CV.boxcut) <- paste("step", 0:(CV.Lm-1), sep="")
  colnames(CV.boxcut) <- colnames(X)
  
  # Get the box membership indicator vector of all observations for each step from all the folds
  # Based on the average box over the folds
  CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
  for (l in 1:CV.Lm) {
    boxcut <- CV.boxcut[l, ] * varsign
    X.cut <- t(t(X) * varsign)
    X.ind <- t(t(X.cut) >= boxcut)
    CV.boxind[l,] <- (rowMeans(X.ind, na.rm = TRUE) == 1)  # Set as TRUE which observations are inside the box boudaries for all axes directions
  }
  rownames(CV.boxind) <- paste("step", 0:(CV.Lm-1), sep="")
  colnames(CV.boxind) <- rownames(X)
  
  # Get the adjusted cross-validated maximum peeling length, thresholded by minimal box support
  CV.Lm <- max(which(apply(CV.boxind, 1, function(x) {length(which(x))/n >= max(1/n, beta)})), na.rm=TRUE)
  
  # Get the adjusted test box definition and membership indicator matrices for each step
  CV.boxind <- CV.boxind[1:CV.Lm,,drop=FALSE]
  CV.boxcut <- CV.boxcut[1:CV.Lm,,drop=FALSE]
  
  # Get the box sample size and support based on `CV.boxind`
  CV.size <- apply(CV.boxind, 1, sum, na.rm=TRUE)
  names(CV.size) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.support <- CV.size/n
  names(CV.support) <- paste("step", 0:(CV.Lm-1), sep="")
  
  # Get the adjusted cross-validated quantities from all folds to the same cross-validated length
  for (k in 1:K) {
    if (!is.list(boxstat.list[[k]])) {
      boxstat.list[[k]] <- rep(list(rep(NA,10)),CV.Lm)
    } else {
      boxstat.list[[k]] <- boxstat.list[[k]][1:CV.Lm]
    }
  }
  
  # Get the variable traces
  # Variable traces are first stacked and truncated in a matrix where folds are by rows and steps by columns
  CV.trace <- lapply.mat(X=trace.list,
                         coltrunc=CV.Lm,
                         fill=NA,
                         MARGIN=2,
                         FUN=function(x) {
                           vote <- table(x, useNA="no")
                           w <- as.numeric(names(which.max(vote)))
                           if(is.empty(w)) {
                             return(NA)
                           } else {
                             return(w)
                           }
                         })
  names(CV.trace) <- paste("step", 0:(CV.Lm-1), sep="")
  
  # Box peeling rules for each step
  CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), colnames(X))))
  for (j in 1:p) {
    if (varsign[j] > 0) {
      ss <- ">="
    } else {
      ss <- "<="
    }
    CV.rules[, j] <- paste(colnames(X)[j], ss, format(x=CV.boxcut[, j], digits=decimals, nsmall=decimals), sep="")
  }
  
  # Compute the averaged box statistics for each step from all the folds
  # Each entry or row signifies a step
  CV.bmp.lhr <- rep(x=NA, times=CV.Lm)
  names(CV.bmp.lhr) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.bmp.lrt <- rep(x=NA, times=CV.Lm)
  names(CV.bmp.lrt) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.bmp.cer <- rep(x=NA, times=CV.Lm)
  names(CV.bmp.cer) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.grp.lhr <- rep(x=NA, times=CV.Lm)
  names(CV.grp.lhr) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.grp.lrt <- rep(x=NA, times=CV.Lm)
  names(CV.grp.lrt) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.grp.cer <- rep(x=NA, times=CV.Lm)
  names(CV.grp.cer) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.time.bar <- rep(x=NA, times=CV.Lm)
  names(CV.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.prob.bar <- rep(x=NA, times=CV.Lm)
  names(CV.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.max.time.bar <- rep(x=NA, times=CV.Lm)
  names(CV.max.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.min.prob.bar <- rep(x=NA, times=CV.Lm)
  names(CV.min.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  for (l in 1:CV.Lm) {
    sumbmplhr <- rep(x=NA, times=K)
    sumbmplrt <- rep(x=NA, times=K)
    sumbmpcer <- rep(x=NA, times=K)
    sumgrplhr <- rep(x=NA, times=K)
    sumgrplrt <- rep(x=NA, times=K)
    sumgrpcer <- rep(x=NA, times=K)
    sumtime <- rep(x=NA, times=K)
    sumprob <- rep(x=NA, times=K)
    summaxtime <- rep(x=NA, times=K)
    summinprob <- rep(x=NA, times=K)
    for (k in 1:K) {
      outbounds <- boxstat.list[[k]][[l]]
      if (!is.null(outbounds)) {
        sumbmplhr[k] <- boxstat.list[[k]][[l]][[1]]
        sumbmplrt[k] <- boxstat.list[[k]][[l]][[2]]
        sumbmpcer[k] <- boxstat.list[[k]][[l]][[3]]
        sumgrplhr[k] <- boxstat.list[[k]][[l]][[4]]
        sumgrplrt[k] <- boxstat.list[[k]][[l]][[5]]
        sumgrpcer[k] <- boxstat.list[[k]][[l]][[6]]
        sumtime[k] <- boxstat.list[[k]][[l]][[7]]
        sumprob[k] <- boxstat.list[[k]][[l]][[8]]
        summaxtime[k] <- boxstat.list[[k]][[l]][[9]]
        summinprob[k] <- boxstat.list[[k]][[l]][[10]]
      }
    }
    CV.bmp.lhr[l] <- mean(sumbmplhr, na.rm=TRUE)
    CV.bmp.lrt[l] <- mean(sumbmplrt, na.rm=TRUE)
    CV.bmp.cer[l] <- mean(sumbmpcer, na.rm=TRUE)
    CV.grp.lhr[l] <- mean(sumgrplhr, na.rm=TRUE)
    CV.grp.lrt[l] <- mean(sumgrplrt, na.rm=TRUE)
    CV.grp.cer[l] <- mean(sumgrpcer, na.rm=TRUE)
    CV.time.bar[l] <- mean(sumtime, na.rm=TRUE)
    CV.prob.bar[l] <- mean(sumprob, na.rm=TRUE)
    CV.max.time.bar[l] <- mean(summaxtime, na.rm=TRUE)
    CV.min.prob.bar[l] <- mean(summinprob, na.rm=TRUE)
  }
  CV.bmp.lhr[is.nan(CV.bmp.lhr)] <- NA
  CV.bmp.lrt[is.nan(CV.bmp.lrt)] <- NA
  CV.bmp.cer[is.nan(CV.bmp.cer)] <- NA
  CV.grp.lhr[is.nan(CV.grp.lhr)] <- NA
  CV.grp.lrt[is.nan(CV.grp.lrt)] <- NA
  CV.grp.cer[is.nan(CV.grp.cer)] <- NA
  CV.time.bar[is.nan(CV.time.bar)] <- NA
  CV.prob.bar[is.nan(CV.prob.bar)] <- NA
  CV.max.time.bar[is.nan(CV.max.time.bar)] <- NA
  CV.min.prob.bar[is.nan(CV.min.prob.bar)] <- NA
  CV.maxsteps <- CV.Lm
  
  # Formatting the results depending on successful cross-validated PRSP algorithm
  if (cvcriterion == "lhr") {
    if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
      if (all(is.na(CV.bmp.lhr)) || all(is.nan(CV.bmp.lhr))) {
        success <- FALSE
        CV.maxsteps <- NA
        CV.support <- rep(x=NA, times=CV.Lm)
        CV.size <- rep(x=NA, times=CV.Lm)
        CV.bmp.lhr <- rep(x=NA, times=CV.Lm)
        CV.bmp.lrt <- rep(x=NA, times=CV.Lm)
        CV.bmp.cer <- rep(x=NA, times=CV.Lm)
        CV.grp.lhr <- rep(x=NA, times=CV.Lm)
        CV.grp.lrt <- rep(x=NA, times=CV.Lm)
        CV.grp.cer <- rep(x=NA, times=CV.Lm)
        CV.time.bar <- rep(x=NA, times=CV.Lm)
        CV.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.max.time.bar <- rep(x=NA, times=CV.Lm)
        CV.min.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p)
        CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
        CV.trace <- rep(x=NA, times=CV.Lm)
        CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p))
      } else {
        # Cross-validated optimal length from all folds by maximization of the LHR (between inbox and outbox test samples)
        success <- TRUE
      }
    } else if (peelcriterion == "grp") {
      if (all(is.na(CV.grp.lhr)) || all(is.nan(CV.grp.lhr))) {
        success <- FALSE
        CV.maxsteps <- NA
        CV.support <- rep(x=NA, times=CV.Lm)
        CV.size <- rep(x=NA, times=CV.Lm)
        CV.bmp.lhr <- rep(x=NA, times=CV.Lm)
        CV.bmp.lrt <- rep(x=NA, times=CV.Lm)
        CV.bmp.cer <- rep(x=NA, times=CV.Lm)
        CV.grp.lhr <- rep(x=NA, times=CV.Lm)
        CV.grp.lrt <- rep(x=NA, times=CV.Lm)
        CV.grp.cer <- rep(x=NA, times=CV.Lm)
        CV.time.bar <- rep(x=NA, times=CV.Lm)
        CV.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.max.time.bar <- rep(x=NA, times=CV.Lm)
        CV.min.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p)
        CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
        CV.trace <- rep(x=NA, times=CV.Lm)
        CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p))
      } else {
        # Cross-validated optimal length from all folds by minimization of the CER (between predicted and observed inbox test samples survival times)
        success <- TRUE
      }
    } else {
      stop("Invalid peeling criterion. Exiting ...\n\n")
    }
  } else if (cvcriterion == "lrt") {
    if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
      if (all(is.na(CV.bmp.lrt)) || all(is.nan(CV.bmp.lrt))) {
        success <- FALSE
        CV.maxsteps <- NA
        CV.support <- rep(x=NA, times=CV.Lm)
        CV.size <- rep(x=NA, times=CV.Lm)
        CV.bmp.lhr <- rep(x=NA, times=CV.Lm)
        CV.bmp.lrt <- rep(x=NA, times=CV.Lm)
        CV.bmp.cer <- rep(x=NA, times=CV.Lm)
        CV.grp.lhr <- rep(x=NA, times=CV.Lm)
        CV.grp.lrt <- rep(x=NA, times=CV.Lm)
        CV.grp.cer <- rep(x=NA, times=CV.Lm)
        CV.time.bar <- rep(x=NA, times=CV.Lm)
        CV.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.max.time.bar <- rep(x=NA, times=CV.Lm)
        CV.min.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p)
        CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
        CV.trace <- rep(x=NA, times=CV.Lm)
        CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p))
      } else {
        # Cross-validated optimal length from all folds by maximization of the LRT (between inbox and outbox test samples)
        success <- TRUE
      }
    } else if (peelcriterion == "grp") {
      if (all(is.na(CV.grp.lrt)) || all(is.nan(CV.grp.lrt))) {
        success <- FALSE
        CV.maxsteps <- NA
        CV.support <- rep(x=NA, times=CV.Lm)
        CV.size <- rep(x=NA, times=CV.Lm)
        CV.bmp.lhr <- rep(x=NA, times=CV.Lm)
        CV.bmp.lrt <- rep(x=NA, times=CV.Lm)
        CV.bmp.cer <- rep(x=NA, times=CV.Lm)
        CV.grp.lhr <- rep(x=NA, times=CV.Lm)
        CV.grp.lrt <- rep(x=NA, times=CV.Lm)
        CV.grp.cer <- rep(x=NA, times=CV.Lm)
        CV.time.bar <- rep(x=NA, times=CV.Lm)
        CV.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.max.time.bar <- rep(x=NA, times=CV.Lm)
        CV.min.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p)
        CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
        CV.trace <- rep(x=NA, times=CV.Lm)
        CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p))
      } else {
        # Cross-validated optimal length from all folds by minimization of the CER (between predicted and observed inbox test samples survival times)
        success <- TRUE
      }
    } else {
      stop("Invalid peeling criterion. Exiting ...\n\n")
    }
  } else if (cvcriterion == "cer") {
    if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
      if (all(is.na(CV.bmp.cer)) || all(is.nan(CV.bmp.cer))) {
        success <- FALSE
        CV.maxsteps <- NA
        CV.support <- rep(x=NA, times=CV.Lm)
        CV.size <- rep(x=NA, times=CV.Lm)
        CV.bmp.lhr <- rep(x=NA, times=CV.Lm)
        CV.bmp.lrt <- rep(x=NA, times=CV.Lm)
        CV.bmp.cer <- rep(x=NA, times=CV.Lm)
        CV.grp.lhr <- rep(x=NA, times=CV.Lm)
        CV.grp.lrt <- rep(x=NA, times=CV.Lm)
        CV.grp.cer <- rep(x=NA, times=CV.Lm)
        CV.time.bar <- rep(x=NA, times=CV.Lm)
        CV.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.max.time.bar <- rep(x=NA, times=CV.Lm)
        CV.min.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p)
        CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
        CV.trace <- rep(x=NA, times=CV.Lm)
        CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p))
      } else {
        # Cross-validated optimal length from all folds by minimization of the CER (between predicted and observed inbox test samples survival times)
        success <- TRUE
      }
    } else if (peelcriterion == "grp") {
      if (all(is.na(CV.grp.cer)) || all(is.nan(CV.grp.cer))) {
        success <- FALSE
        CV.maxsteps <- NA
        CV.support <- rep(x=NA, times=CV.Lm)
        CV.size <- rep(x=NA, times=CV.Lm)
        CV.bmp.lhr <- rep(x=NA, times=CV.Lm)
        CV.bmp.lrt <- rep(x=NA, times=CV.Lm)
        CV.bmp.cer <- rep(x=NA, times=CV.Lm)
        CV.grp.lhr <- rep(x=NA, times=CV.Lm)
        CV.grp.lrt <- rep(x=NA, times=CV.Lm)
        CV.grp.cer <- rep(x=NA, times=CV.Lm)
        CV.time.bar <- rep(x=NA, times=CV.Lm)
        CV.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.max.time.bar <- rep(x=NA, times=CV.Lm)
        CV.min.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p)
        CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
        CV.trace <- rep(x=NA, times=CV.Lm)
        CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p))
      } else {
        # Cross-validated optimal length from all folds by minimization of the CER (between predicted and observed inbox test samples survival times)
        success <- TRUE
      }
    } else {
      stop("Invalid peeling criterion. Exiting ...\n\n")
    }
  } else {
    stop("Invalid CV criterion option. Exiting ... \n\n")
  }
  
  # Box statistics for each step
  CV.stats <-  data.frame("cv.support"=CV.support,
                          "cv.size"=CV.size,
                          "cv.bmp.lhr"=CV.bmp.lhr,
                          "cv.bmp.lrt"=CV.bmp.lrt,
                          "cv.bmp.cer"=CV.bmp.cer,
                          "cv.grp.lhr"=CV.grp.lhr,
                          "cv.grp.lrt"=CV.grp.lrt,
                          "cv.grp.cer"=CV.grp.cer,
                          "cv.time.bar"=CV.time.bar,
                          "cv.prob.bar"=CV.prob.bar,
                          "cv.max.time.bar"=CV.max.time.bar,
                          "cv.min.prob.bar"=CV.min.prob.bar)
  rownames(CV.stats) <- paste("step", 0:(CV.Lm-1), sep="")
  
  # Create the return object 'CV.fit'
  CV.fit <- list("cv.maxsteps"=CV.maxsteps,
                 "cv.boxcut"=CV.boxcut,
                 "cv.boxind"=CV.boxind,
                 "cv.trace"=CV.trace,
                 "cv.rules"=CV.rules,
                 "cv.stats"=CV.stats)
  
  return(list("X"=X,
              "y"=y,
              "delta"=delta,
              "cvfit"=CV.fit,
              "success"=success,
              "seed"=seed))
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                   cv.comb.box (X,
#                                y,
#                                delta,
#                                K,
#                                cv,
#                                cvarg,
#                                varsign,
#                                initcutpts,
#                                groups,
#                                probval,
#                                timeval,
#                                decimals,
#                                verbose,
#                                seed)
#
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

cv.comb.box <- function(X,
                        y,
                        delta,
                        K,
                        cv,
                        cvarg,
                        varsign,
                        initcutpts,
                        groups,
                        decimals,
                        probval,
                        timeval,
                        verbose,
                        seed) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Parsing and evaluating 'cvarg' parameters to evaluate 'beta' and 'peelcriterion'
  alpha <- NULL
  beta <- NULL
  L <- NULL
  peelcriterion <- NULL
  cvcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=cvarg, split=",")) ))
  
  # Creating argument `arg` of `cv.comb.peel` for PRSP parameters
  arg <- paste("alpha=", alpha,
               ",beta=", beta,
               ",L=", L,
               ",peelcriterion=\"", peelcriterion, "\"", sep="")
  
  folds <- cv.folds(y=delta, K=K, seed=seed)
  ord <- folds$foldkey
  
  times.list <- vector(mode="list", length=K)
  status.list <- vector(mode="list", length=K)
  boxind.list <- vector(mode="list", length=K)
  boxcut.list <- vector(mode="list", length=K)
  trace.list <- vector(mode="list", length=K)
  maxsteps <- numeric(K)
  ind <- numeric(0)
 
  for (k in 1:K) {
    if ((verbose) && (cv)) cat("Fold : ", k, "\n", sep="")
    if (K == 1) {
      traindata <- testdata <- X[folds$perm[(folds$which == k)], , drop=FALSE]
      traintime <- testtime <- y[folds$perm[(folds$which == k)]]
      trainstatus <- teststatus <- delta[folds$perm[(folds$which == k)]]
      traingroups <- testgroups <- groups[folds$perm[(folds$which == k)]]
    } else {
      traindata <- X[folds$perm[(folds$which != k)], , drop=FALSE]
      traintime <- y[folds$perm[(folds$which != k)]]
      trainstatus <- delta[folds$perm[(folds$which != k)]]
      traingroups <- groups[folds$perm[(folds$which != k)]]
      testdata <- X[folds$perm[(folds$which == k)], , drop=FALSE]
      testtime <- y[folds$perm[(folds$which == k)]]
      teststatus <- delta[folds$perm[(folds$which == k)]]
      testgroups <- groups[folds$perm[(folds$which == k)]]
    }
    peelobj <- tryCatch(expr={cv.comb.peel(traindata=traindata, trainstatus=trainstatus, traintime=traintime,
                                           testdata=testdata, teststatus=teststatus, testtime=testtime,
                                           varsign=varsign, initcutpts=initcutpts, 
                                           arg=arg, traingroups=traingroups, testgroups=testgroups)},
                        error=function(x){NULL})
    # Store the test set data from each fold (Note: there is a new order of observations from each fold)
    times.list[[k]] <- testtime
    status.list[[k]] <- teststatus
    # Store the test set results from each fold
    if (is.list(peelobj)) {
      maxsteps[k] <- peelobj$maxsteps
      boxind.list[[k]] <- peelobj$boxind
      trace.list[[k]] <- peelobj$trace
    } else {
      if (verbose) cat("CV fold dropped... \n")
      ind <- c(ind, k)
      maxsteps[k] <- NA
      boxind.list[[k]] <- matrix(data=NA, nrow=1, ncol=length(which(folds$which == k)))
      trace.list[[k]] <- NA
    }
  }
  if (!is.empty(ind)) {
    maxsteps <- maxsteps[-ind]
    trace.list <- trace.list[-ind]
  }

  # Cross-validated maximum peeling length from all folds
  CV.Lm <- min(maxsteps, na.rm=TRUE)
   
  # Get the test box membership indicator vector of all observations for each step from all the folds
  # Based on the combined membership indicator vectors over the folds
  # Re-ordered by initial order of observations
  CV.boxind <- cbindlist(boxind.list, trunc=CV.Lm)[,ord,drop=FALSE]
  rownames(CV.boxind) <- paste("step", 0:(CV.Lm-1), sep="")
  colnames(CV.boxind) <- rownames(X)
  
  # Get the combined boxcut (truncated to the same cross-validated length) for each step from all the folds
  # using the box circumscribing the combined test set samples over all the folds
  CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), colnames(X)))
  for (l in 1:CV.Lm) {
    for (j in 1:p) {
      if (varsign[j] > 0) {
        CV.boxcut[l,j] <- min(X[CV.boxind[l,],j], na.rm=TRUE)
      } else {
        CV.boxcut[l,j] <- max(X[CV.boxind[l,],j], na.rm=TRUE)
      }
    }
  }
  
  # Get the box membership indicator vector of all observations (including the dropped ones in any fold) for each step from all the folds
  # Based on the combined box over the folds
  CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
  for (l in 1:CV.Lm) {
    boxcut <- CV.boxcut[l, ] * varsign
    X.cut <- t(t(X) * varsign)
    X.ind <- t(t(X.cut) >= boxcut)
    CV.boxind[l,] <- (rowMeans(X.ind, na.rm = TRUE) == 1)  # Set as TRUE which observations are inside the box boudaries for all axes directions
  }
  rownames(CV.boxind) <- paste("step", 0:(CV.Lm-1), sep="")
  colnames(CV.boxind) <- rownames(X)

  # Get the adjusted cross-validated maximum peeling length, thresholded by minimal box support
  CV.Lm <- max(which(apply(CV.boxind, 1, function(x) {length(which(x))/n >= max(1/n, beta)})), na.rm=TRUE)
  
  # Get the adjusted test box definition and membership indicator matrices for each step
  CV.boxind <- CV.boxind[1:CV.Lm,,drop=FALSE]
  CV.boxcut <- CV.boxcut[1:CV.Lm,,drop=FALSE]
 
  # Get the box sample size and support based on `CV.boxind`
  CV.size <- apply(CV.boxind, 1, sum, na.rm=TRUE)
  names(CV.size) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.support <- CV.size/n
  names(CV.support) <- paste("step", 0:(CV.Lm-1), sep="")
  
  # Concatenates the observations of test times and status from all folds
  # Re-ordered by initial order of observations
  CV.times <- unlist(times.list)[ord]
  CV.status <- unlist(status.list)[ord]

  # Get the variable traces
  # Variable traces are first stacked and truncated in a matrix where folds are by rows and steps by columns
  CV.trace <- lapply.mat(X=trace.list,
                         coltrunc=CV.Lm,
                         fill=NA,
                         MARGIN=2,
                         FUN=function(x) {
                           vote <- table(x, useNA="no")
                           w <- as.numeric(names(which.max(vote)))
                           if(is.empty(w)) {
                             return(NA)
                           } else {
                             return(w)
                           }
                         })
  names(CV.trace) <- paste("step", 0:(CV.Lm-1), sep="")

  # Box peeling rules for each step
  CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), colnames(X))))
  for (j in 1:p) {
    if (varsign[j] > 0) {
      ss <- ">="
    } else {
      ss <- "<="
    }
    CV.rules[, j] <- paste(colnames(X)[j], ss, format(x=CV.boxcut[, j], digits=decimals, nsmall=decimals), sep="")
  }
  
  # Compute the combined test box statistics from all folds for all steps, each entry or row signifies a step
  CV.bmp.lhr <- rep(x=NA, times=CV.Lm)
  names(CV.bmp.lhr) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.bmp.lrt <- rep(x=NA, times=CV.Lm)
  names(CV.bmp.lrt) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.bmp.cer <- rep(x=NA, times=CV.Lm)
  names(CV.bmp.cer) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.grp.lhr <- rep(x=NA, times=CV.Lm)
  names(CV.grp.lhr) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.grp.lrt <- rep(x=NA, times=CV.Lm)
  names(CV.grp.lrt) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.grp.cer <- rep(x=NA, times=CV.Lm)
  names(CV.grp.cer) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.time.bar <- rep(x=NA, times=CV.Lm)
  names(CV.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.prob.bar <- rep(x=NA, times=CV.Lm)
  names(CV.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.max.time.bar <- rep(x=NA, times=CV.Lm)
  names(CV.max.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.min.prob.bar <- rep(x=NA, times=CV.Lm)
  names(CV.min.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  timemat <- matrix(NA, nrow=CV.Lm, ncol=n)
  probmat <- matrix(NA, nrow=CV.Lm, ncol=n)
  ind.rem <- numeric(0)
  grp1 <- 1*(groups == levels(groups)[1])
  wg <- which(grp1 == 1)
  for (l in 1:CV.Lm) {
    box1 <- 1*CV.boxind[l,]
    wb <- which(box1 == 1)
    if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
      if ((sum(box1, na.rm=TRUE) == length(box1[!is.na(box1)])) || (sum(box1, na.rm=TRUE) == 0)) {
        surv.formula <- (survival::Surv(CV.times, CV.status) ~ 1 + box1)
        surv.fit <- survival::survfit(surv.formula, na.action="na.omit")
        timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
        probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
        CV.bmp.lhr[l] <- 0
        CV.bmp.lrt[l] <- 0
        CV.bmp.cer[l] <- 1
        CV.grp.lhr[l] <- NA
        CV.grp.lrt[l] <- NA
        CV.grp.cer[l] <- NA
      } else if ((sum(box1, na.rm=TRUE) != length(box1[!is.na(box1)])) && (sum(box1, na.rm=TRUE) != 0)) {
        surv.formula <- (survival::Surv(CV.times, CV.status) ~ 1 + box1)
        surv.fit <- survival::survfit(surv.formula, na.action="na.omit")
        timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
        probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
        coxobj <- survival::coxph(surv.formula, na.action="na.omit", iter.max=20, timefix=TRUE, method="efron")
        CV.bmp.lhr[l] <- coxobj$coef
        CV.bmp.lrt[l] <- survival::survdiff(surv.formula, na.action="na.omit", rho=0)$chisq
        predobj <- predict(object=coxobj, type="lp", reference="sample", na.action="na.omit")
        CV.bmp.cer[l] <- Hmisc::rcorr.cens(x=predobj, S=survival::Surv(CV.times, CV.status))['C Index']
        CV.grp.lhr[l] <- NA
        CV.grp.lrt[l] <- NA
        CV.grp.cer[l] <- NA
      } else {
        timemat[l, ] <- NA
        probmat[l, ] <- NA
        CV.bmp.lhr[l] <- NA
        CV.bmp.lrt[l] <- NA
        CV.bmp.cer[l] <- NA
        CV.grp.lhr[l] <- NA
        CV.grp.lrt[l] <- NA
        CV.grp.cer[l] <- NA
        ind.rem <- c(ind.rem, l)
      }
    } else if (peelcriterion == "grp") {
      if ((sum((grp1 * box1), na.rm=TRUE) == length((grp1 * box1)[!is.na(grp1 * box1)])) ||
          (sum((grp1 * box1), na.rm=TRUE) == 0)) {
        timemat[l, ] <- NA
        probmat[l, ] <- NA
        CV.bmp.lhr[l] <- NA
        CV.bmp.lrt[l] <- NA
        CV.bmp.cer[l] <- NA
        CV.grp.lhr[l] <- 0
        CV.grp.lrt[l] <- NA
        CV.grp.cer[l] <- 1
      } else if ((sum((grp1 * box1), na.rm=TRUE) != length((grp1 * box1)[!is.na((grp1 * box1))])) &&
                 (sum((grp1 * box1), na.rm=TRUE) != 0)) {
        coxobj <- survival::coxph(survival::Surv(CV.times, CV.status) ~ 1 + grp1 * box1, 
                                  data=list("grp1"=grp1, "box1"=box1),
                                  na.action="na.omit", 
                                  iter.max=20,
                                  timefix=TRUE,
                                  method="efron")
        predobj <- predict(object=coxobj, type="lp", reference="sample", na.action="na.omit")
        timemat[l, ] <- NA
        probmat[l, ] <- NA
        CV.bmp.lhr[l] <- NA
        CV.bmp.lrt[l] <- NA
        CV.bmp.cer[l] <- NA
        CV.grp.lhr[l] <- coxobj$coef["grp1:box1"]
        CV.grp.lrt[l] <- NA
        CV.grp.cer[l] <- Hmisc::rcorr.cens(x=predobj, S=survival::Surv(CV.times, CV.status))['C Index']
      } else {
        timemat[l, ] <- NA
        probmat[l, ] <- NA
        CV.bmp.lhr[l] <- NA
        CV.bmp.lrt[l] <- NA
        CV.bmp.cer[l] <- NA
        CV.grp.lhr[l] <- NA
        CV.grp.lrt[l] <- NA
        CV.grp.cer[l] <- NA
        ind.rem <- c(ind.rem, l)
      }
    } else {
      stop("Invalid peeling criterion. Exiting ...\n\n")
    }
  }
  if ((length(ind.rem) == CV.Lm) || all(is.na(timemat)) || all(is.na(probmat))) {
    CV.time.bar <- rep(x=NA, times=CV.Lm)
    CV.prob.bar <- rep(x=NA, times=CV.Lm)
    CV.max.time.bar <- rep(x=NA, times=CV.Lm)
    CV.min.prob.bar <- rep(x=NA, times=CV.Lm)
  } else {
    endobj <- endpoints(ind=ind.rem, timemat=timemat, probmat=probmat, timeval=timeval, probval=probval)
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
  }
  CV.maxsteps <- CV.Lm
  
  # Formatting the results depending on successful cross-validated PRSP algorithm
  if (cvcriterion == "lhr") {
    if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
      if (all(is.na(CV.bmp.lhr)) || all(is.nan(CV.bmp.lhr))) {
        success <- FALSE
        CV.maxsteps <- NA
        CV.support <- rep(x=NA, times=CV.Lm)
        CV.size <- rep(x=NA, times=CV.Lm)
        CV.bmp.lhr <- rep(x=NA, times=CV.Lm)
        CV.bmp.lrt <- rep(x=NA, times=CV.Lm)
        CV.bmp.cer <- rep(x=NA, times=CV.Lm)
        CV.grp.lhr <- rep(x=NA, times=CV.Lm)
        CV.grp.lrt <- rep(x=NA, times=CV.Lm)
        CV.grp.cer <- rep(x=NA, times=CV.Lm)
        CV.time.bar <- rep(x=NA, times=CV.Lm)
        CV.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.max.time.bar <- rep(x=NA, times=CV.Lm)
        CV.min.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p)
        CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
        CV.trace <- rep(x=NA, times=CV.Lm)
        CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p))
      } else {
        # Cross-validated optimal length from all folds by maximization of the LHR (between inbox and outbox test samples)
        success <- TRUE
      }
    } else if (peelcriterion == "grp") {
      if (all(is.na(CV.grp.lhr)) || all(is.nan(CV.grp.lhr))) {
        success <- FALSE
        CV.maxsteps <- NA
        CV.support <- rep(x=NA, times=CV.Lm)
        CV.size <- rep(x=NA, times=CV.Lm)
        CV.bmp.lhr <- rep(x=NA, times=CV.Lm)
        CV.bmp.lrt <- rep(x=NA, times=CV.Lm)
        CV.bmp.cer <- rep(x=NA, times=CV.Lm)
        CV.grp.lhr <- rep(x=NA, times=CV.Lm)
        CV.grp.lrt <- rep(x=NA, times=CV.Lm)
        CV.grp.cer <- rep(x=NA, times=CV.Lm)
        CV.time.bar <- rep(x=NA, times=CV.Lm)
        CV.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.max.time.bar <- rep(x=NA, times=CV.Lm)
        CV.min.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p)
        CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
        CV.trace <- rep(x=NA, times=CV.Lm)
        CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p))
      } else {
        # Cross-validated optimal length from all folds by minimization of the CER (between predicted and observed inbox test samples survival times)
        success <- TRUE
      }
    } else {
      stop("Invalid peeling criterion. Exiting ...\n\n")
    }
  } else if (cvcriterion == "lrt") {
    if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
      if (all(is.na(CV.bmp.lrt)) || all(is.nan(CV.bmp.lrt))) {
        success <- FALSE
        CV.maxsteps <- NA
        CV.support <- rep(x=NA, times=CV.Lm)
        CV.size <- rep(x=NA, times=CV.Lm)
        CV.bmp.lhr <- rep(x=NA, times=CV.Lm)
        CV.bmp.lrt <- rep(x=NA, times=CV.Lm)
        CV.bmp.cer <- rep(x=NA, times=CV.Lm)
        CV.grp.lhr <- rep(x=NA, times=CV.Lm)
        CV.grp.lrt <- rep(x=NA, times=CV.Lm)
        CV.grp.cer <- rep(x=NA, times=CV.Lm)
        CV.time.bar <- rep(x=NA, times=CV.Lm)
        CV.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.max.time.bar <- rep(x=NA, times=CV.Lm)
        CV.min.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p)
        CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
        CV.trace <- rep(x=NA, times=CV.Lm)
        CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p))
      } else {
        # Cross-validated optimal length from all folds by maximization of the LRT (between inbox and outbox test samples)
        success <- TRUE
      }
    } else if (peelcriterion == "grp") {
      if (all(is.na(CV.grp.lrt)) || all(is.nan(CV.grp.lrt))) {
        success <- FALSE
        CV.maxsteps <- NA
        CV.support <- rep(x=NA, times=CV.Lm)
        CV.size <- rep(x=NA, times=CV.Lm)
        CV.bmp.lhr <- rep(x=NA, times=CV.Lm)
        CV.bmp.lrt <- rep(x=NA, times=CV.Lm)
        CV.bmp.cer <- rep(x=NA, times=CV.Lm)
        CV.grp.lhr <- rep(x=NA, times=CV.Lm)
        CV.grp.lrt <- rep(x=NA, times=CV.Lm)
        CV.grp.cer <- rep(x=NA, times=CV.Lm)
        CV.time.bar <- rep(x=NA, times=CV.Lm)
        CV.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.max.time.bar <- rep(x=NA, times=CV.Lm)
        CV.min.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p)
        CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
        CV.trace <- rep(x=NA, times=CV.Lm)
        CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p))
      } else {
        # Cross-validated optimal length from all folds by minimization of the CER (between predicted and observed inbox test samples survival times)
        success <- TRUE
      }
    } else {
      stop("Invalid peeling criterion. Exiting ...\n\n")
    }
  } else if (cvcriterion == "cer") {
    if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
      if (all(is.na(CV.bmp.cer)) || all(is.nan(CV.bmp.cer))) {
        success <- FALSE
        CV.maxsteps <- NA
        CV.support <- rep(x=NA, times=CV.Lm)
        CV.size <- rep(x=NA, times=CV.Lm)
        CV.bmp.lhr <- rep(x=NA, times=CV.Lm)
        CV.bmp.lrt <- rep(x=NA, times=CV.Lm)
        CV.bmp.cer <- rep(x=NA, times=CV.Lm)
        CV.grp.lhr <- rep(x=NA, times=CV.Lm)
        CV.grp.lrt <- rep(x=NA, times=CV.Lm)
        CV.grp.cer <- rep(x=NA, times=CV.Lm)
        CV.time.bar <- rep(x=NA, times=CV.Lm)
        CV.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.max.time.bar <- rep(x=NA, times=CV.Lm)
        CV.min.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p)
        CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
        CV.trace <- rep(x=NA, times=CV.Lm)
        CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p))
      } else {
        # Cross-validated optimal length from all folds by minimization of the CER (between predicted and observed inbox test samples survival times)
        success <- TRUE
      }
    } else if (peelcriterion == "grp") {
      if (all(is.na(CV.grp.cer)) || all(is.nan(CV.grp.cer))) {
        success <- FALSE
        CV.maxsteps <- NA
        CV.support <- rep(x=NA, times=CV.Lm)
        CV.size <- rep(x=NA, times=CV.Lm)
        CV.bmp.lhr <- rep(x=NA, times=CV.Lm)
        CV.bmp.lrt <- rep(x=NA, times=CV.Lm)
        CV.bmp.cer <- rep(x=NA, times=CV.Lm)
        CV.grp.lhr <- rep(x=NA, times=CV.Lm)
        CV.grp.lrt <- rep(x=NA, times=CV.Lm)
        CV.grp.cer <- rep(x=NA, times=CV.Lm)
        CV.time.bar <- rep(x=NA, times=CV.Lm)
        CV.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.max.time.bar <- rep(x=NA, times=CV.Lm)
        CV.min.prob.bar <- rep(x=NA, times=CV.Lm)
        CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p)
        CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
        CV.trace <- rep(x=NA, times=CV.Lm)
        CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p))
      } else {
        # Cross-validated optimal length from all folds by minimization of the CER (between predicted and observed inbox test samples survival times)
        success <- TRUE
      }
    } else {
      stop("Invalid peeling criterion. Exiting ...\n\n")
    }
  } else {
    stop("Invalid CV criterion option. Exiting ... \n\n")
  }
  
  # Box statistics for each step
  CV.stats <-  data.frame("cv.support"=CV.support,
                          "cv.size"=CV.size,
                          "cv.bmp.lhr"=CV.bmp.lhr,
                          "cv.bmp.lrt"=CV.bmp.lrt,
                          "cv.bmp.cer"=CV.bmp.cer,
                          "cv.grp.lhr"=CV.grp.lhr,
                          "cv.grp.lrt"=CV.grp.lrt,
                          "cv.grp.cer"=CV.grp.cer,
                          "cv.time.bar"=CV.time.bar,
                          "cv.prob.bar"=CV.prob.bar,
                          "cv.max.time.bar"=CV.max.time.bar,
                          "cv.min.prob.bar"=CV.min.prob.bar)
  rownames(CV.stats) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.maxsteps <- CV.Lm
  
  # Create the return object 'CV.fit'
  CV.fit <- list("cv.maxsteps"=CV.maxsteps,
                 "cv.boxcut"=CV.boxcut,
                 "cv.boxind"=CV.boxind,
                 "cv.trace"=CV.trace,
                 "cv.rules"=CV.rules,
                 "cv.stats"=CV.stats)
  
  return(list("X"=X,
              "y"=y,
              "delta"=delta,
              "cvfit"=CV.fit,
              "success"=success,
              "seed"=seed))
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    cv.ave.peel (traindata,
#                                 trainstatus,
#                                 traintime,
#                                 traingroups, 
#                                 testdata,
#                                 teststatus,
#                                 testtime,
#                                 testgroups,
#                                 varsign,
#                                 initcutpts,
#                                 arg,
#                                 probval,
#                                 timeval)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

cv.ave.peel <- function(traindata,
                        trainstatus,
                        traintime,
                        traingroups, 
                        testdata,
                        teststatus,
                        testtime,
                        testgroups,
                        varsign,
                        initcutpts,
                        arg,
                        probval,
                        timeval) {
  
  # Parsing and evaluating 'arg' parameters to evaluate 'peelcriterion'
  alpha <- NULL
  beta <- NULL
  L <- NULL
  peelcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=arg, split=",")) ))
  
  # Training the model
  peelobj <- prsp(traindata=traindata, 
                  traintime=traintime, 
                  trainstatus=trainstatus,
                  varsign=varsign, 
                  initcutpts=initcutpts, 
                  arg=arg, 
                  traingroups=traingroups)
  maxsteps <- peelobj$maxsteps
  
  # Compute the box statistics for all steps, each entry or row signifies a step
  boxstat <- vector(mode="list", length=maxsteps)
  timemat <- matrix(NA, nrow=maxsteps, ncol=nrow(testdata))
  probmat <- matrix(NA, nrow=maxsteps, ncol=nrow(testdata))
  ind.rem <- numeric(0)
  grp1 <- 1*(testgroups == levels(testgroups)[1])
  wg <- which(grp1 == 1)
  for (l in 1:maxsteps) {
    # Extract the rule and sign as one vector
    boxcut <- peelobj$boxcut[l, ] * varsign
    test.cut <- t(t(testdata) * varsign)
    test.ind <- t(t(test.cut) >= boxcut)
    # Set as TRUE which test observations are TRUE for all covariates
    box1 <- 1*(rowMeans(test.ind) == 1)
    wb <- which(box1 == 1)
    if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
      if ((sum(box1, na.rm=TRUE) == length(box1[!is.na(box1)])) || (sum(box1, na.rm=TRUE) == 0)) {
        surv.formula <- (survival::Surv(testtime, teststatus) ~ 1 + box1)
        surv.fit <- survival::survfit(surv.formula, na.action="na.omit")
        timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
        probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
        bmplhr <- 0
        bmplrt <- 0
        bmpcer <- 1
        grplhr <- NA 
        grplrt <- NA
        grpcer <- NA
      } else if ((sum(box1, na.rm=TRUE) != length(box1[!is.na(box1)])) && (sum(box1, na.rm=TRUE) != 0)) {
        surv.formula <- (survival::Surv(testtime, teststatus) ~ 1 + box1)
        surv.fit <- survival::survfit(surv.formula, na.action="na.omit")
        timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
        probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
        coxobj <- survival::coxph(surv.formula, na.action="na.omit", iter.max=20, timefix=TRUE, method="efron")
        bmplhr <- coxobj$coef
        bmplrt <- survival::survdiff(surv.formula, na.action="na.omit", rho=0)$chisq
        predobj <- predict(object=coxobj, type="lp", reference="sample", na.action="na.omit")
        bmpcer <- Hmisc::rcorr.cens(x=predobj, S=Surv(testtime, teststatus))['C Index']
        grplhr <- NA 
        grplrt <- NA
        grpcer <- NA
      } else {
        timemat[l, ] <- NA
        probmat[l, ] <- NA
        bmplhr <- NA
        bmplrt <- NA
        bmpcer <- NA
        grplhr <- NA 
        grplrt <- NA
        grpcer <- NA
        ind.rem <- c(ind.rem, l)
      }
    } else if (peelcriterion == "grp") {
      if ((sum((grp1 * box1), na.rm=TRUE) == length((grp1 * box1)[!is.na((grp1 * box1))])) || 
          (sum((grp1 * box1), na.rm=TRUE) == 0)) {
        timemat[l, ] <- NA
        probmat[l, ] <- NA
        bmplhr <- NA
        bmplrt <- NA
        bmpcer <- NA
        grplhr <- 0
        grplrt <- NA
        grpcer <- 1
      } else  if ((sum((grp1 * box1), na.rm=TRUE) != length((grp1 * box1)[!is.na((grp1 * box1))])) && 
                  (sum((grp1 * box1), na.rm=TRUE) != 0)) {
        coxobj <- survival::coxph(survival::Surv(testtime, teststatus) ~ 1 + grp1 * box1, 
                                  data=list("grp1"=grp1, "box1"=box1),
                                  na.action="na.omit", 
                                  iter.max=20,
                                  timefix=TRUE, 
                                  method="efron")
        predobj <- predict(object=coxobj, type="lp", reference="sample", na.action="na.omit")
        timemat[l, ] <- NA
        probmat[l, ] <- NA
        bmplhr <- NA
        bmplrt <- NA
        bmpcer <- NA
        grplhr <- coxobj$coef["grp1:box1"]
        grplrt <- NA
        grpcer <- Hmisc::rcorr.cens(x=predobj, S=survival::Surv(testtime, teststatus))['C Index']
      } else {
        timemat[l, ] <- NA
        probmat[l, ] <- NA
        bmplhr <- NA
        bmplrt <- NA
        bmpcer <- NA
        grplhr <- NA 
        grplrt <- NA
        grpcer <- NA
        ind.rem <- c(ind.rem, l)
      }
    } else {
      stop("Invalid peeling criterion. Exiting ...\n\n")
    }
    boxstat[[l]] <- c(bmplhr, bmplrt, bmpcer, grplhr, grplrt, grpcer)
    names(boxstat[[l]]) <- NULL
  }
  if ((length(ind.rem) == maxsteps) || all(is.na(timemat)) || all(is.na(probmat))) {
    for (l in 1:maxsteps) {
      boxstat[[l]] <- rep(x=NA, times=10)
    }
  } else {
    endobj <- endpoints (ind=ind.rem, timemat=timemat, probmat=probmat, timeval=timeval, probval=probval)
    time.bar <- endobj$time.bar
    prob.bar <- endobj$prob.bar
    max.time.bar <- endobj$max.time.bar
    min.prob.bar <- endobj$min.prob.bar
    for (l in 1:maxsteps) {
      if (!(l %in% ind.rem)) {
        boxstat[[l]] <- c(boxstat[[l]], time.bar[l], prob.bar[l], max.time.bar[l], min.prob.bar[l])
      } else {
        boxstat[[l]] <- c(boxstat[[l]], NA, NA, NA, NA)
      }
    }
  }
  
  return(list("maxsteps"=maxsteps, "boxstat"=boxstat, "boxcut"=peelobj$boxcut, "trace"=peelobj$vartrace))
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    cv.comb.peel (traindata,
#                                  trainstatus,
#                                  traintime,
#                                  traingroups,
#                                  testdata,
#                                  teststatus,
#                                  testtime,
#                                  testgroups,
#                                  varsign,
#                                  initcutpts,
#                                  arg)
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

cv.comb.peel <- function(traindata,
                         trainstatus,
                         traintime,
                         traingroups,
                         testdata,
                         teststatus,
                         testgroups,
                         testtime,
                         varsign,
                         initcutpts,
                         arg) {
  
  # Training the model
  peelobj <- prsp(traindata=traindata, 
                  traintime=traintime, 
                  trainstatus=trainstatus,
                  varsign=varsign, 
                  initcutpts=initcutpts, 
                  arg=arg, 
                  traingroups=traingroups)
  maxsteps <- peelobj$maxsteps
  
  # Create the indicator matrix of the test data that is within the box for each step
  boxind <- matrix(NA, nrow=maxsteps, ncol=nrow(testdata))
  for (l in 1:maxsteps) {
    # Extract the rule and sign as one vector
    boxcut <- peelobj$boxcut[l, ] * varsign
    test.cut <- t(t(testdata) * varsign)
    test.ind <- t(t(test.cut) >= boxcut)
    # Set as TRUE which test observations are TRUE for all covariates
    boxind[l, ] <- (rowMeans(test.ind) == 1)
  }
  
  return(list("boxind"=boxind, "maxsteps"=maxsteps, "boxcut"=peelobj$boxcut, "trace"=peelobj$vartrace))
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    prsp (traindata,
#                          trainstatus,
#                          traintime,
#                          traingroups,
#                          varsign,
#                          initcutpts,
#                          arg)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

prsp <- function(traindata,
                 traintime,
                 trainstatus,
                 traingroups,
                 varsign,
                 initcutpts,
                 arg) {
  
  # Parsing and evaluating PRSP parameters
  alpha <- NULL
  beta <- NULL
  L <- NULL
  peelcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=arg, split=",")) ))
  
  # Ensures that the training data is a numeric matrix
  traindata <- as.matrix(traindata)
  mode(traindata) <- "numeric"
  
  # Constants
  n <- nrow(traindata)                                   # Number of samples
  p <- ncol(traindata)                                   # Number of initially screened covariates
  Lmax <- ceiling(log(beta) / log(1 - alpha))            # Maximal possible number of peeling steps
  
  # Directions of peeling
  if (is.null(varsign)) {
    varsign <- apply(X=traindata,
                     MARGIN=2,
                     FUN=function(x) {sign(survival::coxph(survival::Surv(traintime, trainstatus) ~ x, 
                                                           eps=0.01, 
                                                           iter.max=20,
                                                           na.action="na.omit", 
                                                           singular.ok=TRUE, 
                                                           timefix=TRUE, 
                                                           method="efron")$coef)})
  }
  names(varsign) <- colnames(traindata)
  
  # Initialization of box boundaries
  if (is.null(initcutpts)) {
    initcutpts <- numeric(p)
    for (j in 1:p) {
      if (varsign[j] == 1) {
        initcutpts[j] <- min(traindata[,j], na.rm=TRUE)
      } else if (varsign[j] == -1) {
        initcutpts[j] <- max(traindata[,j], na.rm=TRUE)
      } else {
        warning("Covariate `", colnames(traindata)[j], "` has no direction of peeling! \n\n", sep="")
      }
    }
  }
  
  # Initialization of returned objects
  vartrace <- numeric(Lmax)
  boxcut <- matrix(data=NA, nrow=Lmax, ncol=p)
  
  # Initialization for the PRSP loop
  lhrlj <- matrix(NA, Lmax, p)
  lrtlj <- matrix(NA, Lmax, p)
  chslj <- matrix(NA, Lmax, p)
  boxcutpts <- initcutpts
  boxmass <- 1
  sel <- rep(x=TRUE, times=n)                            # Initial logical vector of in-box samples
  xsel <- traindata                                      # Initial selection of samples from training data
  l <- 1
  grp1 <- 1 * (traingroups == levels(traingroups)[1])    # Selection of e.g. G1
  wg <- which(grp1 == 1)
  
  # PRSP or PGRSP loop
  while ((boxmass >= beta) & (l <= L)) {
    xsign <- t(t(xsel) * varsign)
    # All potential quantile peelings by dimension
    cutpts.sign <- updatecut(X=xsign, fract=alpha)
    cutpts <- cutpts.sign * varsign
    # Update box membership indicator by dimension
    boxes <- as.matrix(t((t(traindata) * varsign) >= as.vector(cutpts.sign)) & sel)
    vmd <- rep(x=NA, times=p)
    for (j in 1:p) {
      box1j <- 1 * boxes[,j]
      wb <- which(box1j == 1)
      if (peelcriterion == "lhr") {
        # Rate of increase of LHR between in and out bump
        if ((sum(box1j, na.rm=TRUE) != length(box1j[!is.na(box1j)])) &&
            (sum(box1j, na.rm=TRUE) != 0)) {
          fit <- survival::coxph(survival::Surv(traintime, trainstatus) ~ 1 + box1j, 
                                 na.action="na.omit",
                                 iter.max=20,
                                 timefix=TRUE, 
                                 method="efron")
          lhrlj[l,j] <- fit$coef
          if (l == 1) {
            vmd[j] <- (lhrlj[l,j] - 0) / (1 - mean(box1j))
          } else {
            vmd[j] <- (lhrlj[l,j] - lhrlj[l-1,j]) / (boxmass - mean(box1j))
          }
        }
      } else if (peelcriterion == "lrt") {
        # Rate of increase of LRT between in and out bump
        if ((sum(box1j, na.rm=TRUE) != length(box1j[!is.na(box1j)])) &&
            (sum(box1j, na.rm=TRUE) != 0)) {
          fit <- survival::survdiff(survival::Surv(traintime, trainstatus) ~ 1 + box1j, 
                                    na.action="na.omit", 
                                    rho=0)
          lrtlj[l,j] <- fit$chisq
          if (l == 1) {
            vmd[j] <- (lrtlj[l,j] - 0) / (1 - mean(box1j))
          } else {
            vmd[j] <- (lrtlj[l,j] - lrtlj[l-1,j]) / (boxmass - mean(box1j))
          }
        }
      } else if (peelcriterion == "chs") {
        # Rate of increase of CHS between in and out bump
        if ((sum(box1j, na.rm=TRUE) != length(box1j[!is.na(box1j)])) &&
            (sum(box1j, na.rm=TRUE) != 0)) {
          fit <- survival::survfit(formula=survival::Surv(traintime, trainstatus) ~ 1 + box1j, 
                                   na.action="na.omit")
          chslj[l,j] <- sum(cumsum(fit$n.event/fit$n.risk))
          if (l == 1) {
            vmd[j] <- (chslj[l,j] - 0) / (1 - mean(box1j))
          } else {
            vmd[j] <- (chslj[l,j] - chslj[l-1,j]) / (boxmass - mean(box1j))
          }
        }
      } else if (peelcriterion == "grp") {
        # Rate of increase of LHR of interaction effect between the groups and the bumps
        if ((sum((grp1 * box1j), na.rm=TRUE) != length((grp1 * box1j)[!is.na(grp1 * box1j)])) && 
            (sum((grp1 * box1j), na.rm=TRUE) != 0)) {
          fit <- survival::coxph(survival::Surv(traintime, trainstatus) ~ 1 + grp1 * box1j, 
                                 data=list("grp1"=grp1, "box1j"=box1j),
                                 iter.max=20,
                                 na.action="na.omit",
                                 timefix=TRUE, 
                                 method="efron")
          lhrlj[l,j] <- fit$coef["grp1:box1j"]
          if (l == 1) {
            vmd[j] <- (lhrlj[l,j] - 0) / (1 - mean(grp1 * box1j))
          } else {
            vmd[j] <- (lhrlj[l,j] - lhrlj[l-1,j]) / (boxmass - mean(grp1 * box1j))
          }
        }
      } else {
        stop("Invalid peeling criterion. Exiting ...\n\n")
      }
    }
    # If the previous attempted peeling succeeded
    # Maximize the rate of increase of peeling criterion
    # Only one variable has to be selected in case of ties
    if ((!is.empty(vmd[(!is.nan(vmd)) & (!is.infinite(vmd)) & (!is.na(vmd))]))) {
      varj <- which(vmd == max(vmd[(!is.nan(vmd)) & (!is.infinite(vmd)) & (!is.na(vmd))]))
      if (length(varj) > 1) {
        varj <- varj[1] #or varj <- sample(x=varj, size=1)
      }
    } else {
      varj <- sample(x=1:p, size=1)
    }
    # Updating
    sel <- boxes[, varj, drop=TRUE]
    boxmass <- mean(1 * sel, na.rm=TRUE)
    xsel <- traindata[sel, ,drop=FALSE]
    boxcutpts[varj] <- cutpts[varj]
    # Saving trained box quantities of interest for the current peeling step
    boxcut[l, ] <- boxcutpts
    vartrace[l] <- varj
    l <- l + 1
  }
  l <- l - 1
  
  # Prepending the first-step box covering all the training data
  boxcut <- rbind(initcutpts, boxcut[1:l, , drop=FALSE])
  vartrace <- c(0, vartrace[1:l])
  rownames(boxcut) <- paste("step", 0:l, sep="")
  colnames(boxcut) <- colnames(traindata)
  names(vartrace) <- paste("step", 0:l, sep="")
  
  return(list("maxsteps"=l+1,
              "boxcut"=boxcut,
              "vartrace"=vartrace,
              "varsign"=varsign))
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    endpoints (ind, timemat, probmat, timeval, probval)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

endpoints <- function(ind, timemat, probmat, timeval, probval) {
  
  if (!(is.empty(ind))) {
    timemat <- timemat[-ind, , drop=FALSE]
    probmat <- probmat[-ind, , drop=FALSE]
  }
  L <- nrow(timemat) #or L <- nrow(probmat)
  min.prob.bar <- apply(probmat, 1, min, na.rm=TRUE)
  max.time.bar <- apply(timemat, 1, max, na.rm=TRUE)
  if (is.null(probval) && is.null(timeval)) {
    prob.bar <- rep(x=NA, times=L)
    time.bar <- rep(x=NA, times=L)
  } else if (!is.null(probval)) {
    prob.bar <- rep(x=probval, times=L)
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
    time.bar <- rep(x=timeval, times=L)
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
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    updatecut (X, fract)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

updatecut <- function(X, fract) {
  
  p <- dim(X)[2]
  cutpts <- apply(X, 2, "quantile", type=7, probs=fract, na.rm=TRUE)
  
  for (j in 1:p) {
    xunique <- sort(unique(X[, j]))
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
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    cv.folds (y, n, K=5, seed=NULL)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

cv.folds <- function (y, n, K=5, seed=NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  K <- round(rep(x=K, length.out=1))
  
  subsample <- function(index, K) {
    permindex <- sample(x=length(index), replace=FALSE, prob=NULL)
    w <- rep(x=seq_len(K), length.out=length(index))
    out <- list("index"=index[permindex], "which"=w)
    return(out)
  }
  
  regularcv <- function(n, K) {
    if (!isTRUE((K >= 1) && (K <= n)))
      stop(paste("`K` outside allowable range {1,...,", n, "}. Exiting ... \n\n", sep=""))
    if (K == 1) {
      #cat(paste("No cross-validation requested\n", sep=""))
      perm <- seq_len(n)
      fold <- list("index"=perm, "which"=rep(x=1, length.out=n))
    } else if (K == n) {
      #cat(paste("Leave-one-out cross-validation requested\n", sep=""))
      perm <- seq_len(n)
      fold <- list("index"=perm, "which"=seq_len(n))
    } else {
      #cat(paste("Regular ", K, "-fold cross-validation will be carried out\n", sep=""))
      perm <- sample(x=n, size=length(1:n), replace=FALSE, prob=NULL)
      fold <- list("index"=perm, "which"=rep(x=seq_len(K), length.out=n))
    }
    return(fold)
  }
  
  if ((missing(y)) && (missing(n))) {
    
    stop("Both the number of observations `n` and the outcome `y` are missing. Exiting ... \n\n")
    
  } else if ((missing(y)) && (!missing(n))) {
    
    fold <- regularcv(n=round(rep(x=n, length.out=1)), K=K)
    
  } else if ((missing(n)) && (!missing(y))) {
    
    y <- as.numeric(y)
    n <- length(y)
    yf <- factor(x=y, levels=unique(as.character(y)))
    ylev <- levels(yf)
    ynlev <- nlevels(yf)
    
    if ((ynlev == 0) || (ynlev == n)) {
      
      fold <- regularcv(n=n, K=K)
      
    } else {
      
      ytab <- table(y)
      if (length(ytab) == 1)
        warning("One class of the outcome has no records and will be ignored.\n")
      if (K == 1) {
        #cat(paste("No cross-validation requested\n", sep=""))
        fold <- list("index"=seq_len(n), "which"=rep(x=1, length.out=n))
      } else if (K == n) {
        #cat(paste("Leave-one-out cross-validation requested\n", sep=""))
        fold <- list("index"=seq_len(n), "which"=seq_len(n))
      } else {
        #cat(paste("Stratified ", K, "-fold cross-validation will be carried out\n", sep=""))
        index <- seq(along = y)
        indexlist <- vector(mode="list", length=ynlev)
        for (l in 1:ynlev) {
          indexlist[[l]] <- index[y == ylev[l]]
        }
        foldlist <- lapply(X=indexlist, FUN=subsample, K=K)
        fold <- list("index"=numeric(0), "which"=numeric(0))
        for (l in 1:ynlev) {
          fold$index <- c(fold$index, foldlist[[l]]$index)
          fold$which <- c(fold$which, foldlist[[l]]$which)
        }
      }
      
    }
    
  }
  
  permkey <- pmatch(x=1:n, table=fold$index)
  ord <- numeric(0)
  for (k in 1:K) {
    ord <- c(ord, fold$index[(fold$which == k)])
  }
  foldkey <- pmatch(x=1:n, table=ord)
  folds <- list(n=n, K=K, perm=fold$index, permkey=permkey, which=fold$which, foldkey=foldkey, seed=seed)
  
  return(folds)
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    lapply.array (X, rowtrunc=NULL, coltrunc=NULL,
#                                  sub=NULL, fill=NA, MARGIN=1:2, FUN, ...)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

lapply.array <- function (X, rowtrunc=NULL, coltrunc=NULL,
                          sub=NULL, fill=NA, MARGIN=1:2, FUN, ...) {
  
  A <- list2array(list=X, rowtrunc=rowtrunc, coltrunc=coltrunc, sub=sub, fill=fill)
  return(apply(X=A, MARGIN=MARGIN, FUN=FUN, ...))
  
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    lapply.mat (X, coltrunc=NULL,
#                                sub=NULL, fill=NA, MARGIN=2, FUN, ...)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

lapply.mat <- function (X, coltrunc=NULL,
                        sub=NULL, fill=NA, MARGIN=2, FUN, ...) {
  
  M <- list2mat(list=X, coltrunc=coltrunc, sub=sub, fill=fill)
  return(apply(X=M, MARGIN=MARGIN, FUN=FUN, ...))
  
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    list2array (list, rowtrunc=NULL, coltrunc=NULL, sub=NULL, fill=NA)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

list2array <- function (list, rowtrunc=NULL, coltrunc=NULL, sub=NULL, fill=NA) {
  
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
    min.row <- min(sapply(my.list, nrow), na.rm=TRUE)
    max.row <- max(sapply(my.list, nrow), na.rm=TRUE)
    min.col <- min(sapply(my.list, ncol), na.rm=TRUE)
    max.col <- max(sapply(my.list, ncol), na.rm=TRUE)
    if (!is.null(coltrunc)) {
      if (coltrunc == "min") {
        adjusted.list <- lapply(my.list, function(x) {x[,1:min.col,drop=FALSE]})
      } else if (coltrunc == "max") {
        adjusted.list <- lapply(my.list, function(x) {cbind(x, matrix(data=fill, nrow=nrow(x), ncol=max.col - ncol(x)))})
      } else {
        adjusted.list <- lapply(my.list, function(x, coltrunc) {if (coltrunc <= ncol(x)) {
          x[,1:coltrunc,drop=FALSE]
        } else if (coltrunc > ncol(x)) {
          cbind(x, matrix(data=fill, nrow=nrow(x), ncol=coltrunc - ncol(x)))
        }
        }, coltrunc)
      }
    } else {
      adjusted.list <- lapply(my.list, function(x) {cbind(x, matrix(data=fill, nrow=nrow(x), ncol=max.col - ncol(x)))})
    }
    if (!is.null(rowtrunc)) {
      if (rowtrunc == "min") {
        adjusted.list <- lapply(adjusted.list, function(x) {x[1:min.row,,drop=FALSE]})
      } else if (rowtrunc == "max") {
        adjusted.list <- lapply(adjusted.list, function(x) {rbind(x, matrix(data=fill, nrow=max.row - nrow(x), ncol=ncol(x)))})
      } else {
        
        adjusted.list <- lapply(my.list, function(x, rowtrunc) {if (rowtrunc <= nrow(x)) {
          x[1:rowtrunc,,drop=FALSE]
        } else if (rowtrunc > nrow(x)) {
          rbind(x, matrix(data=fill, nrow=rowtrunc - nrow(x), ncol=ncol(x)))
        }
        }, rowtrunc)
      }
    } else {
      adjusted.list <- lapply(adjusted.list, function(x) {rbind(x, matrix(data=fill, nrow=max.row - nrow(x), ncol=ncol(x)))})
    }
    my.array <- array(data=fill, dim=c(nrow(adjusted.list[[1]]), ncol(adjusted.list[[1]]), length(adjusted.list)))
    for(i in 1:length(adjusted.list)) {
      my.array[,,i] <- adjusted.list[[i]]
    }
  } else {
    my.array <- array(data=fill, dim=c(0,0,0))
  }
  
  return(my.array)
  
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    list2mat (list, coltrunc=NULL, sub=NULL, fill=NA)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

list2mat <- function (list, coltrunc=NULL, sub=NULL, fill=NA) {
  
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
    min.col <- min(sapply(my.list, length), na.rm=TRUE)
    max.col <- max(sapply(my.list, length), na.rm=TRUE)
    if (!is.null(coltrunc)) {
      if (coltrunc == "min") {
        adjusted.list <- lapply(my.list, function(x) {x[1:min.col]})
      } else if (coltrunc == "max") {
        adjusted.list <- lapply(my.list, function(x) {c(x, rep(x=fill, times=max.col-length(x)))})
      } else {
        adjusted.list <- lapply(my.list, function(x, coltrunc) {
          if (coltrunc <= length(x)) {
            x[1:coltrunc]
          } else if (coltrunc > length(x)) {
            c(x, rep(x=fill, times=coltrunc-length(x)))
          }
        }, coltrunc)
      }
    } else {
      adjusted.list <- lapply(my.list, function(x) {c(x, rep(x=fill, times=max.col-length(x)))})
    }
    my.mat <- do.call(rbind, adjusted.list)
  } else {
    my.mat <- matrix(data=fill, nrow=0, ncol=0)
  }
  
  return(my.mat)
  
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    cbindlist (list, trunc)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

cbindlist <- function(list, trunc) {
  
  if (!is.empty(list)) {
  
    w <- unlist(lapply(X=list, FUN=function(x) {ifelse(test=is.empty(x), yes=TRUE, no=FALSE)}))
    list <- list[!w]
    max.row <- max(unlist(sapply(X=list, FUN=nrow)), na.rm=TRUE)
    adjusted.list <- lapply(X=list, FUN=function(x) {rbind(x, matrix(data=NA, nrow=max.row - nrow(x), ncol=ncol(x)))})
    my.mat <- adjusted.list[[1]]
    lcl <- length(adjusted.list)
    if (lcl > 1) {
      for(i in 2:lcl){
        my.mat <- cbind(my.mat, adjusted.list[[i]])
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
#===============================================================================================================================#



#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    is.empty(x)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

is.empty <- function(x) {
  
  if (is.null(x)) {
    return(TRUE)
  } else if (is.vector(x, mode="any")) {
    na <- is.na(x)
    if (length(which(na)) != 0) {
      return(FALSE)  #NA is a non-empty value
    } else {
      x <- x[!na]
      if (is.character(x)) {
        if (length(x) == 0) {
          return(TRUE)
        } else if (length(x) > 0) {
          return( all(sapply(as.list(x), function(x) {x == ""})) )
        }
      } else {
        if (length(x) == 0) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      }      
    }
  } else if (is.list(x)) {
    return( is.null(x) )   #NA is a non-empty value
  } else if (is.matrix(x) || is.data.frame(x)) {
    return( ((nrow(x) == 0) || (ncol(x) == 0)) )
  } 
  
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    is.wholenumber(x, tol=.Machine$double.eps^0.5)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

is.wholenumber <- function(x, tol=.Machine$double.eps^0.5) {
  (abs(x - round(x)) < tol)
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    zeroslope(y, x, lag=1, span=0.10, degree=2, family="gaussian", minimum=TRUE)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

zeroslope <- function(y, x, lag, span, degree, family, minimum) {
  n <- length(x)
  if (floor(n * span) == 0) {
    stop("Span 'span' is too small with respect to number of peeling steps. Exiting ...\n\n")
  }
  if (anyNA(x)) {
    stop("Peeling steps cannot contain NA values. Exiting ... \n\n")
  } else {
    y <- y[order(x)]  # reorder the data in ascending values of x
    x <- x[order(x)]  # do the same for x
    na <- (is.na(y) | is.nan(y))
    wa <- which(na)
    if (!is.empty(wa)) {
      xc <- x[-wa]
      yc <- y[-wa]
    } else {
      xc <- x
      yc <- y
    }
    fitc <- loess(yc ~ as.numeric(xc), na.action="na.omit", span=span, degree=degree, family=family)$fitted
    loe <- rep(NA, length(x))
    loe[!na] <- fitc
    d <- diff(x=loe, lag=lag)/diff(x=x, lag=lag)
    d.sign <- diff(x=sign(d), lag=lag)
    zs <- which(d.sign != 0) + lag
    if (minimum) {
      w <- which.min(loe[zs])
    } else {
      w <- which.max(loe[zs])
    }
    return(zs[w])
  }
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    .onAttach (libname, pkgname)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

.onAttach <- function(libname, pkgname) {
  
  SSver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
  packageStartupMessage(paste(pkgname, " ", SSver, " and > 0.7.0 are major releases with significant user-visible changes.", sep=""))
  packageStartupMessage("Type PRIMsrc.news() to see new features, changes, and bug fixes.")
  
}
#===============================================================================================================================#
