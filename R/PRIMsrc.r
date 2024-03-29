#===============================================================================================================================#
# PRIMsrc
#===============================================================================================================================#

#===============================================================================================================================#
# 1. END-USER SURVIVAL BUMP HUNTING FUNCTION
#===============================================================================================================================#

#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                   sbh(X,
#                       y,
#                       groups=NULL,
#                       delta,
#                       B=30,
#                       K=5,
#                       A=1000,
#                       vs=TRUE,
#                       vstype="ppl",
#                       vsarg="alpha=1,
#                              nalpha=1,
#                              nlambda=100",
#                       cv=TRUE,
#                       cvtype="combined",
#                       cvarg="alpha=0.01,
#                              beta=0.10,
#                              peelcriterion=\"lrt\",
#                              cvcriterion=\"cer\"",
#                       pv=FALSE,
#                       control=sbh.control(vscons = 0.5, 
#                                           decimals = 2, 
#                                           onese = FALSE, 
#                                           probval = NULL, 
#                                           timeval = NULL, 
#                                           lag = 2, 
#                                           span = 0.10, 
#                                           degree = 2),
#                       parallel.vs=FALSE,
#                       parallel.rep=FALSE,
#                       parallel.pv=FALSE,
#                       conf=NULL,
#                       verbose=TRUE,
#                       seed=NULL)
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

sbh <- function(X,
                y,
                groups=NULL,
                delta,
                B=30,
                K=5,
                A=1000,
                vs=TRUE,
                vstype="ppl",
                vsarg="alpha=1,
               nalpha=1,
               nlambda=100",
                cv=TRUE,
                cvtype="combined",
                cvarg="alpha=0.01,
               beta=0.10,
               peelcriterion=\"lrt\",
               cvcriterion=\"cer\"",
                pv=FALSE,
                control=sbh.control(vscons = 0.5, 
                                    decimals = 2, 
                                    onese = FALSE, 
                                    probval = NULL, 
                                    timeval = NULL, 
                                    lag = 2, 
                                    span = 0.10, 
                                    degree = 2),
                parallel.vs=FALSE,
                parallel.rep=FALSE,
                parallel.pv=FALSE,
                conf=NULL,
                verbose=TRUE,
                seed=NULL) {
  
  # Parsing and evaluating SBH parameters
  alpha <- NULL
  beta <- NULL
  peelcriterion <- NULL
  cvcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=cvarg, split=",")) ))
  
  # Evaluating ancillary SBH parameters
  vscons <- control$vscons
  decimals <- control$decimals
  onese <- control$onese
  probval <- control$probval
  timeval <- control$timeval
  lag <- control$lag
  span <- control$span
  degree <- control$degree
  
  # Checking Inputs of Parameters
  if ((!is.wholenumber(B)) || (B <= 0)) {
    stop("\nArgument 'B' must be a positive integer. Exiting ... \n\n")
  }
  if ((!is.wholenumber(K)) || (K <= 0)) {
    stop("\nArgument 'K' must be a positive integer. Exiting ... \n\n")
  }
  if ((!is.wholenumber(A)) || (A <= 0)) {
    stop("\nArgument 'A' must be a positive integer. Exiting ... \n\n")
  }
  peelcriterion <- match.arg(arg=peelcriterion, choices=c("lrt", "lhr", "cer", "grp"), several.ok=FALSE)
  if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
    cvcriterion <- match.arg(arg=cvcriterion, choices=c("lrt", "lhr", "cer"), several.ok=FALSE)
    groups <- NULL
  } else if (peelcriterion == "grp") {
    cvcriterion <- match.arg(arg=cvcriterion, choices=c("cer"), several.ok=FALSE)
    if (is.null(groups)) {
      stop("\nArgument `groups` must be specified when used with PRGSP algorithm. Exiting ... \n\n")
    }
    groups <- as.factor(groups)
    if (!is.null(groups) && (length(levels(groups)) != 2)) {
      stop("\nArgument `groups` must have exactly two levels when used with PRGSP algorithm. Exiting ... \n\n")
    }
  } else {
    stop("Invalid peeling criterion option. Exiting ... \n\n")
  }
  if (vs) {
    vstype <- match.arg(arg=vstype, choices=c("ppl", "pcqr", "spca", "prsp"), several.ok=FALSE)
    if ((vscons < 1/K) || (vscons > 1)) {
      stop("\nVariable screening conservativeness must be a real number in [1/K, 1]. Exiting ... \n\n")
    }
  } else {
    vstype <- NA
  }
  if (cv) {
    cvtype <- match.arg(arg=cvtype, choices=c("combined", "averaged"), several.ok=FALSE)
  } else {
    cvtype <- "combined"
    K <- 1
    vscons <- 1
  }
  if (!is.null(probval)) {
    stopifnot (length(probval) == 1L, !is.na(probval), probval >= 0, probval <= 1) 
  } 
  if (!is.null(timeval)) {
    stopifnot (length(timeval) == 1L, !is.na(timeval), timeval >= 0) 
  }
  
  # Checking Inputs of Data
  if (missing(X)) {
    stop("\nNo explanatory variables provided! Exiting ... \n\n")
  }
  if (missing(y)) {
    stop("\nNo outcome variable provided! Exiting ... \n\n")
  }
  if (missing(delta)) {
    stop("\nNo censoring variable provided! Exiting ... \n\n")
  }
  
  # Constants & Metaparameters
  n <- nrow(X)
  p <- ncol(X)
  if (alpha < 1/n) {
    cat("Warning: Parameter `alpha` is less than the minimal allowed value. It is reset to ", 1/n, ".\n\n", sep="")
  }
  alpha <- max(1/n, alpha)                               # Minimal peeling quantile
  if (beta < 1/n) {
    cat("Warning: Parameter `beta` is less than the minimal allowed value. It is reset to ", 1/n, ".\n\n", sep="")
    beta <- max(1/n, beta)                               # Minimal box support
  }
  Lmax <- ceiling(log(beta) / log(1 - alpha))            # Maximal possible number of peeling steps
  
  # Preparing the data
  if (!is.matrix(X)) {
    if (is.data.frame(X)) {
      if (!all(sapply(X, is.numeric))) {
        stop("Input data.frame 'X' is not all numeric.\n\n")
      } else { 
        cat("Warning: Input data.frame `X` was corced to numeric matrix.\n\n", sep="")
        X <- as.matrix(X)
        mode(X) <- "numeric"
      }
    }
  }
  
  digits <- getOption("digits")
  y[y <= 0] <- 10^(-digits)
  
  if (is.null(rownames(X))) {
    rownames(X) <- 1:n
  }
  if (is.null(colnames(X))) {
    colnames(X) <- paste("X", 1:p, sep="")
  }
  
  # Summarizing user choices
  if (!cv) {
    if (B > 1) {
      if (parallel.rep) {
        cat("\nRequested parallel replicated procedure with ", B, " replications without cross-validation. \n", sep="")
      } else {
        cat("\nRequested serial replicated procedure with ", B, " replications without cross-validation. \n", sep="")
      }
    } else {
      cat("\nRequested single procedure without cross-validation and replication. \n", sep="")
    }
  } else {
    if (B > 1) {
      if (parallel.rep) {
        cat("\nRequested parallel replicated ", K, "-fold cross-validated procedure with ", B, " replications. \n", sep="")
      } else {
        cat("\nRequested serial replicated ", K, "-fold cross-validated procedure with ", B, " replications. \n", sep="")
      }
    } else {
      cat("\nRequested single ", K, "-fold cross-validated procedure without replication. \n", sep="")
    }
  }
  cat("Variable screening: ", vs, "\n")
  if (vs) {
    cat("Variable screening technique: ", toupper(x=vstype), "\n")
  }
  cat("Cross-validation: ", cv, "\n")
  if (cv) {
    cat("Cross-validation technique: ", toupper(x=cvtype), "\n")
  }
  cat("Cross-validation criterion: ", toupper(x=cvcriterion), "\n")
  cat("Peeling criterion: ", toupper(x=peelcriterion), "\n")
  cat("Peeling percentile: ", alpha*100, "%\n")
  cat("Minimal box support: ", beta*100, "%\n")
  cat("Computation of p-values: ", pv, "\n")
  cat("Decision rule: ", ifelse(onese, yes="1SE", no="EXTREMUM"), "\n")
  if (vs) {
    cat("Parallelization of computation of variable screening:", parallel.vs, "\n")
  }
  cat("Parallelization of computation of bumps:", parallel.rep, "\n")
  if (pv) {
    cat("Parallelization of computation of p-values:", parallel.pv, "\n")
  }
  cat("\n")
  
  # Screening of informative covariates
  cv.presel.obj <- cv.presel(X=X,
                             y=y,
                             delta=delta,
                             B=B,
                             K=K,
                             vs=vs,
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
  
  success <- cv.presel.obj$success
  seed <- cv.presel.obj$seed
  
  # Survival Bump Hunting model using the algorithm
  if (!success) {
    
    if (vs) cat("Did not find any informative covariates. Exiting ... \n\n", sep="")
    
    # List of CV profiles
    CV.profiles <- list("cv.varprofiles"=NA,
                        "cv.varprofiles.mean"=NA,
                        "cv.varprofiles.se"=NA,
                        "cv.varset.opt"=NA,
                        "cv.varset.1se"=NA,
                        "cv.stepprofiles"=NA,
                        "cv.stepprofiles.mean"=NA,
                        "cv.stepprofiles.se"=NA,
                        "cv.nsteps.opt"=NA,
                        "cv.nsteps.1se"=NA)
    
    # List of CV fitted values
    CV.fit <- list("cv.maxsteps"=NA,
                   "cv.nsteps"=NA,
                   "cv.boxind"=NA,
                   "cv.boxind.size"=NA,
                   "cv.boxind.support"=NA,
                   "cv.rules"=NA,
                   "cv.screened"=NA,
                   "cv.trace"=NA,
                   "cv.sign"=NA,
                   "cv.used"=NA,
                   "cv.stats"=NA,
                   "cv.pval"=NA)
    
  } else {
    
    if (vs) cat("Successfully completed screening of covariates. \n", sep="")
    
    # Screened covariates
    CV.screened <- cv.presel.obj$varsel
    X.sel <- X[, CV.screened, drop=FALSE]
    p.sel <- length(CV.screened)
    
    # Directions of directed peeling of screened covariates
    CV.screened.sign <- cv.presel.obj$varsign
    
    # Initial box boundaries
    initcutpts <- numeric(p.sel)
    for(j in 1:p.sel){
      if (CV.screened.sign[j] == 1) {
        initcutpts[j] <- min(X.sel[,j], na.rm=TRUE)
      } else if (CV.screened.sign[j] == -1) {
        initcutpts[j] <- max(X.sel[,j], na.rm=TRUE)
      } else {
        warning("Covariate `", colnames(X.sel)[j], "` has no direction of peeling and is removed ... \n\n", sep="")
        initcutpts[j] <- NA
      }
    }
    
    # Screened covariates with valid directions of directed peeling
    w <- !is.na(initcutpts)
    CV.screened <- CV.screened[w]
    CV.screened.sign <- CV.screened.sign[w]
    initcutpts <- initcutpts[w]
    X.sel <- X[, CV.screened, drop=FALSE]
    p.sel <- length(CV.screened)
    if (vs) {
      cat("Covariates screened: \n", sep="")
      print(CV.screened)
      cat("Directions of directed peeling of screened covariates: \n", sep="")
      print(CV.screened.sign)
    }
    
    # CV profiles of screening CV criterion from all replicates
    CV.varprofiles <- cv.presel.obj$boxstat.profile
    CV.varprofiles.mean <- cv.presel.obj$boxstat.mean
    CV.varprofiles.se <- cv.presel.obj$boxstat.se
    CV.varset.opt <- cv.presel.obj$boxstat.opt
    CV.varset.1se <- cv.presel.obj$boxstat.1se
    
    # Fitting the Survival Bump Hunting model
    cat("Fitting the Survival Bump Hunting model ... \n")
    CV.box.obj <- cv.box(X=X.sel,
                         y=y,
                         delta=delta,
                         B=B,
                         K=K,
                         cv=cv,
                         cvtype=cvtype,
                         cvarg=paste("alpha=", alpha,
                                     ",beta=", beta,
                                     ",L=", Lmax,
                                     ",peelcriterion=\"", peelcriterion, "\"",
                                     ",cvcriterion=\"", cvcriterion, "\"", sep=""),
                         groups=groups,
                         decimals=decimals,
                         probval=probval,
                         timeval=timeval,
                         varsign=CV.screened.sign,
                         initcutpts=initcutpts,
                         parallel=parallel.rep,
                         conf=conf,
                         verbose=verbose,
                         seed=seed)
    
    success <- CV.box.obj$success
    seed <- CV.box.obj$seed
    
    if (!success) {
      
      cat("Could not find any bump in this dataset. Exiting ... \n\n", sep="")
      
      # List of CV profiles
      CV.profiles <- list("cv.varprofiles"=CV.varprofiles,
                          "cv.varprofiles.mean"=CV.varprofiles.mean,
                          "cv.varprofiles.se"=CV.varprofiles.se,
                          "cv.varset.opt"=CV.varset.opt,
                          "cv.varset.1se"=CV.varset.1se,
                          "cv.stepprofiles"=NA,
                          "cv.stepprofiles.mean"=NA,
                          "cv.stepprofiles.se"=NA,
                          "cv.nsteps.opt"=NA,
                          "cv.nsteps.1se"=NA)
      
      # List of CV fitted values
      CV.fit <- list("cv.maxsteps"=NA,
                     "cv.nsteps"=NA,
                     "cv.boxind"=NA,
                     "cv.boxind.size"=NA,
                     "cv.boxind.support"=NA,
                     "cv.rules"=NA,
                     "cv.screened"=NA,
                     "cv.trace"=NA,
                     "cv.sign"=NA,
                     "cv.used"=NA,
                     "cv.stats"=NA,
                     "cv.pval"=NA)
      
    } else {
      
      CV.maxsteps <- CV.box.obj$cv.maxsteps
      CV.trace.list <- CV.box.obj$cv.trace
      CV.boxind.list <- CV.box.obj$cv.boxind
      CV.boxcut.list <- CV.box.obj$cv.boxcut
      CV.support.list <- CV.box.obj$cv.support
      CV.size.list <- CV.box.obj$cv.size
      CV.bmp.lhr.list <- CV.box.obj$cv.bmp.lhr
      CV.bmp.lrt.list <- CV.box.obj$cv.bmp.lrt
      CV.bmp.cer.list <- CV.box.obj$cv.bmp.cer
      CV.grp.lhr.list <- CV.box.obj$cv.grp.lhr
      CV.grp.lrt.list <- CV.box.obj$cv.grp.lrt
      CV.grp.cer.list <- CV.box.obj$cv.grp.cer
      CV.time.bar.list <- CV.box.obj$cv.time.bar
      CV.prob.bar.list <- CV.box.obj$cv.prob.bar
      CV.max.time.bar.list <- CV.box.obj$cv.max.time.bar
      CV.min.prob.bar.list <- CV.box.obj$cv.min.prob.bar
      
      # Maximum peeling length from all replicates
      CV.maxsteps <- ceiling(mean(CV.maxsteps, na.rm=TRUE))
      
      # Box membership indicator vector of all observations at each step using the modal value or majority vote value over the replicates
      CV.boxind <- lapply.array(X=CV.boxind.list,
                                rowtrunc=CV.maxsteps,
                                MARGIN=c(1,2),
                                FUN=function(x) {
                                  return(mean(x, na.rm=TRUE) >= 0.5)
                                })
      dimnames(CV.boxind) <- list(paste("step", 0:(CV.maxsteps-1), sep=""), rownames(X.sel))
      
      # Adjusted maximum peeling length, thresholded by minimal box support, from all replicates
      CV.maxsteps <- max(which(apply(CV.boxind, 1, function(x) {length(which(x))/n >= beta})), na.rm=TRUE)
      cat("Adjusted maximum peeling length:", CV.maxsteps, "\n")
      
      if (floor(CV.maxsteps * span) == 0) {
        
        cat("Span 'span' is too small with respect to the number of peeling steps (see manual). Exiting ... \n\n", sep="")
        success <- FALSE
        
        # List of CV profiles
        CV.profiles <- list("cv.varprofiles"=CV.varprofiles,
                            "cv.varprofiles.mean"=CV.varprofiles.mean,
                            "cv.varprofiles.se"=CV.varprofiles.se,
                            "cv.varset.opt"=CV.varset.opt,
                            "cv.varset.1se"=CV.varset.1se,
                            "cv.stepprofiles"=NA,
                            "cv.stepprofiles.mean"=NA,
                            "cv.stepprofiles.se"=NA,
                            "cv.nsteps.opt"=NA,
                            "cv.nsteps.1se"=NA)
        
        # List of CV fitted values
        CV.fit <- list("cv.maxsteps"=CV.maxsteps,
                       "cv.nsteps"=NA,
                       "cv.boxind"=NA,
                       "cv.boxind.size"=NA,
                       "cv.boxind.support"=NA,
                       "cv.rules"=NA,
                       "cv.screened"=NA,
                       "cv.trace"=NA,
                       "cv.sign"=NA,
                       "cv.used"=NA,
                       "cv.stats"=NA,
                       "cv.pval"=NA)
        
      } else {
        # Box membership indicator vector of all observations at each step, adjusted to the maximum peeling length
        CV.boxind <- CV.boxind[1:CV.maxsteps,,drop=FALSE]
        
        # Get the box sample size and support based on box membership
        CV.boxind.size <- round(apply(CV.boxind, 1, sum, na.rm=TRUE), digits=0)
        names(CV.boxind.size) <- paste("step", 0:(CV.maxsteps-1), sep="")
        CV.boxind.support <- round(CV.boxind.size/n, digits=decimals)
        names(CV.boxind.support) <- paste("step", 0:(CV.maxsteps-1), sep="")
        
        # Profiles of cross-validation criterion and cross-validated peeling length from all replicates
        if (cvcriterion == "lhr") {
          if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
            CV.stepprofiles <- list2mat(list=CV.bmp.lhr.list, fill=NA, coltrunc=CV.maxsteps)
          } else if (peelcriterion == "grp") {
            CV.stepprofiles <- list2mat(list=CV.grp.lhr.list, fill=NA, coltrunc=CV.maxsteps)
          } else {
            stop("Invalid peeling criterion. Exiting ...\n\n")
          }
          colnames(CV.stepprofiles) <- paste("step", 0:(CV.maxsteps-1), sep="")
          CV.stepprofiles.mean <- apply(CV.stepprofiles, 2, mean, na.rm=TRUE)
          CV.stepprofiles.se <- apply(CV.stepprofiles, 2, sd, na.rm=TRUE)
          if (all(is.na(CV.stepprofiles.mean)) || all(is.nan(CV.stepprofiles.mean)) || is.empty(CV.stepprofiles.mean)) {
            CV.nsteps.opt <- NA
            CV.nsteps.1se <- NA
          } else {
            CV.nsteps.opt <-zeroslope(y=CV.stepprofiles.mean, x=1:CV.maxsteps, lag=lag, span=span, degree=degree, family="gaussian", minimum=FALSE)
            w <- (CV.stepprofiles.mean >= CV.stepprofiles.mean[CV.nsteps.opt]-CV.stepprofiles.se[CV.nsteps.opt])
            if (all(is.na(w)) || is.empty(w)) {
              CV.nsteps.1se <- NA
            } else {
              CV.nsteps.1se <- min(which(w))
            }
          }
        } else if (cvcriterion == "lrt") {
          if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
            CV.stepprofiles <- list2mat(list=CV.bmp.lrt.list, fill=NA, coltrunc=CV.maxsteps)
          } else if (peelcriterion == "grp") {
            CV.stepprofiles <- list2mat(list=CV.grp.lrt.list, fill=NA, coltrunc=CV.maxsteps)
          } else {
            stop("Invalid peeling criterion. Exiting ...\n\n")
          }
          colnames(CV.stepprofiles) <- paste("step", 0:(CV.maxsteps-1), sep="")
          CV.stepprofiles.mean <- apply(CV.stepprofiles, 2, mean, na.rm=TRUE)
          CV.stepprofiles.se <- apply(CV.stepprofiles, 2, sd, na.rm=TRUE)
          if (all(is.na(CV.stepprofiles.mean)) || all(is.nan(CV.stepprofiles.mean)) || is.empty(CV.stepprofiles.mean)) {
            CV.nsteps.opt <- NA
            CV.nsteps.1se <- NA
          } else {
            CV.nsteps.opt <- zeroslope(y=CV.stepprofiles.mean, x=1:CV.maxsteps, lag=lag, span=span, degree=degree, family="gaussian", minimum=FALSE)
            w <- (CV.stepprofiles.mean >= CV.stepprofiles.mean[CV.nsteps.opt]-CV.stepprofiles.se[CV.nsteps.opt])
            if (all(is.na(w)) || is.empty(w)) {
              CV.nsteps.1se <- NA
            } else {
              CV.nsteps.1se <- min(which(w))
            }
          }
        } else if (cvcriterion == "cer") {
          if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
            CV.stepprofiles <- list2mat(list=CV.bmp.cer.list, fill=NA, coltrunc=CV.maxsteps)
          } else if (peelcriterion == "grp") {
            CV.stepprofiles <- list2mat(list=CV.grp.cer.list, fill=NA, coltrunc=CV.maxsteps)
          } else {
            stop("Invalid peeling criterion. Exiting ...\n\n")
          }
          colnames(CV.stepprofiles) <- paste("step", 0:(CV.maxsteps-1), sep="")
          CV.stepprofiles.mean <- apply(CV.stepprofiles, 2, mean, na.rm=TRUE)
          CV.stepprofiles.se <- apply(CV.stepprofiles, 2, sd, na.rm=TRUE)
          if (all(is.na(CV.stepprofiles.mean)) || all(is.nan(CV.stepprofiles.mean)) || is.empty(CV.stepprofiles.mean)) {
            CV.nsteps.opt <- NA
            CV.nsteps.1se <- NA
          } else {
            CV.nsteps.opt <- zeroslope(y=CV.stepprofiles.mean, x=1:CV.maxsteps, lag=lag, span=span, degree=degree, family="gaussian", minimum=TRUE)
            w <- (CV.stepprofiles.mean <= CV.stepprofiles.mean[CV.nsteps.opt]+CV.stepprofiles.se[CV.nsteps.opt])
            if (all(is.na(w)) || is.empty(w)) {
              CV.nsteps.1se <- NA
            } else {
              CV.nsteps.1se <- min(which(w))
            }
          }
        } else {
          stop("Invalid cross-validation criterion option. Exiting ... \n\n")
        }
        
        if (onese) {
          CV.nsteps <- CV.nsteps.1se
        } else {
          CV.nsteps <- CV.nsteps.opt
        }
        
        if (is.na(CV.nsteps)) {
          
          success <- FALSE
          cat("Could not find any bump variables in this dataset. Exiting ... \n\n", sep="")
          
          # List of CV profiles
          CV.profiles <- list("cv.varprofiles"=CV.varprofiles,
                              "cv.varprofiles.mean"=CV.varprofiles.mean,
                              "cv.varprofiles.se"=CV.varprofiles.se,
                              "cv.varset.opt"=CV.varset.opt,
                              "cv.varset.1se"=CV.varset.1se,
                              "cv.stepprofiles"=NA,
                              "cv.stepprofiles.mean"=NA,
                              "cv.stepprofiles.se"=NA,
                              "cv.nsteps.opt"=NA,
                              "cv.nsteps.1se"=NA)
          
          # List of CV fitted values
          CV.fit <- list("cv.maxsteps"=CV.maxsteps,
                         "cv.nsteps"=CV.nsteps,
                         "cv.boxind"=CV.boxind,
                         "cv.boxind.size"=CV.boxind.size,
                         "cv.boxind.support"=CV.boxind.support,
                         "cv.rules"=NA,
                         "cv.screened"=NA,
                         "cv.trace"=NA,
                         "cv.sign"=NA,
                         "cv.used"=NA,
                         "cv.stats"=NA,
                         "cv.pval"=NA)
          
        } else {
          
          cat("Cross-validated peeling length:", CV.nsteps, "\n")
          
          # Box membership indicator vector of all observations at each step, adjusted to the cross-validated peeling length
          CV.boxind <- CV.boxind[1:CV.nsteps,,drop=FALSE]
          
          # Get the box sample size and support based on `CV.boxind`
          CV.boxind.size <- round(apply(CV.boxind, 1, sum), digits=0)
          names(CV.boxind.size) <- paste("step", 0:(CV.nsteps-1), sep="")
          CV.boxind.support <- round(CV.boxind.size/n, digits=decimals)
          names(CV.boxind.support) <- paste("step", 0:(CV.nsteps-1), sep="")
          
          # Trace values of covariate usage at each step for each replicate
          CV.trace <- lapply.mat(X=CV.trace.list,
                                 coltrunc=CV.nsteps,
                                 fill=NA,
                                 MARGIN=2,
                                 FUN=function(x) {
                                   vote <- table(x, useNA="no")
                                   w <- as.numeric(names(which.max(vote)))
                                   if (is.empty(w)) {
                                     return(NA)
                                   } else {
                                     return(w)
                                   }
                                 })
          
          if (is.empty(CV.trace[-1])) {
            
            success <- FALSE
            cat("Could not find any bump variables in this dataset. Exiting ... \n\n", sep="")
            
            # List of CV profiles
            CV.profiles <- list("cv.varprofiles"=CV.varprofiles,
                                "cv.varprofiles.mean"=CV.varprofiles.mean,
                                "cv.varprofiles.se"=CV.varprofiles.se,
                                "cv.varset.opt"=CV.varset.opt,
                                "cv.varset.1se"=CV.varset.1se,
                                "cv.stepprofiles"=NA,
                                "cv.stepprofiles.mean"=NA,
                                "cv.stepprofiles.se"=NA,
                                "cv.nsteps.opt"=NA,
                                "cv.nsteps.1se"=NA)
            
            # List of CV fitted values
            CV.fit <- list("cv.maxsteps"=CV.maxsteps,
                           "cv.nsteps"=CV.nsteps,
                           "cv.boxind"=CV.boxind,
                           "cv.boxind.size"=CV.boxind.size,
                           "cv.boxind.support"=CV.boxind.support,
                           "cv.rules"=NA,
                           "cv.screened"=NA,
                           "cv.trace"=NA,
                           "cv.sign"=NA,
                           "cv.used"=NA,
                           "cv.stats"=NA,
                           "cv.pval"=NA)
            
          } else {
            
            CV.trace <- sapply(X=CV.trace[-1], FUN=function(x) {x[1]})
            m <- pmatch(x=names(CV.screened)[CV.trace], table=colnames(X), nomatch=NA, duplicates.ok=TRUE)
            CV.trace <- c(0, m)
            names(CV.trace) <- paste("step", 0:(CV.nsteps-1), sep="")
            out <- c(1, which(is.na(CV.trace)))
            
            if (length(out) == length(CV.trace)) {
              
              success <- FALSE
              cat("Could not find any bump variables in this dataset. Exiting ... \n\n", sep="")
              
              # List of CV profiles
              CV.profiles <- list("cv.varprofiles"=CV.varprofiles,
                                  "cv.varprofiles.mean"=CV.varprofiles.mean,
                                  "cv.varprofiles.se"=CV.varprofiles.se,
                                  "cv.varset.opt"=CV.varset.opt,
                                  "cv.varset.1se"=CV.varset.1se,
                                  "cv.stepprofiles"=NA,
                                  "cv.stepprofiles.mean"=NA,
                                  "cv.stepprofiles.se"=NA,
                                  "cv.nsteps.opt"=NA,
                                  "cv.nsteps.1se"=NA)
              
              # List of CV fitted values
              CV.fit <- list("cv.maxsteps"=CV.maxsteps,
                             "cv.nsteps"=CV.nsteps,
                             "cv.boxind"=CV.boxind,
                             "cv.boxind.size"=CV.boxind.size,
                             "cv.boxind.support"=CV.boxind.support,
                             "cv.rules"=NA,
                             "cv.screened"=NA,
                             "cv.trace"=NA,
                             "cv.sign"=NA,
                             "cv.used"=NA,
                             "cv.stats"=NA,
                             "cv.pval"=NA)
              
            } else {
              
              success <- TRUE
              cat("Successfully completed fitting of the Survival Bump Hunting model. \n", sep="")
              
              # Covariates used in all replicates:
              CV.used <- sort(unique(CV.trace[-out]))
              names(CV.used) <- colnames(X)[CV.used]
              cat("Covariates used: \n")
              print(CV.used)
              
              # Directions of directed peeling of used covariates
              CV.sign <- CV.screened.sign[names(CV.used)]
              cat("Directions of directed peeling of used covariates: \n", sep="")
              print(CV.sign)
              
              # Box rules of used covariates at each step
              cat("Generating box rules of used covariates ...\n")
              CV.boxcut.mu <- round(lapply.array(X=CV.boxcut.list, rowtrunc=CV.nsteps, FUN=function(x){mean(x, na.rm=TRUE)}, MARGIN=1:2), digits=decimals)
              CV.boxcut.sd <- round(lapply.array(X=CV.boxcut.list, rowtrunc=CV.nsteps, FUN=function(x){sd(x, na.rm=TRUE)}, MARGIN=1:2), digits=decimals)
              rownames(CV.boxcut.mu) <- paste("step", 0:(CV.nsteps-1), sep="")
              rownames(CV.boxcut.sd) <- paste("step", 0:(CV.nsteps-1), sep="")
              colnames(CV.boxcut.mu) <- colnames(X.sel)
              colnames(CV.boxcut.sd) <- colnames(X.sel)
              CV.frame <- as.data.frame(matrix(data=NA, nrow=CV.nsteps, ncol=p.sel, dimnames=list(paste("step", 0:(CV.nsteps-1), sep=""), colnames(X.sel))))
              for (j in 1:p.sel) {
                if (CV.screened.sign[j] > 0) {
                  ss <- ">="
                } else {
                  ss <- "<="
                }
                CV.frame[, j] <- paste(paste(colnames(X.sel)[j], ss, format(x=CV.boxcut.mu[, j], digits=decimals, nsmall=decimals), sep=""),
                                       format(x=CV.boxcut.sd[, j], digits=decimals, nsmall=decimals), sep=" +/- ")
              }
              CV.rules <- list("mean"=CV.boxcut.mu[,names(CV.used),drop=FALSE],
                               "sd"=CV.boxcut.sd[,names(CV.used),drop=FALSE],
                               "frame"=CV.frame[,names(CV.used),drop=FALSE])
              
              # Box statistics at each step
              cat("Generating box statistics ...\n")
              CV.support.mu <- round(CV.boxind.support, digits=decimals)
              CV.support.sd <- round(lapply.mat(X=CV.support.list, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
              CV.size.mu <- round(CV.boxind.size, digits=0)
              CV.size.sd <- round(lapply.mat(X=CV.size.list, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
              
              if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
                CV.bmp.lhr.mu <- round(lapply.mat(X=CV.bmp.lhr.list, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.bmp.lhr.sd <- round(lapply.mat(X=CV.bmp.lhr.list, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.bmp.lrt.mu <- round(lapply.mat(X=CV.bmp.lrt.list, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.bmp.lrt.sd <- round(lapply.mat(X=CV.bmp.lrt.list, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.bmp.cer.mu <- round(lapply.mat(X=CV.bmp.cer.list, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.bmp.cer.sd <- round(lapply.mat(X=CV.bmp.cer.list, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.time.bar.mu <- round(lapply.mat(X=CV.time.bar.list, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.time.bar.sd <- round(lapply.mat(X=CV.time.bar.list, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.prob.bar.mu <- round(lapply.mat(X=CV.prob.bar.list, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.prob.bar.sd <- round(lapply.mat(X=CV.prob.bar.list, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.max.time.bar.mu <- round(lapply.mat(X=CV.max.time.bar.list, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.max.time.bar.sd <- round(lapply.mat(X=CV.max.time.bar.list, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.min.prob.bar.mu <- round(lapply.mat(X=CV.min.prob.bar.list, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.min.prob.bar.sd <- round(lapply.mat(X=CV.min.prob.bar.list, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.stats.mu <- data.frame("Support"=CV.support.mu,
                                          "Size"=CV.size.mu,
                                          "LHR"=CV.bmp.lhr.mu,
                                          "LRT"=CV.bmp.lrt.mu,
                                          "CER"=CV.bmp.cer.mu,
                                          "EFT"=CV.time.bar.mu,
                                          "EFP"=CV.prob.bar.mu,
                                          "MEFT"=CV.max.time.bar.mu,
                                          "MEFP"=CV.min.prob.bar.mu)
                rownames(CV.stats.mu) <- paste("step", 0:(CV.nsteps-1), sep="")
                colnames(CV.stats.mu) <- c("Support", "Size", "LHR", "LRT", "CER", "EFT", "EFP", "MEFT", "MEFP")
                CV.stats.sd <- data.frame("Support"=CV.support.sd,
                                          "Size"=CV.size.sd,
                                          "LHR"=CV.bmp.lhr.sd,
                                          "LRT"=CV.bmp.lrt.sd,
                                          "CER"=CV.bmp.cer.sd,
                                          "EFT"=CV.time.bar.sd,
                                          "EFP"=CV.prob.bar.sd,
                                          "MEFT"=CV.max.time.bar.sd,
                                          "MEFP"=CV.min.prob.bar.sd)
                rownames(CV.stats.sd) <- paste("step", 0:(CV.nsteps-1), sep="")
                colnames(CV.stats.sd) <- c("Support", "Size", "LHR", "LRT", "CER", "EFT", "EFP", "MEFT", "MEFP")
              } else if (peelcriterion == "grp") {
                CV.grp.lhr.mu <- round(lapply.mat(X=CV.grp.lhr.list, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.grp.lhr.sd <- round(lapply.mat(X=CV.grp.lhr.list, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.grp.cer.mu <- round(lapply.mat(X=CV.grp.cer.list, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.grp.cer.sd <- round(lapply.mat(X=CV.grp.cer.list, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
                CV.stats.mu <- data.frame("Support"=CV.support.mu,
                                          "Size"=CV.size.mu,
                                          "GLHR"=CV.grp.lhr.mu,
                                          "GCER"=CV.grp.cer.mu)
                rownames(CV.stats.mu) <- paste("step", 0:(CV.nsteps-1), sep="")
                colnames(CV.stats.mu) <- c("Support", "Size", "GLHR", "GCER")
                CV.stats.sd <- data.frame("Support"=CV.support.sd,
                                          "Size"=CV.size.sd,
                                          "GLHR"=CV.grp.lhr.sd,
                                          "GCER"=CV.grp.cer.sd)
                rownames(CV.stats.sd) <- paste("step", 0:(CV.nsteps-1), sep="")
                colnames(CV.stats.sd) <- c("Support", "Size", "GLHR", "GCER")
              } else {
                stop("Invalid peeling criterion. Exiting ...\n\n")
              }
              CV.stats <- list("mean"=CV.stats.mu, "sd"=CV.stats.sd)
              
              # Computation of p-values at each step
              if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
                obs.stat <- CV.stats$mean$LRT
              } else if (peelcriterion == "grp") {
                obs.stat <- CV.stats$mean$GLHR
              } else {
                stop("Invalid peeling criterion. Exiting ...\n\n")
              }
              CV.pval <- cv.pval(X=X.sel,
                                 y=y,
                                 delta=delta,
                                 K=K,
                                 A=A,
                                 pv=pv,
                                 cv=cv,
                                 cvtype=cvtype,
                                 decimals=decimals,
                                 varsign=CV.screened.sign,
                                 initcutpts=initcutpts,
                                 cvarg=paste("alpha=", alpha,
                                             ",beta=", beta,
                                             ",L=", CV.nsteps,
                                             ",peelcriterion=\"", peelcriterion, "\"",
                                             ",cvcriterion=\"", cvcriterion, "\"", sep=""),
                                 groups=groups,
                                 obs.stat=obs.stat,
                                 parallel=parallel.pv,
                                 conf=conf,
                                 verbose=verbose,
                                 seed=seed)
              
              # Return object 'CV.profiles'
              CV.profiles <- list("cv.varprofiles"=CV.varprofiles,
                                  "cv.varprofiles.mean"=CV.varprofiles.mean,
                                  "cv.varprofiles.se"=CV.varprofiles.se,
                                  "cv.varset.opt"=CV.varset.opt,
                                  "cv.varset.1se"=CV.varset.1se,
                                  "cv.stepprofiles"=CV.stepprofiles,
                                  "cv.stepprofiles.mean"=CV.stepprofiles.mean,
                                  "cv.stepprofiles.se"=CV.stepprofiles.se,
                                  "cv.nsteps.opt"=CV.nsteps.opt,
                                  "cv.nsteps.1se"=CV.nsteps.1se)
              
              # List of CV fitted values
              CV.fit <- list("cv.maxsteps"=CV.maxsteps,
                             "cv.nsteps"=CV.nsteps,
                             "cv.boxind"=CV.boxind,
                             "cv.boxind.size"=CV.boxind.size,
                             "cv.boxind.support"=CV.boxind.support,
                             "cv.rules"=CV.rules,
                             "cv.screened"=CV.screened,
                             "cv.trace"=CV.trace,
                             "cv.sign"=CV.sign,
                             "cv.used"=CV.used,
                             "cv.stats"=CV.stats,
                             "cv.pval"=CV.pval)
            }
          }
        }
      }
    }
  }
  
  # Returning argument 'cvarg' of algorithm parameters in a list format
  cvarg <- list("alpha"=alpha,
                "beta"=beta,
                "peelcriterion"=peelcriterion,
                "cvcriterion"=cvcriterion)
  
  # Returning the final 'sbh' object
  return(structure(list("X"=X,
                        "y"=y,
                        "groups"=groups,
                        "delta"=delta,
                        "B"=B,
                        "K"=K,
                        "A"=A,
                        "vs"=vs,
                        "vstype"=vstype,
                        "vsarg"=cv.presel.obj$vsarg,
                        "cv"=cv,
                        "cvtype"=cvtype,
                        "cvarg"=cvarg,
                        "pv"=pv,
                        "control"=control,
                        "cvprofiles"=CV.profiles,
                        "cvfit"=CV.fit,
                        "success"=success,
                        "seed"=seed),
                   class="sbh"))
}
#===============================================================================================================================#




#===============================================================================================================================#
# 2. END-USER FUNCTION FOR CONTROLLING SBH SECONDARY PARAMETERS
#===============================================================================================================================#

#===============================================================================================================================#
#===============#
#Usage         :
#===============#
#                   sbh.control(vscons = 0.5, 
#                               decimals = 2, 
#                               onese = FALSE, 
#                               probval = NULL, 
#                               timeval = NULL, 
#                               lag = 2, 
#                               span = 0.10, 
#                               degree = 2)
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

sbh.control <- function(vscons = 0.5, 
                        decimals = 2, 
                        onese = FALSE, 
                        probval = NULL, 
                        timeval = NULL, 
                        lag = 2, 
                        span = 0.10, 
                        degree = 2) {
  
  stopifnot (length(decimals) == 1L, !is.na(decimals), as.integer(decimals) >= 1L, 
             length(onese) == 1L, !is.na(onese), is.logical(onese), 
             length(lag) == 1L, !is.na(lag), as.integer(lag) >= 1,
             length(span) == 1L, !is.na(span), span >= 0, span <= 1, 
             length(degree) == 1L, !is.na(degree), as.integer(degree) >= 1)
  
  list("vscons"=vscons, 
       "decimals"=decimals, 
       "onese"=onese, 
       "probval"=probval, 
       "timeval"=timeval, 
       "lag"=lag, 
       "span"=span,
       "degree"=degree)
  
}
#===============================================================================================================================#




#===============================================================================================================================#
# 3. END-USER FUNCTION FOR NEWS
#===============================================================================================================================#

#===============================================================================================================================#
#===============#
#Usage         :
#===============#
#                   PRIMsrc.news(...)
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

PRIMsrc.news <- function(...) {
  
  newsfile <- file.path(system.file(package="PRIMsrc"), "NEWS")
  file.show(newsfile)
  
}
#===============================================================================================================================#




#===============================================================================================================================#
# 4. END-USER S3-GENERIC FUNCTIONS FOR SUMMARY, PRINT, PLOT AND PREDICTION
#===============================================================================================================================#

#===============================================================================================================================#
#===============#
#Usage         :
#===============#
#                   summary(object, ...)
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

summary.sbh <- function(object, ...) {
  
  if (!inherits(object, 'sbh'))
    stop("Argument `object` must be an object of class 'sbh'. Exiting ... \n\n")
  
  cat("\nS3-class object: '", attr(x=object, "class"), "' \n\n")
  if (!object$cv) {
    if (object$B > 1) {
      cat("Replicated procedure with ", object$B, " replications without cross-validation. \n\n", sep="")
    } else {
      cat("Single procedure without cross-validation and replication. \n\n", sep="")
    }
  } else {
    if (object$B > 1) {
      cat("Replicated ", object$K, "-fold cross-validated procedure with ", object$B, " replications. \n\n", sep="")
    } else {
      cat("Single ", object$K, "-fold cross-validated procedure without replication. \n\n", sep="")
    }
  }
  cat("VARIABLE SCREENING:\n")
  cat("\t Variable screening: ", object$vs, "\n")
  if (object$vs) {
    cat("\t Variable screening technique: ", toupper(x=object$vstype), "\n\n")
  }
  cat("CROSS-VALIDATION:\n")
  cat("\t Cross-validation: ", object$cv, "\n")
  if (object$cv) {
    cat("\t Cross-validation technique: ", toupper(x=object$cvtype), "\n")
  }
  cat("ALGORITHM PARAMETERS:\n")
  cat("\t Cross-validation criterion: ", toupper(x=object$cvarg$cvcriterion), "\n")
  cat("\t Peeling criterion: ", toupper(x=object$cvarg$peelcriterion), "\n")
  cat("\t Peeling percentile: ", object$cvarg$alpha*100, "%\n")
  cat("\t Minimal box support: ", object$cvarg$beta*100, "%\n")
  cat("REPORTING:\n")
  cat("\t Decision rule: ", ifelse(object$control$onese, yes="1SE", no="EXTREMUM"), "\n")
  cat("\t Number of decimals: ", object$control$decimals, "\n")
  cat("\t Computation of p-values: ", object$pv, "\n\n")
  
  invisible()
  
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
#Usage         :
#===============#
#                   print(x, ...)
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

print.sbh <- function(x, ...) {
  
  if (!inherits(x, 'sbh'))
    stop("Argument `x` must be an object of class 'sbh'. Exiting ... \n\n")
  
  if (x$success) {
    
    cat("\n")
    cat("Screened (pre-selected) covariates:\n")
    print(x$cvfit$cv.screened)
    cat("\n")
    
    cat("Used (selected) covariates:\n")
    print(x$cvfit$cv.used)
    cat("\n")
    
    cat("Maximum number of peeling steps (counting step #0):\n")
    print(x$cvfit$cv.maxsteps)
    cat("\n")
    
    out <- x$cvfit$cv.nsteps
    names(out) <- NULL
    cat("Optimal number of peeling steps (counting step #0):\n")
    print(out)
    cat("\n")
    
    cat("Trace of covariate usage at each peeling step:\n")
    print(x$cvfit$cv.trace)
    cat("\n")
    
    cat("P-values at each peeling step:\n")
    print(x$cvfit$cv.pval$pval, quote = FALSE)
    cat("\n")
    
    cat("Decision rules for the used covariates (columns) at each peeling step (rows):\n")
    print(x$cvfit$cv.rules$frame, quote = FALSE)
    cat("\n")
    
    cat("Box endpoint quantities of interest (columns) at each peeling step (rows):\n")
    print(x$cvfit$cv.stats$mean)
    cat("\n")
    
    cat("Individual observation box membership indicator (columns) at each peeling step (rows):\n")
    print(x$cvfit$cv.boxind)
    cat("\n")
    
    cat("Variable selection parameters:\n")
    print(x$vsarg)
    cat("\n")
    
    cat("ALGORITHM parameters:\n")
    print(x$cvarg)
    cat("\n")
    
  } else {
    
    cat("Either the covariate screening or the Survival Bump Hunting modeling failed for this dataset.\n
         So, there is nothing to print here.\n")
    
  }
  
  invisible()
  
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    plot(x,
#                         main="Scatter Plot",
#                         proj=x$cvfit$cv.used[c(1,2)], 
#                         steps=x$cvfit$cv.nsteps,
#                         pch=16, 
#                         cex=0.5, 
#                         col=c(1,2),
#                         boxes=TRUE,
#                         asp=NA,
#                         col.box=rep(2,length(steps)), 
#                         lty.box=rep(2,length(steps)), 
#                         lwd.box=rep(1,length(steps)), 
#                         add.caption.box=boxes,
#                         text.caption.box=paste("Step: ", steps, sep=""), 
#                         pch.group=c(1,1),
#                         cex.group=c(1,1),
#                         col.group=c(3,4),
#                         add.caption.group=ifelse(test=(x$cvarg$peelcriterion == "grp"), 
#                                                  yes=TRUE, 
#                                                  no=FALSE),
#                         text.caption.group=levels(x$groups),
#                         device=NULL, 
#                         file="Scatter Plot", 
#                         path=getwd(),
#                         horizontal=FALSE, 
#                         width=5, 
#                         height=5, ...)
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

plot.sbh <- function(x,
                     main="Scatter Plot",
                     proj=x$cvfit$cv.used[c(1,2)],
                     steps=x$cvfit$cv.nsteps,
                     pch=16,
                     cex=0.5,
                     col=c(1,2), 
                     boxes=TRUE,
                     asp=NA,
                     col.box=rep(2,length(steps)),
                     lty.box=rep(2,length(steps)),
                     lwd.box=rep(1,length(steps)),
                     add.caption.box=boxes,
                     text.caption.box=paste("Step: ", steps, sep=""), 
                     pch.group=c(1,1),
                     cex.group=c(1,1),
                     col.group=c(3,4),
                     add.caption.group=ifelse(test=(x$cvarg$peelcriterion == "grp"), 
                                              yes=TRUE, 
                                              no=FALSE),
                     text.caption.group=levels(x$groups), 
                     device=NULL,
                     file="Scatter Plot",
                     path=getwd(),
                     horizontal=FALSE,
                     width=5,
                     height=5, ...) {
  
  if (!inherits(x, 'sbh'))
    stop("Argument `x` must be an object of class 'sbh'. Exiting ... \n\n")
  
  if (x$success) {
    
    if (length(proj) <= 1) {
      stop("Argument `proj` must be a vector of length 2. Exiting ... \n\n")
    } else {
      if (all(proj %in% (1:ncol(x$X)))) {
        if (length(x$cvfit$cv.used) == 1) {
          toadd <- setdiff(proj, x$cvfit$cv.used)
          toplot <- c(toadd, x$cvfit$cv.used)
          toplot <- sort(toplot, decreasing = FALSE)
          cat("Warning: Added dimension ", toadd, " to the used covariates of 'sbh' object `x`. \n", sep="")
        } else {
          if (all(proj %in% x$cvfit$cv.used)) {
            toadd <- NULL
            toplot <- proj
          } else {
            stop("Argument `proj` must be a equal or a subset of the used covariates of 'sbh' object `x`. Exiting ... \n\n")
          }
        }
      } else {
        stop("Argument `proj` must be a equal or a subset of the covariates of 'sbh' object `x`. Exiting ... \n\n")
      }
    }
    
    scatterplot <- function(object,
                            main,
                            toadd, toplot,
                            steps, 
                            pch, cex, col, 
                            boxes, asp,
                            col.box, lty.box, lwd.box, 
                            add.caption.box, text.caption.box, 
                            pch.group, cex.group, col.group, 
                            add.caption.group, text.caption.group, ...) {
      
      if (!is.null(main)) {
        par(mfrow=c(1, 1), oma=c(0, 0, 3, 0), mar=c(2.5, 2.5, 4.0, 1.5), mgp=c(1.5, 0.5, 0))
      } else {
        par(mfrow=c(1, 1), oma=c(0, 0, 0, 0), mar=c(2.5, 2.5, 4.0, 1.5), mgp=c(1.5, 0.5, 0))
      }
      
      peelcriterion <- object$cvarg$peelcriterion
      L <- length(steps)
      y <- object$y
      x <- object$X[,toplot,drop=FALSE]
      x.names <- colnames(x)
      plot(x=x, main=NULL, xlab=x.names[1], ylab=x.names[2], type="n", axes=FALSE, asp=asp, frame.plot=FALSE)
      axis(side=1, at=pretty(range(x[,1])), col=1, col.axis=1, cex.axis=1, line=0)
      axis(side=2, at=pretty(range(x[,2])), col=1, col.axis=1, cex.axis=1, line=0)
      if (peelcriterion == "grp") {
        groups <- as.factor(x=object$groups)
        groups.lev <- levels(groups)
        groups.ng <- nlevels(groups)
        groups.def <- vector(mode="list", length=groups.ng)
        for (g in 1:groups.ng) {
          groups.def[[g]] <- which(groups == groups.lev[g])
        }
        for (g in 1:groups.ng) {
          points(x[groups.def[[g]],], type="p", pch=pch.group[g], cex=cex.group[g], col=col.group[g], ...)
        }
      }
      points(x=x, type="p", pch=pch, cex=cex, col=col[1], ...)
      w <- object$cvfit$cv.boxind[steps[L],]
      points(x=x[w,], type="p", pch=pch, cex=cex, col=col[2], ...)
      if (boxes) {
        x.range <- apply(X=x, MARGIN=2, FUN=range, na.rm=TRUE)
        if (is.null(toadd)) {
          boxcut <- object$cvfit$cv.rules$mean[,x.names,drop=FALSE]
          varsign <- object$cvfit$cv.sign[x.names]
        } else {
          boxcut <- object$cvfit$cv.rules$mean[,1,drop=FALSE]
          varsign <- object$cvfit$cv.sign[1]
          if (toadd > object$cvfit$cv.used) {
            boxcut <- cbind(boxcut, x.range[1,toadd])
            varsign <- c(varsign, 1)
          } else {
            boxcut <- cbind(x.range[1,toadd], boxcut)
            varsign <- c(1, varsign)
          }
        }
        vertices <- vector(mode="list", length=L)
        for (l in 1:L) {
          i <- steps[l]
          vertices[[l]] <- matrix(data=NA, nrow=2, ncol=2, dimnames=list(c("LB","UB"), x.names))
          for (j in 1:2) {
            vertices[[l]][1,j] <- ifelse(test=(varsign[j] > 0),
                                         yes=max(x.range[1,j], boxcut[i,j]),
                                         no=min(x.range[1,j], boxcut[i,j]))
            vertices[[l]][2,j] <- ifelse(test=(varsign[j] < 0),
                                         yes=min(x.range[2,j], boxcut[i,j]),
                                         no=max(x.range[2,j], boxcut[i,j]))
          }
        }
        for (l in 1:L) {
          rect(vertices[[l]][1,1], vertices[[l]][1,2], vertices[[l]][2,1], vertices[[l]][2,2],
               border=col.box, col=NA, lty=lty.box, lwd=lwd.box)
        }
      }
      if (add.caption.box) {
        legend(x="top", inset=-0.18, legend=c("outbox","inbox"), cex=cex, pch=pch, col=col, xpd=TRUE)
        legend(x="top", inset=-0.08, legend=text.caption.box, cex=cex, col=col.box, lty=lty.box, lwd=lwd.box, xpd=TRUE)
      }
      if (add.caption.group) {
        legend(x="right", inset=0, legend=text.caption.group, cex=cex, pch=pch.group, pt.cex=cex.group, col=col.group, xpd=TRUE)
      }
      if (!is.null(main)) {
        mtext(text=main, cex=1, side=3, outer=TRUE)
      }
    }
    
    if (is.null(device)) {
      scatterplot(object=x,
                  main=main,
                  toadd=toadd, toplot=toplot, 
                  steps=steps,
                  pch=pch, cex=cex, col=col,
                  boxes=boxes, asp=asp,
                  col.box=col.box, lty.box=lty.box, lwd.box=lwd.box, 
                  add.caption.box=add.caption.box, text.caption.box=text.caption.box, 
                  pch.group=pch.group, cex.group=cex.group, col.group=col.group,
                  add.caption.group=add.caption.group, text.caption.group=text.caption.group)
    } else if (device == "PS") {
      path <- normalizePath(path=path, winslash="\\", mustWork=FALSE)
      file <- paste(file, ".ps", sep="")
      cat("\nOUTPUT: \n")
      cat("Filename : ", file, "\n")
      cat("Directory: ", path, "\n")
      postscript(file=file.path(path, file, fsep=.Platform$file.sep), width=width, height=height, onefile=TRUE, horizontal=horizontal)
      scatterplot(object=x,
                  main=main,
                  toadd=toadd, toplot=toplot, 
                  steps=steps,
                  pch=pch, cex=cex, col=col,
                  boxes=boxes, asp=asp,
                  col.box=col.box, lty.box=lty.box, lwd.box=lwd.box, 
                  add.caption.box=add.caption.box, text.caption.box=text.caption.box, 
                  pch.group=pch.group, cex.group=cex.group, col.group=col.group,
                  add.caption.group=add.caption.group, text.caption.group=text.caption.group)
      dev.off()
    } else if (device == "PDF") {
      path <- normalizePath(path=path, winslash="\\", mustWork=FALSE)
      file <- paste(file, ".pdf", sep="")
      cat("\nOUTPUT: \n")
      cat("Filename : ", file, "\n")
      cat("Directory: ", path, "\n")
      pdf(file=file.path(path, file, fsep=.Platform$file.sep), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
      scatterplot(object=x,
                  main=main,
                  toadd=toadd, toplot=toplot, 
                  steps=steps,
                  pch=pch, cex=cex, col=col,
                  boxes=boxes, asp=asp,
                  col.box=col.box, lty.box=lty.box, lwd.box=lwd.box, 
                  add.caption.box=add.caption.box, text.caption.box=text.caption.box, 
                  pch.group=pch.group, cex.group=cex.group, col.group=col.group,
                  add.caption.group=add.caption.group, text.caption.group=text.caption.group)
      dev.off()
    } else {
      stop("Currently allowed display devices are \"PS\" (Postscript) or \"PDF\" (Portable Document Format. Exiting ... \n\n) \n")
    }
    
  } else {
    
    cat("Either the covariate screening or the Survival Bump Hunting modeling failed for this dataset.\n
         So, there is nothing to plot here.\n")
    
  }
  
  invisible()
  
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
#Usage         :
#===============#
#                   predict(object,
#                           newdata,
#                           steps=1:object$cvfit$cv.nsteps,
#                           groups=NULL,
#                           na.action=na.omit, ...)
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

predict.sbh <- function (object,
                         newdata,
                         steps=1:object$cvfit$cv.nsteps,
                         groups=NULL,
                         na.action=na.omit, ...) {
  
  if (object$success) {
    
    if (!inherits(object, 'sbh'))
      stop("Argument `object` must be an object of class 'sbh'. Exiting ... \n\n")
    
    if (!is.matrix(newdata)) {
      if (is.data.frame(newdata)) {
        if (!all(sapply(newdata, is.numeric))) {
          stop("Input data.frame 'newdata' is not all numeric.\n\n")
        } else { 
          cat("Warning: Input data.frame `newdata` was corced to numeric matrix.\n\n", sep="")
          newdata <- as.matrix(newdata)
          mode(newdata) <- "numeric"
        }
      }
    }
    
    CV.used <- object$cvfit$cv.used
    if (all(CV.used %in% (1:ncol(newdata)))) {
      x <- newdata[,CV.used,drop=FALSE]
      n <- nrow(x)
      p <- ncol(x)
      x.names <- colnames(x)
      x.range <- apply(X=x, MARGIN=2, FUN=range)
    } else {
      stop("The used covariates of 'sbh' `object` must be equal or a subset of the `newdata` covariates. Exiting ... \n\n")
    }
    
    varnames <- colnames(object$X)
    L <- length(steps)
    boxcut <- object$cvfit$cv.rules$mean[,varnames[CV.used],drop=FALSE]
    
    CV.rules <- list("mean"=object$cvfit$cv.rules$mean[steps,varnames[CV.used],drop=FALSE],
                     "sd"=object$cvfit$cv.rules$sd[steps,varnames[CV.used],drop=FALSE],
                     "frame"=object$cvfit$cv.rules$frame[steps,varnames[CV.used],drop=FALSE])
    
    CV.sign <- object$cvfit$cv.sign[varnames[CV.used]]
    
    CV.boxind <- matrix(NA, nrow=L, ncol=n, dimnames=list(paste("step ", steps, sep=""), rownames(x)))
    for (l in 1:L) {
      i <- steps[l]
      boxcutsign <- boxcut[i,] * CV.sign
      x.cut <- t(t(x) * CV.sign)
      x.ind <- t(t(x.cut) >= boxcutsign)
      CV.boxind[l,] <- (rowMeans(x.ind) == 1)  # Set as TRUE which observations are inside the box boundaries for all axes directions
    }
    
    # Get the box sample size and support based on `CV.boxind`
    CV.boxind.size <- round(apply(CV.boxind, 1, sum), digits=0)
    names(CV.boxind.size) <- paste("step", steps, sep="")
    CV.boxind.support <- round(CV.boxind.size/n, digits=object$control$decimals)
    names(CV.boxind.support) <- paste("step", steps, sep="")
    
    # Return object 'CV.screened'
    CV.screened <- NA
    
    # Return object 'CV.trace'
    CV.trace <- NA
    
    # Return object 'CV.stats'
    CV.stats <- list("mean"=NA, 
                     "sd"=NA)
    
    # Return object 'CV.pval'
    CV.pval <- list("pval"=NA, 
                    "seed"=NA)
    
    # Return object 'CV.profiles'
    CV.profiles <- list("cv.varprofiles"=NA,
                        "cv.varprofiles.mean"=NA,
                        "cv.varprofiles.se"=NA,
                        "cv.varset.opt"=NA,
                        "cv.varset.1se"=NA,
                        "cv.stepprofiles"=NA,
                        "cv.stepprofiles.mean"=NA,
                        "cv.stepprofiles.se"=NA,
                        "cv.nsteps.opt"=NA,
                        "cv.nsteps.1se"=NA)
    
    # List of CV fitted values
    CV.fit <- list("cv.maxsteps"=object$cvfit$cv.maxsteps,
                   "cv.nsteps"=object$cvfit$cv.nsteps,
                   "cv.boxind"=CV.boxind,
                   "cv.boxind.size"=CV.boxind.size,
                   "cv.boxind.support"=CV.boxind.support,
                   "cv.rules"=CV.rules,
                   "cv.sign"=CV.sign,
                   "cv.used"=CV.used,
                   "cv.screened"=CV.screened,
                   "cv.trace"=CV.trace,
                   "cv.stats"=CV.stats,
                   "cv.pval"=CV.pval)
    
    return(structure(list("X"=newdata,
                          "y"=NA,
                          "delta"=NA,
                          "B"=object$B,
                          "K"=object$K,
                          "A"=object$A,
                          "vs"=object$vs,
                          "vstype"=object$vstype,
                          "vsarg"=object$cv.presel.obj$vsarg,
                          "cv"=object$cv,
                          "cvtype"=object$cvtype,
                          "cvarg"=object$cvarg,
                          "groups"=groups,
                          "pv"=object$pv,
                          "control"=object$control,
                          "cvprofiles"=CV.profiles,
                          "cvfit"=CV.fit,
                          "success"=object$success,
                          "seed"=object$seed),
                     class="sbh"))
    
  } else {
    
    cat("Either the covariate screening or the Survival Bump Hunting modeling failed for this dataset.\n")
    
    return(structure(list(NA), class="sbh"))
    
  }
  
}
#===============================================================================================================================#




#===============================================================================================================================#
# 5. END-USER PLOTTING FUNCTIONS FOR MODEL VALIDATION AND VISUALIZATION OF RESULTS
#===============================================================================================================================#

#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    plot_profile(object,
#                                 main="Profile Plots",
#                                 xlim=NULL, 
#                                 ylim=NULL,
#                                 add.sd=TRUE, 
#                                 add.profiles=TRUE,
#                                 add.caption=TRUE, 
#                                 text.caption=c("Mean","Std. Error"),
#                                 pch=20, 
#                                 col=1, 
#                                 lty=1, 
#                                 lwd=0.5, 
#                                 cex=0.5,
#                                 device=NULL, 
#                                 file="Profile Plots", 
#                                 path=getwd(),
#                                 horizontal=FALSE, 
#                                 width=8.5, 
#                                 height=5.0, ...) {
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

plot_profile <- function(object,
                         main="Profile Plots",
                         xlim=NULL,
                         ylim=NULL,
                         add.sd=TRUE,
                         add.profiles=TRUE,
                         add.caption=TRUE,
                         text.caption=c("Mean","Std. Error"),
                         pch=20,
                         col=1,
                         lty=1,
                         lwd=0.5,
                         cex=0.5,
                         device=NULL,
                         file="Profile Plots",
                         path=getwd(),
                         horizontal=FALSE,
                         width=8.5,
                         height=5.0, ...) {
  
  if (!inherits(object, 'sbh'))
    stop("Argument `object` must be an object of class 'sbh'. Exiting ... \n\n")
  
  if (object$success) {
    
    profileplot <- function(object, main, xlim, ylim,
                            add.sd, add.caption, text.caption, add.profiles,
                            pch, col, lty, lwd, cex, ...) {
      
      peelcriterion <- object$cvarg$peelcriterion
      cvcriterion <- object$cvarg$cvcriterion
      if (cvcriterion == "lhr") {
        txt <- "LHR"
      } else if (cvcriterion == "lrt") {
        txt <- "LRT"
      } else if (cvcriterion == "cer") {
        if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
          txt <- "CER"
        } else if (peelcriterion == "grp") {
          txt <- "GCER"
        } else {
          stop("Invalid peeling criterion. Exiting ...\n\n")
        }
      } else {
        stop("Invalid CV criterion. Exiting ... \n\n")
      }
      
      if ((object$vs && object$vstype == "prsp") && (ncol(object$cvprofiles$cv.varprofiles) > 1)) {
        if (!is.null(main)) {
          par(mfrow=c(2, 1), oma=c(0, 0, 3, 0), mar=c(4.0, 3.0, 4.0, 3.0), mgp=c(1.5, 0.5, 0))
        } else {
          par(mfrow=c(2, 1), oma=c(0, 0, 0, 0), mar=c(4.0, 3.0, 4.0, 3.0), mgp=c(1.5, 0.5, 0))
        }
        varprofiles <- object$cvprofiles$cv.varprofiles
        varprofiles.mean <- object$cvprofiles$cv.varprofiles.mean
        varprofiles.se <- object$cvprofiles$cv.varprofiles.se
        if (object$control$onese) {
          cv.varset <- object$cvprofiles$cv.varset.1se
        } else {
          cv.varset <- object$cvprofiles$cv.varset.opt
        }
        msize <- object$vsarg$msize
        
        if (is.null(xlim)) {
          xlimv <- range(msize, na.rm=TRUE)
        } else {
          xlimv <- xlim
        }
        if (is.null(ylim)) {
          ylimv <- range(varprofiles, na.rm=TRUE)
        } else {
          ylimv <- ylim
        }
        if (add.profiles) {
          matplot(t(varprofiles), type="l", axes=FALSE,
                  xlab="", ylab="", main="",
                  xlim=xlimv, ylim=ylimv, pch=pch, lty=1, lwd=lwd/4, cex=cex/4, ...)
          par(new=TRUE)
        }
        plot(x=msize, y=varprofiles.mean, type="b", axes=FALSE,
             main="CV Tuning Profiles of Variable Screening Size", cex.main=0.8,
             xlab="Cardinals of Top-Ranked Variables Subsets", ylab=paste(txt, " Mean Profiles", sep=""),
             xlim=xlimv, ylim=ylimv,
             pch=pch, col=col, lty=lty, lwd=lwd, cex=cex, ...)
        axis(side=1, pos=min(ylimv), at=msize, labels=msize, col=1, cex.axis=1, cex.axis=1, line=NA)
        axis(side=2, pos=min(xlimv), at=pretty(ylimv), col=1, col.axis=1, cex.axis=1, line=NA)
        segments(x0=msize[cv.varset], y0=min(ylimv),
                 x1=msize[cv.varset], y1=varprofiles.mean[cv.varset],
                 col=col, lty=2, lwd=lwd)
        if (add.sd) {
          arrows(msize, varprofiles.mean,
                 msize, varprofiles.mean - varprofiles.se,
                 length=0.03, angle=90, code=2, col=col, lwd=lwd)
          arrows(msize, varprofiles.mean,
                 msize, varprofiles.mean + varprofiles.se,
                 length=0.03, angle=90, code=2, col=col, lwd=lwd)
        }
        if (add.caption) {
          legend(x="top", xpd=TRUE, inset=0, legend=text.caption, pch=pch, col=col, lty=lty, lwd=lwd, cex=cex, pt.cex=cex/2)
        }
      } else {
        if (!is.null(main)) {
          par(mfrow=c(1, 1), oma=c(0, 0, 3, 0), mar=c(4.0, 3.0, 4.0, 3.0), mgp=c(1.5, 0.5, 0))
        } else {
          par(mfrow=c(1, 1), oma=c(0, 0, 0, 0), mar=c(4.0, 3.0, 4.0, 3.0), mgp=c(1.5, 0.5, 0))
        }
      }
      
      if (object$control$onese) {
        cv.nsteps <- object$cvprofiles$cv.nsteps.1se
      } else {
        cv.nsteps <- object$cvprofiles$cv.nsteps.opt
      }
      stepprofiles <- object$cvprofiles$cv.stepprofiles
      stepprofiles.mean <- object$cvprofiles$cv.stepprofiles.mean
      stepprofiles.se <- object$cvprofiles$cv.stepprofiles.se
      Lm <- object$cvfit$cv.maxsteps
      
      if (is.null(xlim)) {
        xlims <- c(0, Lm)
      } else {
        xlims <- xlim
      }
      if (is.null(ylim)) {
        ylims <- range(stepprofiles, na.rm=TRUE)
      } else{
        ylims <- ylim
      }
      if (add.profiles) {
        matplot(t(stepprofiles), type="l", axes=FALSE,
                xlab="", ylab="", main="",
                xlim=xlims, ylim=ylims, pch=pch, lty=1, lwd=lwd/4, cex=cex/4, ...)
        par(new=TRUE)
      }
      plot(0:(Lm-1), stepprofiles.mean, type="b", axes=FALSE,
           main="CV Tuning Profiles of Peeling Length", cex.main=0.8,
           xlab="Peeling Steps", ylab=paste(txt, " Mean Profiles", sep=""),
           xlim=xlims, ylim=ylims,
           pch=pch, col=col, lty=lty, lwd=lwd, cex=cex, ...)
      axis(side=1, pos=min(ylims), at=0:(Lm-1), labels=0:(Lm-1), col=1, cex.axis=1, line=NA)
      axis(side=2, pos=min(xlims), at=pretty(ylims), col=1, cex.axis=1, line=NA)
      segments(x0=(0:(Lm-1))[cv.nsteps], y0=min(ylims),
               x1=(0:(Lm-1))[cv.nsteps], y1=stepprofiles.mean[cv.nsteps],
               col=col, lty=2, lwd=lwd)
      if (add.sd) {
        arrows(0:(Lm-1), stepprofiles.mean,
               0:(Lm-1), stepprofiles.mean - stepprofiles.se,
               length=0.03, angle=90, code=2, col=col, lwd=lwd)
        arrows(0:(Lm-1), stepprofiles.mean,
               0:(Lm-1), stepprofiles.mean + stepprofiles.se,
               length=0.03, angle=90, code=2, col=col, lwd=lwd)
      }
      if (add.caption) {
        legend("top", xpd=TRUE, inset=0, legend=c("Sample Mean", "Std. Error"),
               pch=pch, col=col, lty=lty, lwd=lwd, cex=cex, pt.cex=cex/2)
      }
      if (!is.null(main)) {
        mtext(text=main, cex=1, col=1, side=3, line=1, outer=TRUE)
      }
    }
    
    if (is.null(device)) {
      cat("Device: ",  dev.cur(), "\n")
      profileplot(object=object, main=main, xlim=xlim, ylim=ylim,
                  add.sd=add.sd, add.caption=add.caption, text.caption=text.caption, add.profiles=add.profiles,
                  pch=pch, col=col, lty=lty, lwd=lwd, cex=cex)
    } else if (device == "PS") {
      path <- normalizePath(path=path, winslash="\\", mustWork=FALSE)
      file <- paste(file, ".ps", sep="")
      cat("\nOUTPUT: \n")
      cat("Filename : ", file, "\n")
      cat("Directory: ", path, "\n")
      postscript(file=file.path(path, file, fsep=.Platform$file.sep), width=width, height=height, onefile=TRUE, horizontal=horizontal)
      cat("Device: ",  dev.cur(), "\n")
      profileplot(object=object, main=main, xlim=xlim, ylim=ylim,
                  add.sd=add.sd, add.caption=add.caption, text.caption=text.caption, add.profiles=add.profiles,
                  pch=pch, col=col, lty=lty, lwd=lwd, cex=cex)
      dev.off()
    } else if (device == "PDF") {
      path <- normalizePath(path=path, winslash="\\", mustWork=FALSE)
      file <- paste(file, ".pdf", sep="")
      cat("\nOUTPUT: \n")
      cat("Filename : ", file, "\n")
      cat("Directory: ", path, "\n")
      pdf(file=file.path(path, file, fsep=.Platform$file.sep), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
      cat("Device: ",  dev.cur(), "\n")
      profileplot(object=object, main=main, xlim=xlim, ylim=ylim,
                  add.sd=add.sd, add.caption=add.caption, text.caption=text.caption, add.profiles=add.profiles,
                  pch=pch, col=col, lty=lty, lwd=lwd, cex=cex)
      dev.off()
    } else {
      cat("Currently allowed display devices are \"PS\" (Postscript) or \"PDF\" (Portable Document Format). Exiting ... \n\n")
    }
    
  } else {
    
    cat("Either the covariate screening or the Survival Bump Hunting modeling failed for this dataset.\n
         So, there is nothing to plot here.\n")
    
  }
  
  invisible()
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    plot_traj (object,
#                               main="Trajectory Plots",
#                               toplot=object$cvfit$cv.used,
#                               range=NULL,
#                               col.cov, 
#                               lty.cov, 
#                               lwd.cov,
#                               col=1, 
#                               lty=1, 
#                               lwd=0.5, 
#                               cex=0.5,
#                               add.caption=FALSE, 
#                               text.caption=NULL,
#                               nr=NULL, 
#                               nc=NULL,
#                               device=NULL, 
#                               file="Trajectory Plots", 
#                               path=getwd())
#                               horizontal=FALSE, 
#                               width=8.5, 
#                               height=11, ...)
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

plot_traj <- function(object,
                      main="Trajectory Plots",
                      toplot=object$cvfit$cv.used,
                      range=NULL,
                      col.cov,
                      lty.cov,
                      lwd.cov,
                      col=1,
                      lty=1,
                      lwd=0.5,
                      cex=0.5,
                      add.caption=FALSE,
                      text.caption=NULL,
                      nr=NULL,
                      nc=NULL,
                      device=NULL,
                      file="Trajectory Plots",
                      path=getwd(),
                      horizontal=FALSE,
                      width=8.5,
                      height=11, ...) {
  
  if (!inherits(object, 'sbh'))
    stop("Argument `object` must be an object of class 'sbh'. Exiting ... \n\n")
  
  if (object$success) {
    
    trajplot <- function(object,
                         main,
                         toplot,
                         range,
                         col.cov, lty.cov, lwd.cov,
                         col, lty, lwd, cex, 
                         add.caption, text.caption,
                         nr, nc, ...) {
      p <- length(toplot)
      varnames <- colnames(object$X)
      peelcriterion <- object$cvarg$peelcriterion
      if (is.null(nc))
        nc <- 3
      if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
        if (is.null(nr)) {
          if (p %% nc == 0) {
            nr <- p%/%nc + 2
          } else {
            nr <- p%/%nc + 3
          }
        }
      } else if (peelcriterion == "grp") {
        if (is.null(nr)) {
          if (p %% nc == 0) {
            nr <- p%/%nc + 1
          } else {
            nr <- p%/%nc + 2
          }
        }
      } else {
        stop("Invalid peeling criterion. Exiting ...\n\n")
      }
      if (is.null(range)) {
        if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
          range <- list("lhr"=range(0, object$cvfit$cv.stats$mean$LHR, na.rm=TRUE), 
                        "lrt"=range(0, object$cvfit$cv.stats$mean$LRT, na.rm=TRUE), 
                        "cer"=range(0, object$cvfit$cv.stats$mean$CER, na.rm=TRUE))
        } else if (peelcriterion == "grp") {
          range <- list("glhr"=range(0, object$cvfit$cv.stats$mean$GLHR, na.rm=TRUE),
                        "gcer"=range(0, object$cvfit$cv.stats$mean$GCER, na.rm=TRUE))
        } else {
          stop("Invalid peeling criterion. Exiting ...\n\n")
        }
      }
      if (missing(col.cov)) {
        col.cov <- 2:(p+1)
      }
      if (missing(lty.cov)) {
        lty.cov <- rep(1,p)
      }
      if (missing(lwd.cov)) {
        lwd.cov <- rep(1,p)
      }
      if (!is.null(main)) {
        par(mfrow=c(nr, nc), oma=c(0, 0, 3, 0), mar=c(2.5, 2.5, 2.0, 1.5), mgp=c(1.5, 0.5, 0))
      } else {
        par(mfrow=c(nr, nc), oma=c(0, 0, 0, 0), mar=c(2.5, 2.5, 2.0, 1.5), mgp=c(1.5, 0.5, 0))
      }
      for (j in 1:p) {
        plot(x=object$cvfit$cv.stats$mean$Support,
             y=object$cvfit$cv.rules$mean[,varnames[toplot[j]]],
             type='s', col=col.cov[j], lty=lty.cov[j], lwd=lwd.cov[j],
             main=paste(varnames[toplot[j]], " covariate", sep=""), cex.main=cex,
             xlim=range(0, 1),
             ylim=range(object$X[,toplot[j]], na.rm=TRUE),
             xlab=expression(paste("Box Support (", beta, ")", sep="")),
             ylab="Covariate Range", ...)
      }
      if (add.caption)
        legend(x="bottomleft", inset=0.01, legend=text.caption, cex=cex)
      
      if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
        par(mfg=c(nr-1, 1))
        plot(object$cvfit$cv.stats$mean$Support,
             object$cvfit$cv.stats$mean$Size,
             type='s', col=col, lty=lty, lwd=lwd,
             main="Box Sample Size", cex.main=cex,
             xlim=range(0, 1),
             ylim=range(0, object$cvfit$cv.stats$mean$Size),
             xlab=expression(paste("Box Support (", beta, ")", sep="")),
             ylab=expression(paste("Size (", italic(n), ")", sep="")), ...)
        if (add.caption)
          legend(x="bottomright", inset=0.01, legend=text.caption, cex=cex)
        par(mfg=c(nr-1, 2))
        plot(object$cvfit$cv.stats$mean$Support,
             object$cvfit$cv.stats$mean$MEFT,
             type='s', col=col, lty=lty, lwd=lwd,
             main="Maximum Event-Free Time", cex.main=cex,
             xlim=range(0, 1),
             ylim=range(0, object$cvfit$cv.stats$mean$MEFT, na.rm=TRUE),
             xlab=expression(paste("Box Support (", beta, ")", sep="")),
             ylab="MEFT", ...)
        if (add.caption)
          legend(x="bottomright", inset=0.01, legend=text.caption, cex=cex)
        par(mfg=c(nr-1, 3))
        plot(object$cvfit$cv.stats$mean$Support,
             object$cvfit$cv.stats$mean$MEFP,
             type='s', col=col, lty=lty, lwd=lwd,
             main="Minimum Event-Free Probability", cex.main=cex,
             xlim=range(0, 1),
             ylim=range(0, 1),
             xlab=expression(paste("Box Support (", beta, ")", sep="")),
             ylab="MEFP", ...)
        if (add.caption)
          legend(x="bottomright", inset=0.01, legend=text.caption, cex=cex)
        par(mfg=c(nr, 1))
        plot(object$cvfit$cv.stats$mean$Support,
             object$cvfit$cv.stats$mean$LHR,
             type='s', col=col, lty=lty, lwd=lwd,
             main="Log-Hazard Ratio", cex.main=cex,
             xlim=range(0, 1),
             ylim=range$lhr,
             xlab=expression(paste("Box Support (", beta, ")", sep="")),
             ylab=expression(paste("LHR (", lambda,")", sep="")), ...)
        if (add.caption)
          legend(x="top", inset=0.01, legend=text.caption, cex=cex)
        par(mfg=c(nr, 2))
        plot(object$cvfit$cv.stats$mean$Support,
             object$cvfit$cv.stats$mean$LRT,
             type='s', col=col, lty=lty, lwd=lwd,
             main="Log-Rank Test", cex.main=cex,
             xlim=range(0, 1),
             ylim=range$lrt,
             xlab=expression(paste("Box Support (", beta, ")", sep="")),
             ylab=expression(paste("LRT (", chi^2 ,")", sep="")), ...)
        if (add.caption)
          legend(x="top", inset=0.01, legend=text.caption, cex=cex)
        par(mfg=c(nr, 3))
        plot(object$cvfit$cv.stats$mean$Support,
             object$cvfit$cv.stats$mean$CER,
             type='s', col=col, lty=lty, lwd=lwd,
             main="Concordance Error Rate", cex.main=cex,
             xlim=range(0, 1),
             ylim=range$cer,
             xlab=expression(paste("Box Support (", beta, ")", sep="")),
             ylab=expression(paste("CER (", theta,")", sep="")), ...)
        if (add.caption)
          legend(x="top", inset=0.01, legend=text.caption, cex=cex)
        
      } else if (peelcriterion == "grp") {
        par(mfg=c(nr, 1))
        plot(object$cvfit$cv.stats$mean$Support,
             object$cvfit$cv.stats$mean$Size,
             type='s', col=col, lty=lty, lwd=lwd,
             main="Box Sample Size", cex.main=cex,
             xlim=range(0, 1),
             ylim=range(0, object$cvfit$cv.stats$mean$Size),
             xlab=expression(paste("Box Support (", beta, ")", sep="")),
             ylab=expression(paste("Size (", italic(n), ")", sep="")), ...)
        if (add.caption)
          legend(x="bottomright", inset=0.01, legend=text.caption, cex=cex)
        par(mfg=c(nr, 2))
        plot(object$cvfit$cv.stats$mean$Support,
             object$cvfit$cv.stats$mean$GLHR,
             type='s', col=col, lty=lty, lwd=lwd,
             main="Group Log-Hazard Ratio", cex.main=cex,
             xlim=range(0, 1),
             ylim=range$glhr,
             xlab=expression(paste("Box Support (", beta, ")", sep="")),
             ylab=expression(paste("GLHR (", lambda,")", sep="")), ...)
        if (add.caption)
          legend(x="top", inset=0.01, legend=text.caption, cex=cex)
        par(mfg=c(nr, 3))
        plot(object$cvfit$cv.stats$mean$Support,
             object$cvfit$cv.stats$mean$GCER,
             type='s', col=col, lty=lty, lwd=lwd,
             main="Group Concordance Error Rate", cex.main=cex,
             xlim=range(0, 1),
             ylim=range$gcer,
             xlab=expression(paste("Box Support (", beta, ")", sep="")),
             ylab=expression(paste("GCER (", theta,")", sep="")), ...)
        if (add.caption)
          legend(x="top", inset=0.01, legend=text.caption, cex=cex)
      } else {
        stop("Invalid peeling criterion. Exiting ...\n\n")
      }
      
      if (!is.null(main)) {
        mtext(text=main, cex=1, side=3, outer=TRUE)
      }
    }
    
    if (is.null(device)) {
      cat("Device: ",  dev.cur(), "\n")
      trajplot(object=object,
               main=main,
               toplot=toplot,
               range=range,
               col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
               col=col, lty=lty, lwd=lwd, cex=cex, 
               add.caption=add.caption, text.caption=text.caption,
               nr=nr, nc=nc)
    } else if (device == "PS") {
      path <- normalizePath(path=path, winslash="\\", mustWork=FALSE)
      file <- paste(file, ".ps", sep="")
      cat("\nOUTPUT: \n")
      cat("Filename : ", file, "\n")
      cat("Directory: ", path, "\n")
      postscript(file=file.path(path, file, fsep=.Platform$file.sep), width=width, height=height, onefile=TRUE, horizontal=horizontal)
      cat("Device: ",  dev.cur(), "\n")
      trajplot(object=object,
               main=main,
               toplot=toplot,
               range=range,
               col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
               col=col, lty=lty, lwd=lwd, cex=cex, 
               add.caption=add.caption, text.caption=text.caption,
               nr=nr, nc=nc)
      dev.off()
    } else if (device == "PDF") {
      path <- normalizePath(path=path, winslash="\\", mustWork=FALSE)
      file <- paste(file, ".pdf", sep="")
      cat("\nOUTPUT: \n")
      cat("Filename : ", file, "\n")
      cat("Directory: ", path, "\n")
      pdf(file=file.path(path, file, fsep=.Platform$file.sep), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
      cat("Device: ",  dev.cur(), "\n")
      trajplot(object=object,
               main=main,
               toplot=toplot,
               range=range,
               col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
               col=col, lty=lty, lwd=lwd, cex=cex, 
               add.caption=add.caption, text.caption=text.caption,
               nr=nr, nc=nc)
      dev.off()
    } else {
      cat("Currently allowed display devices are \"PS\" (Postscript) or \"PDF\" (Portable Document Format). Exiting ... \n\n")
    }
    
  } else {
    
    cat("Either the covariate screening or the Survival Bump Hunting modeling failed for this dataset.\n
         So, there is nothing to plot here.\n")
    
  }
  
  invisible()
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    plot_trace (object,
#                                main="Covariate Trace Plots",
#                                xlab=expression(paste("Box Support (", beta, ")", sep="")), 
#                                ylab="Covariate Range (centered)",
#                                toplot=object$cvfit$cv.used,
#                                center=TRUE, 
#                                scale=FALSE,
#                                col.cov, 
#                                lty.cov, 
#                                lwd.cov,
#                                col=1, 
#                                lty=1, 
#                                lwd=0.5, 
#                                cex=0.5,
#                                add.caption=FALSE, 
#                                text.caption=NULL,
#                                device=NULL, 
#                                file="Covariate Trace Plots", 
#                                path=getwd(),
#                                horizontal=FALSE, 
#                                width=8.5, 
#                                height=8.5, ...)
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

plot_trace <- function(object,
                       main="Covariate Trace Plots",
                       xlab=expression(paste("Box Support (", beta, ")", sep="")),
                       ylab="Covariate Range (centered)",
                       toplot=object$cvfit$cv.used,
                       center=TRUE,
                       scale=FALSE,
                       col.cov,
                       lty.cov,
                       lwd.cov,
                       col=1,
                       lty=1,
                       lwd=0.5,
                       cex=0.5,
                       add.caption=FALSE,
                       text.caption=NULL,
                       device=NULL,
                       file="Covariate Trace Plots",
                       path=getwd(),
                       horizontal=FALSE,
                       width=8.5,
                       height=8.5, ...) {
  
  if (!inherits(object, 'sbh'))
    stop("Argument `object` must be an object of class 'sbh'. Exiting ... \n\n")
  
  if (object$success) {
    
    traceplot <- function(object,
                          main, xlab, ylab,
                          toplot,
                          center, scale,
                          col.cov, lty.cov, lwd.cov,
                          col, lty, lwd,
                          cex, add.caption, text.caption, ...) {
      p <- length(toplot)
      varnames <- colnames(object$X)
      maxlength <- max(sapply(X=varnames, FUN=function(x){nchar(x, type="chars", allowNA=TRUE)}))
      
      if (missing(col.cov)) {
        col.cov <- 2:(p+1)
      }
      if (missing(lty.cov)) {
        lty.cov <- rep(1,p)
      }
      if (missing(lwd.cov)) {
        lwd.cov <- rep(1,p)
      }
      if (!is.null(main)) {
        par(mfrow=c(2, 1), oma=c(0, 0, 2, 0), mar=c(2.5, 2+maxlength/2, 2.0, 0.0), mgp=c(1.5, 0.5, 0))
      } else {
        par(mfrow=c(2, 1), oma=c(0, 0, 0, 0), mar=c(2.5, 2+maxlength/2, 2.0, 0.0), mgp=c(1.5, 0.5, 0))
      }
      
      boxcut.scaled <- scale(x=object$cvfit$cv.rules$mean[,varnames[toplot],drop=FALSE], center=center, scale=scale)
      plot(x=object$cvfit$cv.stats$mean$Support,
           y=boxcut.scaled[,1], type='n',
           xlim=range(0,1),
           ylim=range(boxcut.scaled, na.rm=TRUE),
           main="Covariate Importance (average values)", cex.main=cex,
           xlab="",
           ylab="", ...)
      for (j in 1:p) {
        lines(x=object$cvfit$cv.stats$mean$Support,
              y=boxcut.scaled[,j],
              type='l', col=col.cov[j], lty=lty.cov[j], lwd=lwd.cov[j], ...)
      }
      legend(x="topleft", inset=0.01, legend=varnames[toplot], col=col.cov, lty=lty.cov, lwd=lwd.cov, cex=cex)
      if (center)
        abline(h=0, lty=2, col=1, lwd=0.3, xpd=FALSE)
      if (add.caption)
        legend(x="bottom", inset=0.01, legend=text.caption, cex=cex)
      mtext(text=xlab, cex=cex, side=1, line=1, outer=FALSE)
      mtext(text=ylab, cex=cex, side=2, line=2, outer=FALSE)
      
      ticknames <- paste(varnames[toplot], " -", sep="")
      pointtrace <- c(object$cvfit$cv.trace[2], object$cvfit$cv.trace[-1])
      matchtrace <- pmatch(x=pointtrace, table=toplot, duplicates.ok = TRUE)
      plot(x=object$cvfit$cv.stats$mean$Support,
           y=matchtrace,
           type='S', yaxt="n", col=col, lty=lty, lwd=lwd,
           xlim=range(0, 1),
           ylim=range(0, p),
           main="Covariate Usage (modal values)", cex.main=cex,
           xlab="",
           ylab="", ...)
      par(mgp=c(1.5, 0, 0))
      axis(side=2, at=1:p, labels=ticknames, tick=FALSE, las=1, line=NA, cex.axis=cex, outer=FALSE)
      if (add.caption)
        legend(x="bottom", inset=0.01, legend=text.caption, cex=cex)
      mtext(text=xlab, cex=cex, side=1, line=1, outer=FALSE)
      mtext(text="Covariates Used", cex=cex, side=2, line=1+maxlength/2, outer=FALSE)
      if (!is.null(main)) {
        mtext(text=main, cex=1, side=3, outer=TRUE)
      }
    }
    
    if (is.null(device)) {
      cat("Device: ",  dev.cur(), "\n")
      traceplot(object=object,
                main=main, xlab=xlab, ylab=ylab,
                toplot=toplot,
                center=center, scale=scale,
                col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
                col=col, lty=lty, lwd=lwd,
                cex=cex, add.caption=add.caption, text.caption=text.caption)
    } else if (device == "PS") {
      path <- normalizePath(path=path, winslash="\\", mustWork=FALSE)
      file <- paste(file, ".ps", sep="")
      cat("\nOUTPUT: \n")
      cat("Filename : ", file, "\n")
      cat("Directory: ", path, "\n")
      postscript(file=file.path(path, file, fsep=.Platform$file.sep), width=width, height=height, onefile=TRUE, horizontal=horizontal)
      cat("Device: ",  dev.cur(), "\n")
      traceplot(object=object,
                main=main, xlab=xlab, ylab=ylab,
                toplot=toplot,
                center=center, scale=scale,
                col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
                col=col, lty=lty, lwd=lwd,
                cex=cex, add.caption=add.caption, text.caption=text.caption)
      dev.off()
    } else if (device == "PDF") {
      path <- normalizePath(path=path, winslash="\\", mustWork=FALSE)
      file <- paste(file, ".pdf", sep="")
      cat("\nOUTPUT: \n")
      cat("Filename : ", file, "\n")
      cat("Directory: ", path, "\n")
      pdf(file=file.path(path, file, fsep=.Platform$file.sep), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
      cat("Device: ",  dev.cur(), "\n")
      traceplot(object=object,
                main=main, xlab=xlab, ylab=ylab,
                toplot=toplot,
                center=center, scale=scale,
                col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
                col=col, lty=lty, lwd=lwd,
                cex=cex, add.caption=add.caption, text.caption=text.caption)
      dev.off()
    } else {
      cat("Currently allowed display devices are \"PS\" (Postscript) or \"PDF\" (Portable Document Format). Exiting ... \n\n")
    }
    
  } else {
    
    cat("Either the covariate screening or the Survival Bump Hunting modeling failed for this dataset.\n
         So, there is nothing to plot here.\n")
    
  }
  
  invisible()
  
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    plot_km (object,
#                             main="Survival KM Plots",
#                             xlab="Time", 
#                             ylab="Probability",
#                             ci=TRUE,
#                             precision=1e-3, 
#                             mark=3,
#                             col=c(1,2),
#                             lty=1, 
#                             lwd=0.5, 
#                             cex=0.5,
#                             steps=1:object$cvfit$cv.nsteps,
#                             plot.type="bumps",
#                             bump.reference="in-bump",
#                             group.reference=levels(object$groups)[1],
#                             add.caption=TRUE,
#                             text.caption=c("out-bump","in-bump"), 
#                             nr=3, 
#                             nc=4,
#                             device=NULL, 
#                             file="Survival KM Plots", 
#                             path=getwd(),
#                             horizontal=TRUE, 
#                             width=11, 
#                             height=8.5, ...)
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

plot_km <- function(object,
                    main="Survival KM Plots",
                    xlab="Time",
                    ylab="Probability",
                    ci=TRUE,	
                    precision=1e-3,
                    mark=3,
                    col=c(1,2), 
                    lty=1,
                    lwd=0.5,
                    cex=0.5,
                    steps=1:object$cvfit$cv.nsteps,
                    plot.type="bumps",
                    bump.reference="in-bump",
                    group.reference=levels(object$groups)[1],
                    add.caption=TRUE,
                    text.caption=c("out-bump","in-bump"), 
                    nr=3,
                    nc=4,
                    device=NULL,
                    file="Survival KM Plots",
                    path=getwd(),
                    horizontal=TRUE,
                    width=11,
                    height=8.5, ...) {
  
  if (!inherits(object, 'sbh'))
    stop("Argument `object` must be an object of class 'sbh'. Exiting ... \n\n")
  
  if (object$success) {
    
    kmplot <- function(object,
                       main, xlab, ylab,
                       precision, mark,
                       col, lty, lwd, cex,
                       steps, plot.type, 
                       bump.reference, group.reference,
                       add.caption, text.caption, 
                       nr, nc, ...) {
      
      if (!is.null(main)) {
        par(mfrow=c(nr, nc), oma=c(0, 0, 3, 0), mar=c(2.5, 2.5, 1.5, 1.5), mgp=c(1.5, 0.5, 0))
      } else {
        par(mfrow=c(nr, nc), oma=c(0, 0, 0, 0), mar=c(2.5, 2.5, 0.0, 1.5), mgp=c(1.5, 0.5, 0))
      }
      
      peelcriterion <- object$cvarg$peelcriterion
      y <- object$y
      delta <- object$delta
      L <- length(steps)
      for (l in 1:L) {
        i <- steps[l]
        if (bump.reference == "in-bump") {
          box1 <- 1 * object$cvfit$cv.boxind[i,]
        } else if (bump.reference == "out-bump") {
          box1 <- 1 * !object$cvfit$cv.boxind[i,]
        } else {
          stop("Invalid bump reference. Exiting ...\n\n")
        }
        wb <- which(box1 == 1)
        if ((peelcriterion == "lhr") || (peelcriterion == "lrt") || (peelcriterion == "cer")) {
          ng <- length(unique(box1[!is.na(box1)]))
          if (ng == 1) {
            box1 <- 1*box1
          } else {
            box1 <- 2 - 1*box1
          }
          if ((sum(box1, na.rm=TRUE) != length(box1[!is.na(box1)])) && (sum(box1, na.rm=TRUE) != 0)) {
            surv <- survival::survfit(survival::Surv(time=y, event=delta) ~ 1 + box1)
            if (ci) {
              plot(surv, xlim=range(0,y), ylim=range(0,1), main="", conf.int=TRUE, mark.time=FALSE, mark=NA, lty=c(2,2), lwd=c(lwd,lwd), col=c(col[2],col[1]), cex=cex, xlab=xlab, ylab=ylab, ...)
              par(new=TRUE)
            }
            plot(surv, xlim=range(0,y), ylim=range(0,1), main="", conf.int=FALSE, mark.time=TRUE, mark=mark, lty=c(lty,lty), lwd=c(lwd,lwd), col=c(col[2],col[1]), cex=cex, xlab=xlab, ylab=ylab, ...)
            if (add.caption) {
              legend(x="topright", inset=0.01, legend=rev(text.caption), lty=c(lty,lty), lwd=c(lwd,lwd), col=c(col[2],col[1]), cex=0.9*cex)
            }
            if (object$pv) {
              if (object$cvfit$cv.pval$pval[i] <= precision) {
                legend(x="bottom", inset=0.11, col="black", cex=0.9*cex, bty="n",
                       legend=bquote(italic(p) <= .(precision)))
              } else {
                legend(x="bottom", inset=0.11, col="black", cex=0.9*cex, bty="n",
                       legend=bquote(italic(p) == .(format(x=object$cvfit$cv.pval$pval[i], scientific=FALSE, digits=4, nsmall=4))))
              }
            }
            legend(x="bottom", inset=0.01, col="black", cex=0.9*cex, bty="n",
                   legend=substitute(group("", list(paste(italic(LHR) == x, sep="")), ""), list(x=format(x=object$cvfit$cv.stats$mean$LHR[i], digits=3, nsmall=3))))
            legend(x="bottom", inset=0.06, col="black", cex=0.9*cex, bty="n",
                   legend=substitute(group("", list(paste(italic(LRT) == x, sep="")), ""), list(x=format(x=object$cvfit$cv.stats$mean$LRT[i], digits=3, nsmall=3))))
            legend(x="bottom", inset=0.16, legend=paste("Step ", i, sep=""), col="black", cex=0.9*cex, bty="n")
            if (!is.null(main)) {
              mtext(text=main, cex=1, side=3, outer=TRUE)
            }
          }
        } else if (peelcriterion == "grp") {
          groups <- object$groups
          grp1 <- 1*(groups == group.reference)
          wg <- which(grp1 == 1)
          if (plot.type == "groups") {
            if ((sum(grp1[wb], na.rm=TRUE) != length(grp1[wb][!is.na(grp1[wb])])) && 
                (sum(grp1[-wb], na.rm=TRUE) != length(grp1[-wb][!is.na(grp1[-wb])])) && 
                (sum(grp1[wb], na.rm=TRUE) != 0) &&
                (sum(grp1[-wb], na.rm=TRUE) != 0)) {
              surv <- survival::survfit(survival::Surv(time=y[wb], event=delta[wb]) ~ 1 + grp1[wb])
              if (ci) {
                plot(surv, xlim=range(0,y), ylim=range(0,1), main="", conf.int=TRUE, mark.time=FALSE, mark=NA, lty=c(2,2), lwd=c(lwd,lwd), col=c(col[2],col[1]), cex=cex, xlab=xlab, ylab=ylab, ...)
                par(new=TRUE)
              }
              plot(surv, xlim=range(0,y), ylim=range(0,1), main="", conf.int=FALSE, mark.time=TRUE, mark=mark, lty=c(lty,lty), lwd=c(lwd,lwd), col=c(col[2],col[1]), cex=cex, xlab=xlab, ylab=ylab, ...)
              if (add.caption) {
                legend(x="topright", inset=0.01, legend=rev(text.caption), lty=c(lty,lty), lwd=c(lwd,lwd), col=c(col[2],col[1]), cex=0.9*cex)
              }
              legend(x="bottom", inset=0.01, legend=paste("Step ", i, sep=""), col="black", cex=0.9*cex, bty="n")
              if (!is.null(main)) {
                mtext(text=main, cex=1, side=3, outer=TRUE)
              }
            }
          } else if (plot.type == "bumps") {
            if ((sum(box1[wg], na.rm=TRUE) != length(box1[wg][!is.na(box1[wg])])) &&
                (sum(box1[-wg], na.rm=TRUE) != length(box1[-wg][!is.na(box1[-wg])])) &&
                (sum(box1[wg], na.rm=TRUE) != 0) &&
                (sum(box1[-wg], na.rm=TRUE) != 0)) {
              surv <- survival::survfit(formula=survival::Surv(time=y[wg], event=delta[wg]) ~ 1 + box1[wg], na.action="na.omit")
              if (ci) {
                plot(surv, xlim=range(0,y), ylim=range(0,1), main="", conf.int=TRUE, mark.time=FALSE, mark=NA, lty=c(2,2), lwd=c(lwd,lwd), col=c(col[2],col[1]), cex=cex, xlab=xlab, ylab=ylab, ...)
                par(new=TRUE)
              }
              plot(surv, xlim=range(0,y), ylim=range(0,1), main="", conf.int=FALSE, mark.time=TRUE, mark=mark, lty=c(lty,lty), lwd=c(lwd,lwd), col=c(col[2],col[1]), cex=cex, xlab=xlab, ylab=ylab, ...)
              if (add.caption) {
                legend(x="topright", inset=0.01, legend=rev(text.caption), lty=c(lty,lty), lwd=c(lwd,lwd), col=c(col[2],col[1]), cex=0.9*cex)
              }
              legend(x="bottom", inset=0.01, legend=paste("Step ", i, sep=""), col="black", cex=0.9*cex, bty="n")
              if (!is.null(main)) {
                mtext(text=main, cex=1, side=3, outer=TRUE)
              }
            }
          } else {
            stop("Invalid plot type. Exiting ...\n\n")
          }
        } else {
          stop("Invalid peeling criterion. Exiting ...\n\n")
        }
      }
    }
    
    if (is.null(device)) {
      cat("Device: ",  dev.cur(), "\n")
      kmplot(object=object,
             main=main, xlab=xlab, ylab=ylab,
             precision=precision, mark=mark,
             col=col, lty=lty, lwd=lwd, cex=cex,
             steps=steps, plot.type=plot.type, 
             bump.reference=bump.reference, group.reference=group.reference,
             add.caption=add.caption, text.caption=text.caption,
             nr=nr, nc=nc)
    } else if (device == "PS") {
      path <- normalizePath(path=path, winslash="\\", mustWork=FALSE)
      file <- paste(file, ".ps", sep="")
      cat("\nOUTPUT: \n")
      cat("Filename : ", file, "\n")
      cat("Directory: ", path, "\n")
      postscript(file=file.path(path, file, fsep=.Platform$file.sep), width=width, height=height, onefile=TRUE, horizontal=horizontal)
      cat("Device: ",  dev.cur(), "\n")
      kmplot(object=object,
             main=main, xlab=xlab, ylab=ylab,
             precision=precision, mark=mark,
             col=col, lty=lty, lwd=lwd, cex=cex,
             steps=steps, plot.type=plot.type, 
             bump.reference=bump.reference, group.reference=group.reference,
             add.caption=add.caption, text.caption=text.caption,
             nr=nr, nc=nc)
      dev.off()
    } else if (device == "PDF") {
      path <- normalizePath(path=path, winslash="\\", mustWork=FALSE)
      file <- paste(file, ".pdf", sep="")
      cat("\nOUTPUT: \n")
      cat("Filename : ", file, "\n")
      cat("Directory: ", path, "\n")
      pdf(file=file.path(path, file, fsep=.Platform$file.sep), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
      cat("Device: ",  dev.cur(), "\n")
      kmplot(object=object,
             main=main, xlab=xlab, ylab=ylab,
             precision=precision, mark=mark,
             col=col, lty=lty, lwd=lwd, cex=cex,
             steps=steps, plot.type=plot.type, 
             bump.reference=bump.reference, group.reference=group.reference,
             add.caption=add.caption, text.caption=text.caption,
             nr=nr, nc=nc)
      dev.off()
    } else {
      cat("Currently allowed display devices are \"PS\" (Postscript) or \"PDF\" (Portable Document Format). Exiting ... \n\n")
    }
    
  } else {
    
    cat("Either the covariate screening or the Survival Bump Hunting modeling failed for this dataset.\n
         So, there is nothing to plot here.\n")
    
  }
  
  invisible()
  
}
#===============================================================================================================================#
