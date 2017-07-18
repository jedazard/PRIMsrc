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
#                       delta,
#                       B=30,
#                       K=5,
#                       A=1000,
#                       vs=TRUE,
#                       vstype="ppl",
#                       vsarg="alpha=1,
#                              nalpha=1,
#                              nlambda=100,
#                              vscons=0.5",
#                       cv=TRUE,
#                       cvtype="combined",
#                       cvarg="alpha=0.01,
#                              beta=0.05,
#                              minn=5,
#                              L=NULL,
#                              peelcriterion=\"lrt\"",
#                              cvcriterion=\"cer\"",
#                       pv=FALSE,
#                       decimals=2,
#                       onese=FALSE,
#                       probval=NULL,
#                       timeval=NULL,
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
                delta,
                B=30,
                K=5,
                A=1000,
                vs=TRUE,
                vstype="ppl",
                vsarg="alpha=1,
               nalpha=1,
               nlambda=100,
               vscons=0.5",
                cv=TRUE,
                cvtype="combined",
                cvarg="alpha=0.01,
               beta=0.05,
               minn=5,
               L=NULL,
               peelcriterion=\"lrt\",
               cvcriterion=\"cer\"",
                pv=FALSE,
                decimals=2,
                onese=FALSE,
                probval=NULL,
                timeval=NULL,
                parallel.vs=FALSE,
                parallel.rep=FALSE,
                parallel.pv=FALSE,
                conf=NULL,
                verbose=TRUE,
                seed=NULL) {

   # Parsing and evaluating PRSP parameters of SBH
   alpha <- NULL
   beta <- NULL
   minn <- NULL
   L <- NULL
   peelcriterion <- NULL
   cvcriterion <- NULL
   eval(parse( text=unlist(strsplit(x=cvarg, split=",")) ))

   # Checking matching of parameters
   if (vs) {
      vstype <- match.arg(arg=vstype, choices=c("pcqr", "ppl", "spca", "prsp"), several.ok=FALSE)
   } else {
      vstype <- NA
   }
   if (cv) {
      cvtype <- match.arg(arg=cvtype, choices=c("combined", "averaged"), several.ok=FALSE)
   } else {
      cvtype <- NA
   }

   # Checking inputs
   B <- as.integer(B)
   K <- as.integer(K)
   A <- as.integer(A)
   if (B <= 0) {
      stop("\n'B' must be a positive integer. Exiting ... \n\n")
   }
   if (K <= 0) {
      stop("\n'K' must be a positive integer. Exiting ... \n\n")
   }
   if (A <= 0) {
      stop("\n'A' must be a positive integer. Exiting ... \n\n")
   }
   if (missing(X)) {
      stop("\nNo dataset provided! Exiting ... \n\n")
   }
   if (!(is.matrix(X))) {
      X <- as.matrix(X)
   }
   mode(X) <- "numeric"
   n <- nrow(X)
   p <- ncol(X)
   if (is.null(rownames(X))) {
      rownames(X) <- 1:n
   }
   if (is.null(colnames(X))) {
      colnames(X) <- paste("X", 1:p, sep="")
   }
   mini <- max(minn, ceiling(beta*n))
   if (n < mini) {
      stop("The number of data points must be greater than the threshold of ", mini, " points. Exiting ... \n\n")
   }
   digits <- getOption("digits")
   y[y <= 0] <- 10^(-digits)

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
   cat("PRSP cross-validation criterion: ", toupper(x=cvcriterion), "\n")
   cat("PRSP Peeling criterion: ", toupper(x=peelcriterion), "\n")
   cat("PRSP Peeling percentile: ", alpha*100, "%\n")
   cat("PRSP Minimal box support: ", beta*100, "%\n")
   cat("PRSP Minimal box sample size: ", minn, "\n")
   cat("PRSP Maximum peeling length: ", ifelse(is.null(L), yes="AUTOMATIC", no=L), "\n")
   cat("Computation of p-values: ", pv, "\n")
   cat("Decision rule: ", ifelse(onese, yes="1SE", no="EXTREMUM"), "\n")
   if (vs) {
      cat("Parallelization of computation of variable screening:", parallel.vs, "\n")
   }
   cat("Parallelization of computation of PRSP:", parallel.rep, "\n")
   if (pv) {
      cat("Parallelization of computation of p-values:", parallel.pv, "\n")
   }

   # Setting the cluster up
   if ((parallel.rep) || (parallel.vs) || (parallel.pv)) {
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
   }

   # Screening of informative covariates
   cv.presel.obj <- cv.presel(X=X,
                              y=y,
                              delta=delta,
                              B=B,
                              K=K,
                              vs=vs,
                              vstype=vstype,
                              vsarg=vsarg,
                              cv=cv,
                              cvtype=cvtype,
                              onese=onese,
                              parallel.vs=parallel.vs,
                              parallel.rep=parallel.rep,
                              clus.vs=cluster,
                              clus.rep=cluster,
                              verbose=verbose,
                              seed=seed)

   success <- cv.presel.obj$success
   seed <- cv.presel.obj$seed

   # Survival Bump Hunting model using the PRSP algorithm
   if (!success) {

      if (vs) cat("Did not find any informative covariates. Exiting ... \n\n", sep="")

      # List of CV profiles
      CV.profiles <- NULL
      # List of CV fitted values
      CV.fit <- NULL

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

      # Fitting the Survival Bump Hunting model using the PRSP algorithm
      cat("Fitting the Survival Bump Hunting model using the PRSP algorithm ... \n")
      CV.box.obj <- cv.box(X=X.sel,
                           y=y,
                           delta=delta,
                           B=B,
                           K=K,
                           cv=cv,
                           cvtype=cvtype,
                           cvarg=paste("alpha=", alpha,
                                       ",beta=", beta,
                                       ",minn=", minn,
                                       ",L=", L,
                                       ",peelcriterion=\"", peelcriterion, "\"", sep=""),
                           decimals=decimals,
                           probval=probval,
                           timeval=timeval,
                           varsign=CV.screened.sign,
                           initcutpts=initcutpts,
                           parallel.rep=parallel.rep,
                           clus.rep=cluster,
                           verbose=verbose,
                           seed=seed)

      success <- CV.box.obj$success
      seed <- CV.box.obj$seed

      if (!success) {

         cat("The PRSP algorithm could not find any bump in this dataset. Exiting ... \n\n", sep="")

         # List of CV profiles
         CV.profiles <- NULL
         # List of CV fitted values
         CV.fit <- NULL

      } else {

         CV.maxsteps <- CV.box.obj$cv.maxsteps
         CV.trace.list <- CV.box.obj$cv.trace
         CV.boxind <- CV.box.obj$cv.boxind
         CV.boxcut <- CV.box.obj$cv.boxcut
         CV.support <- CV.box.obj$cv.support
         CV.size <- CV.box.obj$cv.size
         CV.lhr <- CV.box.obj$cv.lhr
         CV.lrt <- CV.box.obj$cv.lrt
         CV.cer <- CV.box.obj$cv.cer
         CV.time.bar <- CV.box.obj$cv.time.bar
         CV.prob.bar <- CV.box.obj$cv.prob.bar
         CV.max.time.bar <- CV.box.obj$cv.max.time.bar
         CV.min.prob.bar <- CV.box.obj$cv.min.prob.bar

         # Maximum peeling length from all replicates
         CV.maxsteps <- ceiling(mean(CV.maxsteps))

         # Box membership indicator vector of all observations at each step over the replicates
         CV.boxind <- lapply.array(X=CV.boxind,
                                   rowtrunc=CV.maxsteps,
                                   MARGIN=c(1,2),
                                   FUN=function(x) {
                                      return(mean(x, na.rm=TRUE) >= 0.5)
                                   })
         dimnames(CV.boxind) <- list(paste("step", 0:(CV.maxsteps-1), sep=""), rownames(X.sel))

         # Adjusted maximum peeling length, thresholded by minimal box support, from all replicates
         CV.maxsteps <- max(which(apply(CV.boxind, 1, function(x) {length(which(x))/n >= max(minn/n, beta)})))

         # Adjusted cross-validated profiles of peeling steps and optimal peeling lengths from all replicates
         if (!cv) {
            CV.stepprofiles <- NULL
            CV.stepprofiles.mean <- NULL
            CV.stepprofiles.se <- NULL
            CV.nsteps.opt <- CV.maxsteps
            CV.nsteps.1se <- CV.maxsteps
         } else if ((cvtype == "averaged") || (cvtype == "combined")) {
            cat("Generating cross-validated profiles of peeling steps and optimal peeling lengths from all replicates ...\n")
            CV.lhr.mat <- list2mat(list=CV.lhr, fill=NA, coltrunc=CV.maxsteps)
            CV.lrt.mat <- list2mat(list=CV.lrt, fill=NA, coltrunc=CV.maxsteps)
            CV.cer.mat <- list2mat(list=CV.cer, fill=NA, coltrunc=CV.maxsteps)
            colnames(CV.lhr.mat) <- paste("step", 0:(CV.maxsteps-1), sep="")
            colnames(CV.lrt.mat) <- paste("step", 0:(CV.maxsteps-1), sep="")
            colnames(CV.cer.mat) <- paste("step", 0:(CV.maxsteps-1), sep="")
            if (cvcriterion=="lhr") {
               CV.stepprofiles <- CV.lhr.mat
               CV.stepprofiles.mean <- apply(CV.lhr.mat, 2, mean, na.rm=TRUE)
               CV.stepprofiles.se <- apply(CV.lhr.mat, 2, sd, na.rm=TRUE)
               if (all(is.na(CV.stepprofiles.mean)) || is.empty(CV.stepprofiles.mean)) {
                  CV.nsteps.opt <- NA
                  CV.nsteps.1se <- NA
               } else {
                  CV.nsteps.opt <- which.max(CV.stepprofiles.mean)
                  w <- CV.stepprofiles.mean >= CV.stepprofiles.mean[CV.nsteps.opt]-CV.stepprofiles.se[CV.nsteps.opt]
                  if (all(is.na(w)) || is.empty(w)) {
                     CV.nsteps.1se <- NA
                  } else {
                     CV.nsteps.1se <- min(which(w))
                  }
               }
            } else if (cvcriterion=="lrt") {
               CV.stepprofiles <- CV.lrt.mat
               CV.stepprofiles.mean <- apply(CV.lrt.mat, 2, mean, na.rm=TRUE)
               CV.stepprofiles.se <- apply(CV.lrt.mat, 2, sd, na.rm=TRUE)
               if (all(is.na(CV.stepprofiles.mean)) || is.empty(CV.stepprofiles.mean)) {
                  CV.nsteps.opt <- NA
                  CV.nsteps.1se <- NA
               } else {
                  CV.nsteps.opt <- which.max(CV.stepprofiles.mean)
                  w <- CV.stepprofiles.mean >= CV.stepprofiles.mean[CV.nsteps.opt]-CV.stepprofiles.se[CV.nsteps.opt]
                  if (all(is.na(w)) || is.empty(w)) {
                     CV.nsteps.1se <- NA
                  } else {
                     CV.nsteps.1se <- min(which(w))
                  }
               }
            } else if (cvcriterion=="cer") {
               CV.stepprofiles <- CV.cer.mat
               CV.stepprofiles.mean <- apply(CV.cer.mat, 2, mean, na.rm=TRUE)
               CV.stepprofiles.se <- apply(CV.cer.mat, 2, sd, na.rm=TRUE)
               if (all(is.na(CV.stepprofiles.mean)) || is.empty(CV.stepprofiles.mean)) {
                  CV.nsteps.opt <- NA
                  CV.nsteps.1se <- NA
               } else {
                  CV.nsteps.opt <- which.min(CV.stepprofiles.mean)
                  w <- CV.stepprofiles.mean <= CV.stepprofiles.mean[CV.nsteps.opt]+CV.stepprofiles.se[CV.nsteps.opt]
                  if (all(is.na(w)) || is.empty(w)) {
                     CV.nsteps.1se <- NA
                  } else {
                     CV.nsteps.1se <- min(which(w))
                  }
               }
            } else {
               stop("Invalid CV criterion option. Exiting ... \n\n")
            }
         } else {
            stop("Invalid CV type option. Exiting ... \n\n")
         }
         if (onese) {
            CV.nsteps <- CV.nsteps.1se
         } else {
            CV.nsteps <- CV.nsteps.opt
         }

         # Adjusted box membership indicator of all observations at each step using the modal or majority vote value over the replicates
         cat("Generating box memberships ...\n")
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
         CV.trace <- sapply(X=CV.trace[-1], FUN=function(x) {x[1]})
         m <- pmatch(x=names(CV.screened)[CV.trace], table=colnames(X), nomatch=NA, duplicates.ok=TRUE)
         CV.trace <- c(0, m)
         names(CV.trace) <- paste("step", 0:(CV.nsteps-1), sep="")
         out <- c(1, which(is.na(CV.trace)))

         if (length(out) == length(CV.trace)) {

            success <- FALSE
            cat("The PRSP algorithm could not find any bump variables in this dataset. Exiting ... \n\n", sep="")

            # List of CV profiles
            CV.profiles <- NULL
            # List of CV fitted values
            CV.fit <- NULL

         } else {

            success <- TRUE
            cat("Successfully completed PRSP algorithm. \n", sep="")

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
            CV.boxcut.mu <- round(lapply.array(X=CV.boxcut, rowtrunc=CV.nsteps, FUN=function(x){mean(x, na.rm=TRUE)}, MARGIN=1:2), digits=decimals)
            CV.boxcut.sd <- round(lapply.array(X=CV.boxcut, rowtrunc=CV.nsteps, FUN=function(x){sd(x, na.rm=TRUE)}, MARGIN=1:2), digits=decimals)
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
            CV.support.sd <- round(lapply.mat(X=CV.support, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.size.mu <- round(CV.boxind.size, digits=0)
            CV.size.sd <- round(lapply.mat(X=CV.size, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.lhr.mu <- round(lapply.mat(X=CV.lhr, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.lhr.sd <- round(lapply.mat(X=CV.lhr, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.lrt.mu <- round(lapply.mat(X=CV.lrt, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.lrt.sd <- round(lapply.mat(X=CV.lrt, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.cer.mu <- round(lapply.mat(X=CV.cer, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.cer.sd <- round(lapply.mat(X=CV.cer, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.time.bar.mu <- round(lapply.mat(X=CV.time.bar, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.time.bar.sd <- round(lapply.mat(X=CV.time.bar, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.prob.bar.mu <- round(lapply.mat(X=CV.prob.bar, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.prob.bar.sd <- round(lapply.mat(X=CV.prob.bar, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.max.time.bar.mu <- round(lapply.mat(X=CV.max.time.bar, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.max.time.bar.sd <- round(lapply.mat(X=CV.max.time.bar, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.min.prob.bar.mu <- round(lapply.mat(X=CV.min.prob.bar, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.min.prob.bar.sd <- round(lapply.mat(X=CV.min.prob.bar, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
            CV.stats.mu <- data.frame("Support"=CV.support.mu,
                                      "Size"=CV.size.mu,
                                      "LHR"=CV.lhr.mu,
                                      "LRT"=CV.lrt.mu,
                                      "CER"=CV.cer.mu,
                                      "EFT"=CV.time.bar.mu,
                                      "EFP"=CV.prob.bar.mu,
                                      "MEFT"=CV.max.time.bar.mu,
                                      "MEFP"=CV.min.prob.bar.mu)
            rownames(CV.stats.mu) <- paste("step", 0:(CV.nsteps-1), sep="")
            colnames(CV.stats.mu) <- c("Support", "Size", "LHR", "LRT", "CER", "EFT", "EFP", "MEFT", "MEFP")
            CV.stats.sd <- data.frame("Support"=CV.support.sd,
                                      "Size"=CV.size.sd,
                                      "LHR"=CV.lhr.sd,
                                      "LRT"=CV.lrt.sd,
                                      "CER"=CV.cer.sd,
                                      "EFT"=CV.time.bar.sd,
                                      "EFP"=CV.prob.bar.sd,
                                      "MEFT"=CV.max.time.bar.sd,
                                      "MEFP"=CV.min.prob.bar.sd)
            rownames(CV.stats.sd) <- paste("step", 0:(CV.nsteps-1), sep="")
            colnames(CV.stats.sd) <- c("Support", "Size", "LHR", "LRT", "CER", "EFT", "EFP", "MEFT", "MEFP")
            CV.stats <- list("mean"=CV.stats.mu, "sd"=CV.stats.sd)

            # Computation of p-values at each step
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
                                           ",minn=", minn,
                                           ",L=", CV.nsteps,
                                           ",peelcriterion=\"", peelcriterion, "\"", sep=""),
                               obs.chisq=CV.stats$mean$LRT,
                               parallel.pv=parallel.pv,
                               clus.pv=cluster,
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

            # Return object 'CV.fit'
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

            cat("Finished!\n")
         }
      }
   }

   # Stopping the cluster and cleaning all MPI states
   if ((parallel.rep) || (parallel.vs) || (parallel.pv)) {
      parallel::stopCluster(cl=cluster)
   }

   # Returning the final 'sbh' object
   return(structure(list("X"=X,
                         "y"=y,
                         "delta"=delta,
                         "B"=B,
                         "K"=K,
                         "A"=A,
                         "vs"=vs,
                         "vstype"=vstype,
                         "vsarg"=vsarg,
                         "cv"=cv,
                         "cvtype"=cvtype,
                         "cvarg"=cvarg,
                         "pv"=pv,
                         "decimals"=decimals,
                         "onese"=onese,
                         "probval"=probval,
                         "timeval"=timeval,
                         "cvprofiles"=CV.profiles,
                         "cvfit"=CV.fit,
                         "success"=success,
                         "seed"=seed),
                    class="sbh"))
}
#===============================================================================================================================#




#===============================================================================================================================#
# 2. END-USER FUNCTION FOR NEWS
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
# 3. END-USER S3-GENERIC FUNCTIONS FOR SUMMARY, PRINT, PLOT AND PREDICTION
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
      stop("Primary argument `object` must be an object of class 'sbh'. Exiting ... \n\n")

   alpha <- NULL
   beta <- NULL
   minn <- NULL
   L <- NULL
   peelcriterion <- NULL
   cvcriterion <- NULL
   eval(parse( text=unlist(strsplit(x=object$cvarg, split=",")) ))

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
   cat("\t Variable screening technique: ", toupper(x=object$vstype), "\n\n")
   cat("CROSS-VALIDATION:\n")
   cat("\t Cross-validation: ", object$cv, "\n")
   cat("\t Cross-validation technique: ", toupper(x=object$cvtype), "\n\n")
   cat("PRSP PARAMETERS:\n")
   cat("\t Cross-validation criterion: ", toupper(x=cvcriterion), "\n")
   cat("\t Peeling criterion: ", toupper(x=peelcriterion), "\n")
   cat("\t Peeling percentile: ", alpha*100, "%\n")
   cat("\t Minimal box support: ", beta*100, "%\n")
   cat("\t Minimal box sample size: ", minn, "\n")
   cat("\t Maximum peeling length: ", ifelse(is.null(L), yes="AUTOMATIC", no=L), "\n\n")
   cat("REPORTING:\n")
   cat("\t Decision rule: ", ifelse(object$onese, yes="1SE", no="EXTREMUM"), "\n")
   cat("\t Number of decimals: ", object$decimals, "\n")
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
      stop("Primary argument `x` must be an object of class 'sbh'. Exiting ... \n\n")

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
#                         main=NULL,
#                         proj=c(1,2), splom=TRUE, boxes=FALSE,
#                         steps=x$cvfit$cv.nsteps,
#                         pch=16, cex=0.5, col=, col=2:(length(steps)+1),
#                         col.box=2:(length(steps)+1), lty.box=rep(2,length(steps)), lwd.box=rep(1,length(steps)),
#                         add.legend=TRUE,
#                         device=NULL, file="Scatter Plot", path=getwd(),
#                         horizontal=FALSE, width=5, height=5, ...)
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
                     main=NULL,
                     proj=c(1,2), splom=TRUE, boxes=FALSE,
                     steps=x$cvfit$cv.nsteps,
                     pch=16, cex=0.5, col=2:(length(steps)+1),
                     col.box=2:(length(steps)+1), lty.box=rep(2,length(steps)), lwd.box=rep(1,length(steps)),
                     add.legend=TRUE,
                     device=NULL, file="Scatter Plot", path=getwd(),
                     horizontal=FALSE, width=5, height=5, ...) {

   if (!inherits(x, 'sbh'))
      stop("Primary argument `x` must be an object of class 'sbh'. Exiting ... \n\n")

   if (x$success) {

      scatterplot <- function(object,
                              main,
                              proj, splom, boxes,
                              steps,
                              add.legend, pch, cex, col,
                              col.box, lty.box, lwd.box, ...) {

         if (!is.null(main)) {
            par(mfrow=c(1, 1), oma=c(0, 0, 3, 0), mar=c(2.5, 2.5, 4.0, 1.5), mgp=c(1.5, 0.5, 0))
         } else {
            par(mfrow=c(1, 1), oma=c(0, 0, 0, 0), mar=c(2.5, 2.5, 4.0, 1.5), mgp=c(1.5, 0.5, 0))
         }

         toplot <- object$cvfit$cv.used[proj]
         varnames <- colnames(object$X)
         x <- object$X[,varnames[toplot],drop=FALSE]
         x.names <- colnames(x)

         if (is.null(steps))
            steps <- object$cvfit$cv.nsteps

         L <- length(steps)
         plot(x=x, type="p", pch=pch, cex=cex, col=1, main=NULL, xlab=x.names[1], ylab=x.names[2], ...)
         if (splom) {
            for (i in 1:L) {
               w <- object$cvfit$cv.boxind[steps[i],]
               points(x=object$X[w,varnames[toplot],drop=FALSE], type="p", pch=pch, cex=cex, col=col[i], ...)
            }
         }
         if (boxes) {
            x.range <- apply(X=x, MARGIN=2, FUN=range)
            boxcut <- object$cvfit$cv.rules$mean[steps,varnames[toplot],drop=FALSE]
            varsign <- object$cvfit$cv.sign[varnames[toplot]]
            vertices <- vector(mode="list", length=L)
            for (i in 1:L) {
               vertices[[i]] <- matrix(data=NA, nrow=2, ncol=2, dimnames=list(c("LB","UB"), x.names))
               for (j in 1:2) {
                  vertices[[i]][1,j] <- ifelse(test=(varsign[j] > 0),
                                               yes=max(x.range[1,j], boxcut[i,j]),
                                               no=min(x.range[1,j], boxcut[i,j]))
                  vertices[[i]][2,j] <- ifelse(test=(varsign[j] < 0),
                                               yes=min(x.range[2,j], boxcut[i,j]),
                                               no=max(x.range[2,j], boxcut[i,j]))
               }
            }
            for (i in 1:L) {
               rect(vertices[[i]][1,1], vertices[[i]][1,2], vertices[[i]][2,1], vertices[[i]][2,2],
                    border=col.box[i], col=NA, lty=lty.box[i], lwd=lwd.box[i])
            }
         }
         if (!is.null(main)) {
            mtext(text=main, cex=1, side=3, outer=TRUE)
         }
         if (add.legend) {
            legend("topleft", xpd=TRUE, inset=0.01, legend=paste("Step: ", steps, sep=""), pch=pch, col=col, cex=cex)
         }
      }

      if (is.null(device)) {
         scatterplot(object=x,
                     main=main,
                     proj=proj, splom=splom, boxes=boxes, steps=steps,
                     add.legend=add.legend, pch=pch, cex=cex, col=col,
                     col.box=col.box, lty.box=lty.box, lwd.box=lwd.box)
      } else if (device == "PS") {
         path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
         file <- paste(file, ".ps", sep="")
         cat("\nOUTPUT: \n")
         cat("Filename : ", file, "\n")
         cat("Directory: ", path, "\n")
         postscript(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, horizontal=horizontal)
         scatterplot(object=x,
                     main=main,
                     proj=proj, splom=splom, boxes=boxes, steps=steps,
                     add.legend=add.legend, pch=pch, cex=cex, col=col,
                     col.box=col.box, lty.box=lty.box, lwd.box=lwd.box)
         dev.off()
      } else if (device == "PDF") {
         path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
         file <- paste(file, ".pdf", sep="")
         cat("\nOUTPUT: \n")
         cat("Filename : ", file, "\n")
         cat("Directory: ", path, "\n")
         pdf(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
         scatterplot(object=x,
                     main=main,
                     proj=proj, splom=splom, boxes=boxes, steps=steps,
                     add.legend=add.legend, pch=pch, cex=cex, col=col,
                     col.box=col.box, lty.box=lty.box, lwd.box=lwd.box)
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
#                           steps,
#                           na.action = na.omit, ...)
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
                         steps,
                         na.action = na.omit, ...) {

   if (object$success) {

      if (!inherits(object, 'sbh'))
         stop("Primary argument `object` must be an object of class 'sbh'. Exiting ... \n\n")

      x <- as.matrix(newdata)
      x.names <- colnames(x)
      x.range <- apply(X=x, MARGIN=2, FUN=range)
      n <- nrow(x)
      p <- ncol(x)

      toplot <- object$cvfit$cv.used
      varnames <- colnames(object$X)

      if (length(toplot) != p) {
         stop("Non-matching dimensionality of newdata to 'sbh' object of used covariates. Exiting ... \n\n")
      }

      if (missing(steps) || is.null(steps))
         steps <- object$cvfit$cv.nsteps

      L <- length(steps)
      boxcut <- object$cvfit$cv.rules$mean[steps,varnames[toplot],drop=FALSE]
      varsign <- object$cvfit$cv.sign[varnames[toplot]]

      pred.boxind <- matrix(NA, nrow=L, ncol=n, dimnames=list(paste("step ", steps, sep=""), rownames(x)))
      for (l in 1:L) {
         boxcutsign <- boxcut[l, ] * varsign
         x.cut <- t(t(x) * varsign)
         x.ind <- t(t(x.cut) >= boxcutsign)
         pred.boxind[l,] <- (rowMeans(x.ind) == 1)  # Set as TRUE which observations are inside the box boundaries for all axes directions
      }

      pred.vertices <- vector(mode="list", length=L)
      names(pred.vertices) <- paste("step ", steps, sep="")
      for (i in 1:L) {
         pred.vertices[[i]] <- matrix(data=NA, nrow=2, ncol=p, dimnames=list(c("LB","UB"), x.names))
         for (j in 1:p) {
            pred.vertices[[i]][1,j] <- ifelse(test=(varsign[j] > 0),
                                              yes=max(x.range[1,j], boxcut[i,j]),
                                              no=min(x.range[1,j], boxcut[i,j]))
            pred.vertices[[i]][2,j] <- ifelse(test=(varsign[j] < 0),
                                              yes=min(x.range[2,j], boxcut[i,j]),
                                              no=max(x.range[2,j], boxcut[i,j]))
         }
      }

   } else {

      cat("Either the covariate screening or the Survival Bump Hunting modeling failed for this dataset.\n
           So, there is nothing to predict here.\n")

   }

   return(list("boxind"=pred.boxind, "vertices"=pred.vertices))
}
#===============================================================================================================================#




#===============================================================================================================================#
# 4. END-USER PLOTTING FUNCTIONS FOR MODEL VALIDATION AND VISUALIZATION OF RESULTS
#===============================================================================================================================#

#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    plot_profile(object,
#                                 main=NULL,
#                                 xlim=NULL, ylim=NULL,
#                                 add.sd=TRUE, add.legend=TRUE, add.profiles=TRUE,
#                                 pch=20, col=1, lty=1, lwd=0.5, cex=0.5,
#                                 device=NULL, file="Profile Plots", path=getwd(),
#                                 horizontal=FALSE, width=8.5, height=5.0, ...) {
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
                         main=NULL,
                         xlim=NULL, ylim=NULL,
                         add.sd=TRUE, add.legend=TRUE, add.profiles=TRUE,
                         pch=20, col=1, lty=1, lwd=0.5, cex=0.5,
                         device=NULL, file="Profile Plots", path=getwd(),
                         horizontal=FALSE, width=8.5, height=5.0, ...) {

   if (!inherits(object, 'sbh'))
      stop("Primary argument `object` must be an object of class 'sbh'. Exiting ... \n\n")

   if (object$success) {

      if (object$cv) {

         profileplot <- function(object, main, xlim, ylim,
                                 add.sd, add.legend, add.profiles,
                                 pch, col, lty, lwd, cex, ...) {

            # Parsing and evaluating PRSP parameters of SBH to retrieve 'cvcriterion' for CV of peeling length
            alpha <- NULL
            beta <- NULL
            minn <- NULL
            L <- NULL
            peelcriterion <- NULL
            cvcriterion <- NULL
            eval(parse( text=unlist(strsplit(x=object$cvarg, split=",")) ))

            if (cvcriterion == "lhr") {
               txt <- "LHR"
            } else if (cvcriterion == "lrt") {
               txt <- "LRT"
            } else if (cvcriterion == "cer") {
               txt <- "CER"
            } else {
               stop("Invalid CV criterion. Exiting ... \n\n")
            }

            if (object$vs && object$vstype == "prsp" && ncol(object$cvprofiles$cv.varprofiles) > 1) {
               if (!is.null(main)) {
                  par(mfrow=c(2, 1), oma=c(0, 0, 3, 0), mar=c(4.0, 3.0, 4.0, 3.0), mgp=c(1.5, 0.5, 0))
               } else {
                  par(mfrow=c(2, 1), oma=c(0, 0, 0, 0), mar=c(4.0, 3.0, 4.0, 3.0), mgp=c(1.5, 0.5, 0))
               }
               varprofiles <- object$cvprofiles$cv.varprofiles
               varprofiles.mean <- object$cvprofiles$cv.varprofiles.mean
               varprofiles.se <- object$cvprofiles$cv.varprofiles.se
               if (object$onese) {
                  cv.varset <- object$cvprofiles$cv.varset.1se
               } else {
                  cv.varset <- object$cvprofiles$cv.varset.opt
               }

               # Parsing and evaluating PRSP parameters of SBH to retrieve 'S' for CV of variable selection
               alpha <- NULL
               beta <- NULL
               minn <- NULL
               L <- NULL
               S <- NULL
               peelcriterion <- NULL
               cvcriterion <- NULL
               vscons <- NULL
               eval(parse( text=unlist(strsplit(x=object$vsarg, split=",")) ))

               # Grid of cardinals of top-ranked variables subsets
               n <- nrow(object$X)
               p <- ncol(object$X)
               Lmax <- ceiling(log(1/n) / log(1 - (1/n)))
               if (L > Lmax) {
                  L <- pmin(Lmax, L)
               }
               L <- max(1,floor(L))
               if (is.null(S)) {
                  S <- max(1,floor(p/10))
                  size <- unique(ceiling(seq(from=max(1,floor(S/100)), to=S, length=min(S,floor(100*S/p)))))
               } else {
                  if (S > p) {
                     S <- pmin(p, S)
                  }
                  size <- max(1,floor(S))
               }

               if (is.null(xlim)) {
                  xlimv <- range(size, na.rm=TRUE)
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
                          ylim=ylimv, pch=pch, lty=1, lwd=lwd/4, cex=cex/4, ...)
                  par(new=TRUE)
               }
               plot(x=size, y=varprofiles.mean, type="b", axes=FALSE,
                    main="CV Tuning Profiles of Variable Screening Size", cex.main=0.8,
                    xlab="Cardinals of Top-Ranked Variables Subsets", ylab=paste(txt, " Mean Profiles", sep=""),
                    xlim=xlimv, ylim=ylimv,
                    pch=pch, col=col, lty=lty, lwd=lwd, cex=cex, ...)
               axis(side=1, pos=min(ylimv), at=size, labels=size, col=1, cex.axis=1, cex.axis=1, line=NA)
               axis(side=2, pos=min(xlimv), at=pretty(ylimv), col=1, col.axis=1, cex.axis=1, line=NA)
               segments(x0=size[cv.varset], y0=min(ylimv),
                        x1=size[cv.varset], y1=varprofiles.mean[cv.varset],
                        col=col, lty=2, lwd=lwd)
               if (add.sd) {
                  arrows(size, varprofiles.mean,
                         size, varprofiles.mean - varprofiles.se,
                         length=0.03, angle=90, code=2, col=col, lwd=lwd)
                  arrows(size, varprofiles.mean,
                         size, varprofiles.mean + varprofiles.se,
                         length=0.03, angle=90, code=2, col=col, lwd=lwd)
               }
               if (add.legend) {
                  legend("top", xpd=TRUE, inset=0, legend=c("Sample Mean", "Std. Error"),
                         pch=pch, col=col, lty=lty, lwd=lwd, cex=cex, pt.cex=cex/2)
               }
            } else {
               if (!is.null(main)) {
                  par(mfrow=c(1, 1), oma=c(0, 0, 3, 0), mar=c(4.0, 3.0, 4.0, 3.0), mgp=c(1.5, 0.5, 0))
               } else {
                  par(mfrow=c(1, 1), oma=c(0, 0, 0, 0), mar=c(4.0, 3.0, 4.0, 3.0), mgp=c(1.5, 0.5, 0))
               }
            }

            if (object$onese) {
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
                       ylim=ylims, pch=pch, lty=1, lwd=lwd/4, cex=cex/4, ...)
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
            if (add.legend) {
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
                        add.sd=add.sd, add.legend=add.legend, add.profiles=add.profiles,
                        pch=pch, col=col, lty=lty, lwd=lwd, cex=cex)
         } else if (device == "PS") {
            path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
            file <- paste(file, ".ps", sep="")
            cat("\nOUTPUT: \n")
            cat("Filename : ", file, "\n")
            cat("Directory: ", path, "\n")
            postscript(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, horizontal=horizontal)
            cat("Device: ",  dev.cur(), "\n")
            profileplot(object=object, main=main, xlim=xlim, ylim=ylim,
                        add.sd=add.sd, add.legend=add.legend, add.profiles=add.profiles,
                        pch=pch, col=col, lty=lty, lwd=lwd, cex=cex)
            dev.off()
         } else if (device == "PDF") {
            path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
            file <- paste(file, ".pdf", sep="")
            cat("\nOUTPUT: \n")
            cat("Filename : ", file, "\n")
            cat("Directory: ", path, "\n")
            pdf(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
            cat("Device: ",  dev.cur(), "\n")
            profileplot(object=object, main=main, xlim=xlim, ylim=ylim,
                        add.sd=add.sd, add.legend=add.legend, add.profiles=add.profiles,
                        pch=pch, col=col, lty=lty, lwd=lwd, cex=cex)
            dev.off()
         } else {
            stop("Currently allowed display devices are \"PS\" (Postscript) or \"PDF\" (Portable Document Format). Exiting ... \n\n")
         }

      } else {

         cat("No CV here, so no cross-validated tuning profile to plot \n")

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
#                    plot_boxtraj (object,
#                                  main=NULL,
#                                  toplot=object$cvfit$cv.used,
#                                  col.cov, lty.cov, lwd.cov,
#                                  col=1, lty=1, lwd=0.5, cex=0.5,
#                                  add.legend=FALSE, text.legend=NULL,
#                                  nr=NULL, nc=NULL,
#                                  device=NULL, file="Trajectory Plots", path=getwd())
#                                  horizontal=FALSE, width=8.5, height=11, ...)
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

plot_boxtraj <- function(object,
                         main=NULL,
                         toplot=object$cvfit$cv.used,
                         col.cov, lty.cov, lwd.cov,
                         col=1, lty=1, lwd=0.5, cex=0.5,
                         add.legend=FALSE, text.legend=NULL,
                         nr=NULL, nc=NULL,
                         device=NULL, file="Trajectory Plots", path=getwd(),
                         horizontal=FALSE, width=8.5, height=11, ...) {

   if (!inherits(object, 'sbh'))
      stop("Primary argument must be an object of class 'sbh' \n")

   if (object$success) {

      boxtrajplot <- function(object,
                              main,
                              toplot,
                              col.cov, lty.cov, lwd.cov,
                              col, lty, lwd,
                              cex, add.legend, text.legend,
                              nr, nc, ...) {
         p <- length(toplot)
         varnames <- colnames(object$X)
         if (is.null(nc))
            nc <- 3
         if (is.null(nr)) {
            if (p %% nc == 0) {
               nr <- p%/%nc + 2
            } else {
               nr <- ((p+(1:nc))[which((p+(1:nc)) %% nc == 0)])%/%nc + 2
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
                 main=paste(varnames[toplot[j]], " covariate trajectory", sep=""), cex.main=cex,
                 xlim=range(0,1),
                 #ylim=range(object$X[,toplot[j]], na.rm=TRUE),
                 xlab="Box Mass",
                 ylab="Covariate Range", ...)
         }
         if (add.legend)
            legend("bottomleft", inset=0.01, legend=text.legend, cex=cex)
         par(mfg=c(nr-1, 1))
         plot(object$cvfit$cv.stats$mean$Support,
              object$cvfit$cv.stats$mean$Support,
              type='s', col=col, lty=lty, lwd=lwd,
              main="Box support trajectory", cex.main=cex,
              xlim=range(0,1),
              ylim=range(0,1),
              xlab="Box Mass",
              ylab=expression(paste("Support (", beta, ")", sep="")), ...)
         if (add.legend)
            legend("bottomright", inset=0.01, legend=text.legend, cex=cex)
         par(mfg=c(nr-1, 2))
         plot(object$cvfit$cv.stats$mean$Support,
              object$cvfit$cv.stats$mean$MEFT,
              type='s', col=col, lty=lty, lwd=lwd,
              main="MEFT trajectory", cex.main=cex,
              xlim=range(0,1),
              ylim=range(0, object$cvfit$cv.stats$mean$MEFT, na.rm=TRUE),
              xlab="Box Mass",
              ylab="Time", ...)
         if (add.legend)
            legend("bottomright", inset=0.01, legend=text.legend, cex=cex)
         par(mfg=c(nr-1, 3))
         plot(object$cvfit$cv.stats$mean$Support,
              object$cvfit$cv.stats$mean$MEFP,
              type='s', col=col, lty=lty, lwd=lwd,
              main="MEFP trajectory", cex.main=cex,
              xlim=range(0,1),
              ylim=range(0,1),
              xlab="Box Mass",
              ylab="Probability", ...)
         if (add.legend)
            legend("bottomright", inset=0.01, legend=text.legend, cex=cex)
         par(mfg=c(nr, 1))
         plot(object$cvfit$cv.stats$mean$Support,
              object$cvfit$cv.stats$mean$LHR,
              type='s', col=col, lty=lty, lwd=lwd,
              main="LHR trajectory", cex.main=cex,
              xlim=range(0,1),
              ylim=range(0, object$cvfit$cv.stats$mean$LHR, na.rm=TRUE),
              xlab="Box Mass",
              ylab=expression(paste("Log-Hazard Ratio (", lambda,")", sep="")), ...)
         if (add.legend)
            legend("top", inset=0.01, legend=text.legend, cex=cex)
         par(mfg=c(nr, 2))
         plot(object$cvfit$cv.stats$mean$Support,
              object$cvfit$cv.stats$mean$LRT,
              type='s', col=col, lty=lty, lwd=lwd,
              main="LRT trajectory", cex.main=cex,
              xlim=range(0,1),
              ylim=range(0, object$cvfit$cv.stats$mean$LRT, na.rm=TRUE),
              xlab="Box Mass",
              ylab=expression(paste("Log-rank test (", chi^2 ,")", sep="")), ...)
         if (add.legend)
            legend("top", inset=0.01, legend=text.legend, cex=cex)
         par(mfg=c(nr, 3))
         plot(object$cvfit$cv.stats$mean$Support,
              object$cvfit$cv.stats$mean$CER,
              type='s', col=col, lty=lty, lwd=lwd,
              main="CER trajectory", cex.main=cex,
              xlim=range(0,1),
              ylim=range(0,1),
              xlab="Box Mass",
              ylab=expression(paste("1-C (", theta,")", sep="")), ...)
         if (add.legend)
            legend("top", inset=0.01, legend=text.legend, cex=cex)
         if (!is.null(main)) {
            mtext(text=main, cex=1, side=3, outer=TRUE)
         }
      }

      if (is.null(device)) {
         cat("Device: ",  dev.cur(), "\n")
         boxtrajplot(object=object,
                     main=main,
                     toplot=toplot,
                     col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
                     col=col, lty=lty, lwd=lwd,
                     cex=cex, add.legend=add.legend, text.legend=text.legend,
                     nr=nr, nc=nc)
      } else if (device == "PS") {
         path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
         file <- paste(file, ".ps", sep="")
         cat("\nOUTPUT: \n")
         cat("Filename : ", file, "\n")
         cat("Directory: ", path, "\n")
         postscript(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, horizontal=horizontal)
         cat("Device: ",  dev.cur(), "\n")
         boxtrajplot(object=object,
                     main=main,
                     toplot=toplot,
                     col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
                     col=col, lty=lty, lwd=lwd,
                     cex=cex, add.legend=add.legend, text.legend=text.legend,
                     nr=nr, nc=nc)
         dev.off()
      } else if (device == "PDF") {
         path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
         file <- paste(file, ".pdf", sep="")
         cat("\nOUTPUT: \n")
         cat("Filename : ", file, "\n")
         cat("Directory: ", path, "\n")
         pdf(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
         cat("Device: ",  dev.cur(), "\n")
         boxtrajplot(object=object,
                     main=main,
                     toplot=toplot,
                     col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
                     col=col, lty=lty, lwd=lwd,
                     cex=cex, add.legend=add.legend, text.legend=text.legend,
                     nr=nr, nc=nc)
         dev.off()
      } else {
         stop("Currently allowed display devices are \"PS\" (Postscript) or \"PDF\" (Portable Document Format). Exiting ... \n\n")
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
#                    plot_boxtrace (object,
#                                   main=NULL,
#                                   xlab="Box Mass", ylab="Covariate Range (centered)",
#                                   toplot=object$cvfit$cv.used,
#                                   center=TRUE, scale=FALSE,
#                                   col.cov, lty.cov, lwd.cov,
#                                   col=1, lty=1, lwd=0.5, cex=0.5,
#                                   add.legend=FALSE, text.legend=NULL,
#                                   device=NULL, file="Covariate Trace Plots", path=getwd(),
#                                   horizontal=FALSE, width=8.5, height=8.5, ...)
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

plot_boxtrace <- function(object,
                          main=NULL,
                          xlab="Box Mass", ylab="Covariate Range (centered)",
                          toplot=object$cvfit$cv.used,
                          center=TRUE, scale=FALSE,
                          col.cov, lty.cov, lwd.cov,
                          col=1, lty=1, lwd=0.5, cex=0.5,
                          add.legend=FALSE, text.legend=NULL,
                          device=NULL, file="Covariate Trace Plots", path=getwd(),
                          horizontal=FALSE, width=8.5, height=8.5, ...) {

   if (!inherits(object, 'sbh'))
      stop("Primary argument must be an object of class 'sbh' \n")

   if (object$success) {

      boxtraceplot <- function(object,
                               main, xlab, ylab,
                               toplot,
                               center, scale,
                               col.cov, lty.cov, lwd.cov,
                               col, lty, lwd,
                               cex, add.legend, text.legend, ...) {
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
         legend("topleft", inset=0.01, legend=varnames[toplot], col=col.cov, lty=lty.cov, lwd=lwd.cov, cex=cex)
         if (center)
            abline(h=0, lty=2, col=1, lwd=0.3, xpd=FALSE)
         if (add.legend)
            legend("bottom", inset=0.01, legend=text.legend, cex=cex)
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
         if (add.legend)
            legend("bottom", inset=0.01, legend=text.legend, cex=cex)
         mtext(text=xlab, cex=cex, side=1, line=1, outer=FALSE)
         mtext(text="Covariates Used", cex=cex, side=2, line=1+maxlength/2, outer=FALSE)
         if (!is.null(main)) {
            mtext(text=main, cex=1, side=3, outer=TRUE)
         }
      }

      if (is.null(device)) {
         cat("Device: ",  dev.cur(), "\n")
         boxtraceplot(object=object,
                      main=main, xlab=xlab, ylab=ylab,
                      toplot=toplot,
                      center=center, scale=scale,
                      col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
                      col=col, lty=lty, lwd=lwd,
                      cex=cex, add.legend=add.legend, text.legend=text.legend)
      } else if (device == "PS") {
         path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
         file <- paste(file, ".ps", sep="")
         cat("\nOUTPUT: \n")
         cat("Filename : ", file, "\n")
         cat("Directory: ", path, "\n")
         postscript(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, horizontal=horizontal)
         cat("Device: ",  dev.cur(), "\n")
         boxtraceplot(object=object,
                      main=main, xlab=xlab, ylab=ylab,
                      toplot=toplot,
                      center=center, scale=scale,
                      col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
                      col=col, lty=lty, lwd=lwd,
                      cex=cex, add.legend=add.legend, text.legend=text.legend)
         dev.off()
      } else if (device == "PDF") {
         path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
         file <- paste(file, ".pdf", sep="")
         cat("\nOUTPUT: \n")
         cat("Filename : ", file, "\n")
         cat("Directory: ", path, "\n")
         pdf(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
         cat("Device: ",  dev.cur(), "\n")
         boxtraceplot(object=object,
                      main=main, xlab=xlab, ylab=ylab,
                      toplot=toplot,
                      center=center, scale=scale,
                      col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
                      col=col, lty=lty, lwd=lwd,
                      cex=cex, add.legend=add.legend, text.legend=text.legend)
         dev.off()
      } else {
         stop("Currently allowed display devices are \"PS\" (Postscript) or \"PDF\" (Portable Document Format). Exiting ... \n\n")
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
#                    plot_boxkm (object,
#                                main=NULL,
#                                xlab="Time", ylab="Probability",
#                                precision=1e-3, mark=3,
#                                col=2, lty=1, lwd=0.5, cex=0.5,
#                                steps=1:object$cvfit$cv.nsteps,
#                                nr=3, nc=4,
#                                device=NULL, file="Survival Plots", path=getwd(),
#                                horizontal=TRUE, width=11, height=8.5, ...)
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

plot_boxkm <- function(object,
                       main=NULL,
                       xlab="Time", ylab="Probability",
                       precision=1e-3, mark=3,
                       col=2, lty=1, lwd=0.5, cex=0.5,
                       steps=1:object$cvfit$cv.nsteps,
                       nr=3, nc=4,
                       device=NULL, file="Survival Plots", path=getwd(),
                       horizontal=TRUE, width=11, height=8.5, ...) {

   if (!inherits(object, 'sbh'))
      stop("Primary argument must be an object of class 'sbh' \n")

   if (object$success) {

      boxkmplot <- function(object,
                            main, xlab, ylab,
                            precision, mark,
                            col, lty, lwd, cex,
                            steps, nr, nc, ...) {
         if (!is.null(main)) {
            par(mfrow=c(nr, nc), oma=c(0, 0, 3, 0), mar=c(2.5, 2.5, 1.5, 1.5), mgp=c(1.5, 0.5, 0))
         } else {
            par(mfrow=c(nr, nc), oma=c(0, 0, 0, 0), mar=c(2.5, 2.5, 0.0, 1.5), mgp=c(1.5, 0.5, 0))
         }
         y <- object$y
         delta <- object$delta
         L <- object$cvfit$cv.nsteps
         for (l in steps) {
            boxind <- object$cvfit$cv.boxind[l,]
            ng <- length(unique(boxind[!is.na(boxind)]))
            if (ng == 1) {
               boxind <- 1*boxind
            } else {
               boxind <- 2 - 1*boxind
            }
            surv <- survival::survfit(survival::Surv(time=y, event=delta) ~ 1 + boxind)
            if (l == 1) {
               plot(surv, main="", conf.int=TRUE, mark.time=FALSE, mark=NA, lty=2, lwd=lwd, col=col, cex=cex, xlab=xlab, ylab=ylab, ...)
               par(new=TRUE)
               plot(surv, main="", conf.int=FALSE, mark.time=TRUE, mark=mark, lty=lty, lwd=lwd, col=col, cex=cex, xlab=xlab, ylab=ylab, ...)
            } else {
               plot(surv, main="", conf.int=TRUE, mark.time=FALSE, mark=NA, lty=c(2,2), lwd=c(lwd,lwd), col=c(col,1), cex=cex, xlab=xlab, ylab=ylab, ...)
               par(new=TRUE)
               plot(surv, main="", conf.int=FALSE, mark.time=TRUE, mark=mark, lty=c(lty,lty), lwd=c(lwd,lwd), col=c(col,1), cex=cex, xlab=xlab, ylab=ylab, ...)
            }
            legend("topright", inset=0.01, legend=c("outbox", "inbox"), lty=c(lty,lty), lwd=c(lwd,lwd), col=c(1,2), cex=0.9*cex)
            if (object$pv) {
               if (object$cvfit$cv.pval$pval[l] <= precision) {
                  legend("bottom", inset=0.11, col="black", cex=0.9*cex, bty="n",
                         legend=bquote(italic(p) <= .(precision)))
               } else {
                  legend("bottom", inset=0.11, col="black", cex=0.9*cex, bty="n",
                         legend=bquote(italic(p) == .(format(x=object$cvfit$cv.pval$pval[l], scientific=FALSE, digits=4, nsmall=4))))
               }
            }
            legend("bottom", inset=0.01, col="black", cex=0.9*cex, bty="n",
                   legend=substitute(group("", list(paste(italic(LHR) == x, sep="")), ""), list(x=format(x=object$cvfit$cv.stats$mean$LHR[l], digits=3, nsmall=3))))
            legend("bottom", inset=0.06, col="black", cex=0.9*cex, bty="n",
                   legend=substitute(group("", list(paste(italic(LRT) == x, sep="")), ""), list(x=format(x=object$cvfit$cv.stats$mean$LRT[l], digits=3, nsmall=3))))
            legend("bottom", inset=0.16, legend=paste("Step ", l-1, sep=""), col=1, cex=0.9*cex, bty="n")
         }
         if (!is.null(main)) {
            mtext(text=main, cex=1, side=3, outer=TRUE)
         }
      }

      if (is.null(device)) {
         cat("Device: ",  dev.cur(), "\n")
         boxkmplot(object=object,
                   main=main, xlab=xlab, ylab=ylab,
                   precision=precision, mark=mark,
                   col=col, lty=lty, lwd=lwd, cex=cex,
                   steps=steps,
                   nr=nr, nc=nc)
      } else if (device == "PS") {
         path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
         file <- paste(file, ".ps", sep="")
         cat("\nOUTPUT: \n")
         cat("Filename : ", file, "\n")
         cat("Directory: ", path, "\n")
         postscript(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, horizontal=horizontal)
         cat("Device: ",  dev.cur(), "\n")
         boxkmplot(object=object,
                   main=main, xlab=xlab, ylab=ylab,
                   precision=precision, mark=mark,
                   col=col, lty=lty, lwd=lwd, cex=cex,
                   steps=steps,
                   nr=nr, nc=nc)
         dev.off()
      } else if (device == "PDF") {
         path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
         file <- paste(file, ".pdf", sep="")
         cat("\nOUTPUT: \n")
         cat("Filename : ", file, "\n")
         cat("Directory: ", path, "\n")
         pdf(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
         cat("Device: ",  dev.cur(), "\n")
         boxkmplot(object=object,
                   main=main, xlab=xlab, ylab=ylab,
                   precision=precision, mark=mark,
                   col=col, lty=lty, lwd=lwd, cex=cex,
                   steps=steps,
                   nr=nr, nc=nc)
         dev.off()
      } else {
         stop("Currently allowed display devices are \"PS\" (Postscript) or \"PDF\" (Portable Document Format). Exiting ... \n\n")
      }

   } else {

      cat("Either the covariate screening or the Survival Bump Hunting modeling failed for this dataset.\n
           So, there is nothing to plot here.\n")

   }

   invisible()

}
#===============================================================================================================================#


