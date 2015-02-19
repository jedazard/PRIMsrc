\name{cv.pval}

\alias{cv.pval}

\title{Internal Subroutine (never to be called by end-user)}

\description{Internal Subroutine of:
    \code{\link[PrimSRC]{sbh}}
}

\usage{
    cv.pval(x, times, status, 
            cvtype, 
            varsign, varnz, initcutpts, 
            A, K, arg, obs.chisq, 
            parallel, conf)
}

\arguments{
  \item{x}{
    \code{Numeric} \code{matrix} of original covariates.
  }
  \item{times}{
    \code{Numeric} \code{vector} of original survival times.
  }
  \item{status}{
    \code{Logical} \code{vector} of original event indicators.
  }
  \item{cvtype}{
    \code{Character} \code{vector} describing the cross-validation technique in \{"none", "averaged", "combined"\}.
  }
  \item{varsign}{
    \code{numeric} \code{vector} in \{-1,+1\} of directions of peeling for all variables.
  } 
  \item{varnz}{
    \code{numeric} \code{vector} giving the selected variable by regularized (Elastic-Net) Cox-regression.
  }
  \item{initcutpts}{
    \code{numeric} \code{vector} of initial box boundaries for each dimension.
  }
  \item{A}{
    Positive \code{integer} scalar of the number of permutations for the computation of cross-validated p-values. Defaults to 1000.
  }
  \item{K}{
    Positive \code{integer} scalar of the number of folds for the cross-validation procedure.
  }
  \item{arg}{
    \code{Character} \code{vector} describing the parameters to use:
    \itemize{
        \item{\code{alpha}} = fraction to peel off at each step.
        \item{\code{beta}} = minimum support size resulting from the peeling sequence.
        \item{\code{minn}} = minimum number of observation in a box.
        \item{\code{L}} = fixed peeling length.
        \item{\code{peelcriterion}} in \{"hr" (LHR), "lr" (LRT)\}.
    }
  }
  \item{obs.chisq}{
    \code{Numeric} \code{vector} of observed LRT values at each peeling step.
  }
  \item{parallel}{
    \code{Logical}. Is parallel computing to be performed? Optional. Defaults to \code{FALSE}.
  }
  \item{conf}{
    \code{List} of parameters for cluster configuration.
    Inputs for R package \pkg{parallel} function \code{makeCluster} (R package \pkg{parallel}) for cluster setup.
    Optional, defaults to \code{NULL}. See details for usage.
  }
}



