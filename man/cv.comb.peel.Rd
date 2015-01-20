\name{cv.comb.peel}

\alias{cv.comb.peel}

\title{Internal Subroutine (never to be called by end-user)}

\description{Internal Subroutine of:
    \code{\link[PrimSRC]{cv.comb.fold}}
}

\usage{
    cv.comb.peel(x, times, status, 
                 traindata, trainstatus, traintime, 
                 testdata, teststatus, testtime, 
                 K, arg, seed)
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
  \item{traindata}{
      \code{Numeric} \code{matrix} of training covariates.
  }
  \item{trainstatus}{
      \code{Logical} \code{vector} of training event indicators.
  }
  \item{traintime}{
    \code{Numeric} \code{vector} of training survival times.
  }
  \item{testdata}{
        \code{Numeric} \code{matrix} of test covariates.
    }
  \item{teststatus}{
    \code{Logical} \code{vector} of test event indicators.
  }
  \item{testtime}{
      \code{Numeric} \code{vector} of test survival times.
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
  \item{seed}{
    Positive \code{integer} scalar for the seed.
 }

}
