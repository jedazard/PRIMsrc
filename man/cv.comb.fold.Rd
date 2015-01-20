\name{cv.comb.fold}

\alias{cv.comb.fold}

\title{Internal Subroutine (never to be called by end-user)}

\description{Internal Subroutine of:
    \code{\link[PrimSRC]{cv.comb.box}}
}

\usage{
    cv.comb.fold(x, times, status, K, arg, seed)
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
