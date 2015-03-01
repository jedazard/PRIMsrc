\name{cv.ave.peel}

\alias{cv.ave.peel}

\title{Internal Subroutine (never to be called by end-user)}

\description{Internal Subroutine of:
    \code{\link[PRIMsrc]{cv.ave.fold}}
}

\usage{
    cv.ave.peel(traindata, trainstatus, traintime, 
                testdata, teststatus, testtime, 
                probval, timeval, 
                varsign, selected, initcutpts,
                K, arg, seed)
}

\arguments{
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
  \item{varsign}{
    \code{numeric} \code{vector} in \{-1,+1\} of directions of peeling for all variables.
  } 
  \item{selected}{
    \code{numeric} \code{vector} giving the selected variable by regularized (Elastic-Net) Cox-regression.
  }
  \item{initcutpts}{
    \code{numeric} \code{vector} of initial box boundaries for each dimension.
  }
  \item{probval}{
    \code{Numeric} scalar of the survival probability at which we want to get the endpoint box survival time.
  }
  \item{timeval}{
    \code{Numeric} scalar of the survival time at which we want to get the endpoint box survival probability.
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

