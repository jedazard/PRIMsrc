\name{peel.box}

\alias{peel.box}

\title{Internal Subroutine (never to be called by end-user)}

\description{
    Internal Subroutine of:
    \code{\link[PRIMsrc]{cv.ave.peel}} and 
    \code{\link[PRIMsrc]{cv.comb.peel}}.
}

\usage{
    peel.box(traindata, traintime, trainstatus, 
             varsign, selected, initcutpts, 
             arg, seed)
}

\arguments{
  \item{traindata}{
    \code{Numeric} \code{matrix} of training covariates.
  }
  \item{traintime}{
   \code{Numeric} \code{vector} of training survival times.
  }
  \item{trainstatus}{
    \code{Logical} \code{vector} of training event indicators.
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
    Positive \code{integer} scalar of the user seed to reproduce the results.
 }
}

\details{
    The maximal number of peeling steps (ncut) is determined either by \code{alpha} and \code{beta} 
    metaparameters or the smallest possible fraction of the training data, i.e. \eqn{\frac{1}{n}}:
    \itemize{
        \item{\eqn{ceiling(\frac{\log(\code{beta})} {\log(1 - \code{alpha})})}}  if \code{alpha} and \code{beta} are fixed by user
        \item{\eqn{ceiling(\frac{\log(1/n)} {\log(1 - \code{alpha})})}}          if \code{alpha} is fixed by user and \code{beta} is fixed by data
        \item{\eqn{ceiling(\frac{\log(\code{beta})} {\log(1 - (1/n))})}}         if \code{alpha} is fixed by data and \code{beta} is fixed by user
        \item{\eqn{ceiling(\frac{\log(1/n)} {\log(1 - (1/n))})}}                 if \code{alpha} and \code{beta} are fixed by data
    }
    If \code{L} is not used to specify a fixed number of peeling steps (i.e. \code{NULL}),
    then \code{beta} is used as the stopping rule instead.
}
