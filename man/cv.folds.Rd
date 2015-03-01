\name{cv.folds}

\alias{cv.folds}

\title{Internal Subroutine (never to be called by end-user)}

\description{Internal Subroutine of:
    \code{\link[PRIMsrc]{cv.ave.fold}} and
    \code{\link[PRIMsrc]{cv.comb.fold}}
}

\usage{
    cv.folds(n, K, seed = NULL)
}

\arguments{
  \item{n}{
    \code{Integer} scalar giving the number of observations to be split into groups.
  }
  \item{K}{
    \code{Integer} giving the number of groups into which the observations should be randomly split.
    Setting K also specifies the type of folds to be generated.
    Possible types are 'random cross-validation', 'leave-one-out cross-validation', or no cross-validation:
    \itemize{
        \item{\code{K}} in \{2,...,\eqn{n-1}\} yields random cross-validation.
        \item{\code{K}} = 1 yields no cross-validation.
        \item{\code{K}} = \eqn{n} yields leave-one-out cross-validation.
    }
  }
  \item{seed}{
    \code{Integer} scalar seed for RNG of random splitting of the data.
  }
}
