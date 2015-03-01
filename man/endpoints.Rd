\name{endpoints}

\alias{endpoints}

\title{Internal Subroutine (never to be called by end-user)}

\description{Internal Subroutine of:
    \code{\link[PRIMsrc]{cv.ave.peel}} and
    \code{\link[PRIMsrc]{cv.comb.peel}}. 
    Returns the maximum time value and the corresponding minimum survival probability
    for every box of the trajectory (box peeling sequence). Also returns the time value
    and/or the corresponding survival probability depending on what is specified by the user.
}

\usage{
    endpoints(ind, timemat, probmat, timeval, probval)
}

\arguments{
  \item{ind}{
    \code{Numeric} \code{vector} taking on values in \{0,1\}. 
  }
  \item{timemat}{
    \code{Numeric} \code{matrix} of survival times 
    for each peeling step (by rows) and test samples (by columns).
  }
  \item{probmat}{
    \code{Numeric} \code{matrix} of survival probabilities 
    for each peeling step (by rows) and test samples (by columns).
  }
  \item{probval}{
    \code{Numeric} scalar of the survival probability at which we want to get the endpoint box survival time. Defaults to \code{NULL}.
  }
  \item{timeval}{
    \code{Numeric} scalar of the survival time at which we want to get the endpoint box survival probability. Defaults to \code{NULL}.
  }
}
