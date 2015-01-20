\name{is.empty}

\alias{is.empty}

\title{Internal Subroutine (never to be called by end-user)}

\description{
    Internal Subroutine of:
    \code{\link[PrimSRC]{endpoints}}, 
    \code{\link[PrimSRC]{list2array}}, 
    \code{\link[PrimSRC]{list2mat}} and 
    \code{\link[PrimSRC]{cbindlist}}. 
    Checks if object is empty. Represents the empty array, matrix, or vector of zero dimension or length. 
    Often returned by expressions and functions whose value is undefined.
    It returns a \code{ogical}: \code{TRUE} if its argument is empty and \code{FALSE} otherwise.

}

\usage{
    is.empty(x)
}

\arguments{
  \item{x}{\code{array}, \code{matrix} or \code{vector} of any type.}
}

