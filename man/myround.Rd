\name{myround}

\alias{myround}

\title{Internal Subroutine (never to be called by end-user)}

\description{
    Internal Subroutine of:
    \code{\link[PrimSRC]{sbh}}. 
    Go to the nearest digit rounding.
    Note that for rounding off a 0.5, the "go to the even digit" standard is NOT used here.
}

\usage{
    myround(x, digits = 0)
}

\arguments{
  \item{x}{
    \code{Numeric} \code{vector}
  }
  \item{digits}{
    \code{Integer} scalar indicating the number of decimal points (precision) to be used. Defaults to 0.
  }
}
