\name{cbindlist}

\alias{cbindlist}

\title{Internal Subroutine (never to be called by end-user)}

\description{Internal Subroutine of:
    \code{\link[PRIMsrc]{cv.pval}} and
    \code{\link[PRIMsrc]{cv.comb.box}}.
    It is used to bind a list of matrices by columns (even of different number of rows)
    into a single matrix.
}

\usage{
    cbindlist(list, trunc)
}

\arguments{
  \item{list}{
    \code{List} of \code{numeric} \code{matrices}.
  }
  \item{trunc}{
    \code{Integer} scalar for truncation to the same number of rows.
  }
}
