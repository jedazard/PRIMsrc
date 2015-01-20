\name{lapply.array}

\alias{lapply.array}

\title{Internal Subroutine (never to be called by end-user)}

\description{Internal Subroutine of:
    \code{\link[PrimSRC]{sbh}} and 
    \code{\link[PrimSRC]{cv.comb.box}}. 
    Computes FUN element-wise on entries of a list of matrices (even of different row or column numbers).
}

\usage{
    lapply.array(X, trunc = NULL, sub = NULL, 
                 fill = NA, MARGIN = 1:2, FUN, ...)
}

\arguments{
  \item{X}{
    \code{List} of \code{numeric} \code{matrices}.
  }
  \item{trunc}{
    \code{Integer} scalar for truncation to the same number of rows.
  }
  \item{sub}{
    \code{Integer} scalar for reaching out to one sub-level of the list. Defaults to \code{NULL}.
  }
  \item{fill}{
    \code{Numeric} scalar to fill in the missing values. Defaults to \code{NA}.  
  }
  \item{MARGIN}{
    \code{Integer} \code{vector} giving the subscripts which the function will be applied over. 
    E.g., for a matrix 1 indicates rows, 2 indicates columns, c(1, 2) indicates rows and columns. 
    Where \code{X} has named dimnames, it can be a \code{character} \code{vector} selecting dimension names.
    Defaults to 1:2.
  }
  \item{FUN}{
    The function to be applied. 
    In the case of functions like +, \%*\%, etc., the function name must be backquoted or quoted.
  }
  \item{\dots}{
    Optional arguments to FUN.
  }
}
