\name{list2array}

\alias{list2array}

\title{Internal Subroutine (never to be called by end-user)}

\description{
    Internal Subroutine of:
    \code{\link[PrimSRC]{sbh}}, 
    \code{\link[PrimSRC]{cv.null}} and 
    \code{\link[PrimSRC]{lapply.mat}}. 
    Used to bind a list of matrices (even of different row or column numbers) 
    by third dimension of a single 3D array.
}

\usage{
    list2array(list, trunc = NULL, sub = NULL, fill = NA)
}

\arguments{
  \item{list}{
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
}
