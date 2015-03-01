\name{list2mat}

\alias{list2mat}

\title{Internal Subroutine (never to be called by end-user)}

\description{
    Internal Subroutine of:
    \code{\link[PRIMsrc]{sbh}}, 
    \code{\link[PRIMsrc]{cv.null}} and 
    \code{\link[PRIMsrc]{lapply.mat}}. 
    Used to bind a list of vectors (even of different length) 
    by rows into a single matrix.
}
             
\usage{
    list2mat(list, trunc = NULL, sub = NULL, fill = NA)
}

\arguments{
  \item{list}{
    \code{List} of \code{numeric} \code{vectors}.
  }
  \item{trunc}{
    \code{Integer} scalar for truncation to the same length.
  }
  \item{sub}{
    \code{Integer} scalar for reaching out to one sub-level of the list. Defaults to \code{NULL}.
  }
  \item{fill}{
    \code{Numeric} scalar to fill in the missing values. Defaults to \code{NA}.  
  }
}
