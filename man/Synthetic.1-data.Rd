\name{Synthetic.1}

\docType{data}

\alias{Synthetic.1}

\title{\code{Synthetic.1} Dataset}

\description{
  Modeling survival model #1 as described in Dazard et al. (2015) with censoring.
  Here, the regression function uses all of the predictors, which are also part of the design matrix. 
  Survival time was generated from an exponential model with rate parameter \eqn{\lambda} (and mean \eqn{\frac{1}{\lambda}})
  according to a Cox-PH model with hazard exp(eta), where eta(.) is the regression function.
  Censoring indicator were generated from a uniform distribution on [0, 3].
  In this synthetic example, all covariates are continuous, i.i.d. from a multivariate uniform distribution on [0, 1].
}

\usage{Synthetic.1}

\format{
  Each dataset consists of a \code{numeric} \code{matrix} containing \eqn{n=250} observations (samples) by rows
  and \eqn{p=5} variables by columns (\eqn{p < n} case),
  including the censoring indicator and (censored) time-to-event variables.
  Compressed Rda data file.
}

\source{See simulated survival model #1 in Dazard et al., 2015.}

\references{
    \itemize{
        \item Dazard J-E., Choe M., LeBlanc M. and Rao J.S. (2015).
              "\emph{Cross-validation and Peeling Strategies for Survival Bump Hunting using Recursive Peeling Methods.}"
              (Submitted).
        \item Dazard J-E., Choe M., LeBlanc M. and Rao J.S. (2014).
              "\emph{Cross-Validation of Survival Bump Hunting by Recursive Peeling Methods.}"
              In JSM Proceedings, Survival Methods for Risk Estimation/Prediction Section. Boston, MA, USA.
              American Statistical Association IMS - JSM, p. 3366-3380.
        \item Dazard J-E. and J. S. Rao (2010).
              "\emph{Local Sparse Bump Hunting.}"
              J. Comp Graph. Statistics, 19(4):900-92.
    }
}

\author{
    \itemize{
        \item "Jean-Eudes Dazard, Ph.D." \email{jxd101@case.edu}
        \item "Michael Choe, M.D." \email{mjc206@case.edu}
        \item "Michael LeBlanc, Ph.D." \email{mleblanc@fhcrc.org}
        \item "Alberto Santana, MBA." \email{ahs4@case.edu}
    }
    Maintainer: "Jean-Eudes Dazard, Ph.D." \email{jxd101@case.edu}

    Acknowledgments: This project was partially funded by the National Institutes of Health
    NIH - National Cancer Institute (R01-CA160593) to J-E. Dazard and J.S. Rao.
}

\keyword{datasets}

