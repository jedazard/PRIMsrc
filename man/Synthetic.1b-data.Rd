\name{Synthetic.1b}

\docType{data}

\alias{Synthetic.1b}

\title{Synthetic Dataset #1b: \eqn{p < n} case}

\description{
  Dataset from simulated regression survival model #1b as described in Dazard et al. (2015).
  Here, the regression function uses all of the predictors, which are also part of the design matrix. 
  In this example, the signal is limited to a box-shaped region R of the predictor space. 
  Survival time was generated from an exponential model with rate parameter \eqn{\lambda} (and mean 1/\eqn{\lambda})
  according to a Cox-PH model with hazard exp(eta), where eta(.) is the regression function.
  Censoring indicator were generated from a uniform distribution on [0, 3].
  In this synthetic example, all covariates are continuous, i.i.d. from a multivariate uniform distribution on [0, 1].
}

\usage{
    data("Synthetic.1b", package="PRIMsrc")
}

\format{
  Each dataset consists of a \code{numeric} \code{matrix} containing \eqn{n=250} observations (samples) 
  by rows and \eqn{p=3} variables by columns, not including the censoring indicator and (censored) time-to-event variables.
  It comes as a compressed Rda data file.
}

\source{See simulated survival model #1b in Dazard et al., 2015.}

\author{
    \itemize{
        \item "Jean-Eudes Dazard, Ph.D." \email{jean-eudes.dazard@case.edu}
        \item "Michael Choe, M.D." \email{mjc206@case.edu}
        \item "Michael LeBlanc, Ph.D." \email{mleblanc@fhcrc.org}
        \item "Alberto Santana, MBA." \email{ahs4@case.edu}
        \item "J. Sunil Rao, Ph.D." \email{Rao@biostat.med.miami.edu}
    }
    Maintainer: "Jean-Eudes Dazard, Ph.D." \email{jean-eudes.dazard@case.edu}
}

\references{
    \itemize{
        \item Dazard J-E. and Rao J.S. (2021a).
              "\emph{Variable Selection Strategies for High-Dimensional Recursive Peeling-Based Survival Bump Hunting Models.}"
              (in prep).
        \item Dazard J-E. and Rao J.S. (2021b).
              "\emph{Group Bump Hunting by Recursive Peeling-Based Methods: Application to Survival/Risk Predictive Models.}"
              (in prep).
        \item Dazard J-E., Choe M., Pawitan Y., and Rao J.S. (2021c).
              "\emph{Identification and Characterization of Informative Prognostic Subgroups by Survival Bump Hunting.}"
              (in prep).
        \item Rao J.S., Huilin Y., and Dazard J-E. (2020).
              "\emph{Disparity Subtyping: Bringing Precision Medicine Closer to Disparity Science.}"
              Cancer Epidemiology Biomarkers & Prevention, 29(6 Suppl):C018.
        \item Yi C. and Huang J. (2017).
              "\emph{Semismooth Newton Coordinate Descent Algorithm for Elastic-Net Penalized Huber Loss Regression and Quantile Regression}."
              J. Comp Graph. Statistics, 26(3):547-557.
        \item Dazard J-E., Choe M., LeBlanc M., and Rao J.S. (2016).
              "\emph{Cross-validation and Peeling Strategies for Survival Bump Hunting using Recursive Peeling Methods.}"
              Statistical Analysis and Data Mining, 9(1):12-42.
        \item Dazard J-E., Choe M., LeBlanc M., and Rao J.S. (2015).
              "\emph{R package PRIMsrc: Bump Hunting by Patient Rule Induction Method for Survival, Regression and Classification.}"
              In JSM Proceedings, Statistical Programmers and Analysts Section. Seattle, WA, USA.
              American Statistical Association IMS - JSM, p. 650-664.
        \item Dazard J-E., Choe M., LeBlanc M., and Rao J.S. (2014).
              "\emph{Cross-Validation of Survival Bump Hunting by Recursive Peeling Methods.}"
              In JSM Proceedings, Survival Methods for Risk Estimation/Prediction Section. Boston, MA, USA.
              American Statistical Association IMS - JSM, p. 3366-3380.
    }
}

\section{Acknowledgments}{
    This work made use of the High Performance Computing Resource in the Core Facility for Advanced Research Computing at Case Western Reserve University. 
    This project was partially funded by the National Institutes of Health NIH - National Cancer Institute (R01-CA160593) to J-E. Dazard and J.S. Rao.
}

\keyword{datasets}

