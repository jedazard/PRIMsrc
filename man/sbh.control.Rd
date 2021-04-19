\name{sbh.control}

\alias{sbh.control}

\title{Parameters Control Function}

\description{
  End-user function to set ancillary parameters of main end-user function \code{sbh} for fitting a Survival Bump Hunting (SBH) model.
  It is used to set some variable screening parameters, optional formats and outputs of \code{sbh}, as well as internally to tune 
  the scatterplot smoother used for finding cross-validated model selection/tuning profile extremum. 
}

\usage{
  sbh.control(vscons = 0.5, 
              decimals = 2, 
              onese = FALSE, 
              probval = NULL, 
              timeval = NULL, 
              lag = 2, 
              span = 0.10, 
              degree = 2)
}

\arguments{
    \item{vscons}{
        \code{numeric} scalar in [1/\code{K}, 1], specifying the conservativeness of the variable screening (pre-selection) 
        procedure, where 1/\code{K} is the least conservative and 1 is the most. Defaults to 0.5.
    }
    \item{decimals}{
        Positive \code{integer} of the number of user-specified significant decimals to output results. Defaults to 2.
    }
    \item{onese}{
        \code{logical} scalar. Flag for using the 1-standard error rule instead of extremum value of the cross-validation criterion 
        when tuning/optimizing model parameters. Defaults to \code{FALSE}.
    }
    \item{probval}{
        \code{numeric} scalar in [0, 1] of the survival probability at which we want to get the endpoint box survival time. 
        Defaults to \code{NULL} (i.e. maximal survival probability value is used).
    }
    \item{timeval}{
        \code{numeric} scalar of the survival time at which we want to get the endpoint box survival probability.
        Defaults to \code{NULL} (i.e. maximal survival time value is used).
    }
    \item{lag}{
        Positive \code{integer} indicating which lag to use in the lagged and iterated difference function. Defaults to 2. 
    }
    \item{span}{
        \code{numeric} scalar in [0, 1], specifying the degree of smoothing in the internal \code{stats::loess} function. 
        Defaults to 0.10. 
    }
    \item{degree}{
        Positive \code{integer} indicating the degree of the polynomials (normally 1 or 2) to be used in the internal 
        \code{stats::loess} function. Here, degree 0 is not also allowed unlike in \code{stats::loess}). Defaults to 2. 
    }
}

\details{
    Example of \code{vscons} values for pre-selection are as follows:
    \itemize{
        \item{`1.0` represents a presence in all the folds (unanimity vote)}
        \item{`0.5` represents a presence in at least half of the folds (majority vote)}
        \item{`1/\code{K}` represents a presence in at least one of the folds (minority vote)}
    }
    Although any value in the interval [1/\code{K},1] is accepted, we recommand using the interval [1/\code{K}, 1/2] to avoid 
    excessive conservativeness. Final variable usage (selection) is done at the time of fitting the Survival Bump Hunting (SBH) model 
    itself using our PRSP algorithm on previously screened variables by collecting those variables that have the maximum occurrence 
    frequency in each peeling step over cross-validation folds and replicates.
}

\value{
  A list of 8 components.
}

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

\note{
    End-user function to be used with \code{sbh}.
}

\seealso{
    \itemize{
        \item \code{\link{sbh}}
        \item \code{diff} (R package \pkg{base})       
        \item \code{loess} (R package \pkg{stats})
    }
}

\keyword{Exploratory Survival/Risk Analysis}
\keyword{Survival/Risk Estimation & Prediction}
\keyword{Non-Parametric Method}
\keyword{Cross-Validation}
\keyword{Bump Hunting}
\keyword{Rule-Induction Method}

