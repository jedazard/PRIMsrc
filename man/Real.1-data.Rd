\name{Real.1}

\docType{data}

\alias{Real.1}

\title{Real Dataset #1: Clinical Dataset (\eqn{p < n} case)}

\description{
  Publicly available HIV clinical data from the Women's Interagency HIV cohort Study (WIHS). 
  The entire study enrolled 1164 women. Inclusion criteria of the study are: women at enrolment must be
  (i) alive, (ii) HIV-1 infected, and (iii) free of clinical AIDS symptoms. Women were followed until the 
  first of the following occurred: (i) treatment initiation (HAART), (ii) AIDS diagnosis, (iii) death, 
  or administrative censoring. The studied outcomes were the competing risks "AIDS/Death (before HAART)" 
  and "Treatment Initiation (HAART)". However, for simplification purposes, only the first of the two competing events 
  (i.e. the time to AIDS/Death), was used. Likewise, for simplification in this clinical dataset example, 
  only complete cases were used. Variables included history of Injection Drug Use ("IDU") 
  at enrollment, African American ethnicity ('Race'), age ('Age'), and baseline CD4 count ('CD4') for a total of 
  \eqn{p=4} clinical covariates. The question in this dataset example was whether it is possible to achieve a 
  prognostication of patients for AIDS and HAART. See Bacon et al. (2005) and the 
  \href{https://statepi.jhsph.edu/wihs/wordpress/}{WIHS} website for more details.
}

\usage{
    data("Real.1", package="PRIMsrc")
}

\format{
  Dataset consists of a \code{numeric} \code{data.frame} containing \eqn{n=485} complete observations (samples) 
  by rows and \eqn{p=4} clinical covariates by columns, not including the censoring indicator and (censored) time-to-event variables.
  It comes as a compressed Rda data file.
}

\source{See real data application in Dazard et al., 2015.}

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
        \item Rao J.S., Huilin Y., and Dazard J-E. (2021d).
              "\emph{Disparity Subtyping: Bringing Precision Medicine Closer to Disparity Science.}"
              (in prep).
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
        \item Bacon M.C., von Wyl V., Alden C. et al. (2005).
              "\emph{The Women's Interagency HIV Study: an observational cohort brings clinical sciences to the bench.}" 
              Clin. Diagn. Lab. Immunol., 12(9):1013-1019.
    }
}

\section{Acknowledgments}{
    This work made use of the High Performance Computing Resource in the Core Facility for Advanced Research Computing at Case Western Reserve University. 
    This project was partially funded by the National Institutes of Health NIH - National Cancer Institute (R01-CA160593) to J-E. Dazard and J.S. Rao.
}

\keyword{Real Dataset}
\keyword{AIDS Prognostication}

\seealso{
    \itemize{
        Women's Interagency HIV cohort Study website: \url{https://statepi.jhsph.edu/wihs/wordpress/}
    }
}
