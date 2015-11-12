\name{Real.2}

\docType{data}

\alias{Real.2}

\title{Real Dataset #2: Genomic Dataset (\eqn{p >> n} case)}

\description{
  Publicly available lung cancer genomic data from the Chemores Cohort Study. This data is part of an integrated study of mRNA, miRNA 
  and clinical variables to characterize the molecular distinctions between squamous cell carcinoma (SCC) and adenocarcinoma (AC)
  in Non Small Cell Lung Cancer (NSCLC) aside large cell lung carcinoma (LCC). Tissue samples were analysed from a cohort of 
  123 patients who underwent complete surgical resection at the Institut Mutualiste Montsouris (Paris, France) between 
  30 January 2002 and 26 June 2006. In this genomic dataset, the expression levels of Agilent miRNA probes (\eqn{p=939}) 
  were included from the \eqn{n=123} samples of the Chemores cohort. The data contains normalized expression levels. 
  See below the paper by Lazar et al. (2013) and Array Express data repository for complete description of the samples, tissue preparation, 
  Agilent array technology, data normalization, etc. This dataset represents a situation where the number of covariates dominates 
  the number of complete observations, or \eqn{p >> n} case. 
}

\usage{Real.2}

\format{
  Dataset consists of a \code{numeric} \code{data.frame} containing \eqn{n=123} complete observations (samples) 
  by rows and \eqn{p=939} genomic covariates by columns, not including the censoring indicator and (censored) time-to-event variables. 
  It comes as a compressed Rda data file.
}

\source{See real data application in Dazard et al., 2015.}

\references{
    \itemize{
        \item Dazard J-E., Choe M., LeBlanc M. and Rao J.S. (2015).
              "\emph{Cross-validation and Peeling Strategies for Survival Bump Hunting using Recursive Peeling Methods.}"
              Statistical Analysis and Data Mining (in press).
        \item Dazard J-E., Choe M., LeBlanc M. and Rao J.S. (2014).
              "\emph{Cross-Validation of Survival Bump Hunting by Recursive Peeling Methods.}"
              In JSM Proceedings, Survival Methods for Risk Estimation/Prediction Section. Boston, MA, USA.
              American Statistical Association IMS - JSM, p. 3366-3380.
        \item Dazard J-E., Choe M., LeBlanc M. and Rao J.S. (2015).
              "\emph{R package PRIMsrc: Bump Hunting by Patient Rule Induction Method for Survival, Regression and Classification.}"
              In JSM Proceedings, Statistical Programmers and Analysts Section. Seattle, WA, USA.
              American Statistical Association IMS - JSM, (in press).
        \item Dazard J-E. and J.S. Rao (2010).
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

\keyword{Real Dataset}
\keyword{Tumor sample comparisons}

\seealso{
    \itemize{
       Array Express data repository at the European Bioinformatics Institute.
       Accession number: #E-MTAB-1134 (MIR).
       \url{www.ebi.ac.uk/arrayexpress/}
    }
    \itemize{
       CHEMORES Consortium and website.
       \url{www.chemores.ki.se/index.html}
    }
}
