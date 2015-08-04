\name{Real.2}

\docType{data}

\alias{Real.2}

\title{Real Dataset #2: Large Gene Expression Dataset (\eqn{p >> n} case)}

\description{
  Publicly available breast cancer gene expression profiling dataset from the Uppsala cohort study.
  It is entitled: "Genetic Reclassification of Histologic Grade Delineates New Clinical Subtypes of Breast Cancer".   
  The goal of the study was to provide a more objective measure of grade with prognostic benefit for patients with 
  moderately differentiated grade II (G2) tumors. To that end, expression profiles of primary invasive breast tumors 
  were analyzed on microarrays to find a gene expression signature capable of discerning tumors of grade I (G1) and grade III (G3) histology.
  In this dataset, only the Uppsala cohort (\eqn{n=249}) and only the gene expression data was included 
  although other clinical covariates are available as well. It contains \eqn{p=22647} mRNA measurements from 
  the Affymetrix Human Genome U133A Array platform on \eqn{n=249} samples. The data was left 
  after removal of outliers and incomplete observations and after removal of Affymetrix controls. 
  See Ivshina et al. (2005) and Gene Expression Omnibus database repository (Accession number: #GSE4922) for more details.
}

\usage{Real.2}

\format{
  Dataset consists of a \code{numeric} \code{data.frame} containing \eqn{n=177} complete observations (samples) 
  by rows and \eqn{p=22577} variables by columns, including the censoring indicator and (censored) time-to-event variables. 
  It comes as a compressed Rda data file.
}

\source{See real clinical data application in Dazard et al., 2015.}

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
        \item  Ivshina AV, George J, Senko O, Mow B et al. (2006).
              "\emph{Genetic reclassification of histologic grade delineates new clinical subtypes of breast cancer.}
              Cancer Res 66(21):10292-301. PMID: 17079448
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
       Gene Expression Omnibus (GEO) database. Accession number: #GSE4922
       \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4922}
    }
}
