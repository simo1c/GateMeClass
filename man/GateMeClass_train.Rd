\name{GateMeClass_train}
\alias{GateMeClass_train}
\title{Extracting the marker table from a reference annotated dataset}
\description{The training module of GateMeClass inputs an annotated cytometry dataset in the form of an expression matrix and outputs a marker table with a pseudo gating strategy for each cell type}
\usage{
GateMeClass_train(reference,
    labels,
    RSS = T,
    GMM_parameterization = "V",
    verbose = T, 
    seed = 1)
}

\arguments{
  \item{reference}{
  Expression matrix of the reference annotated dataset, class = matrix
  }
  \item{labels}{
  Labels of the reference dataset, class = character
  }
  \item{RSS}{
  RSS (Ranked Set Sampling). It is particularly advised in combination with GMM_parameterization = "V" to have a better resolution of the marker distribution, class = logical
  }
  \item{GMM_parameterization}{
  GMM (Gaussian-Mixture-Model) parameter: "V" (Varying) or "E" (Equal), class = character
  }
  \item{verbose}{
  Show output information
  }
  \item{seed}{
  Seed for randomization
  }
}
\examples{
For usage examples refer to https://github.com/simo1c/GateMeClass
}
