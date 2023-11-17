\name{GateMeClass_annotate}
\alias{GateMeClass_annotate}
\title{Extracting the marker table from a reference annotated dataset}
\description{It inputs an annotated cytometry dataset in the form of an expression matrix and outputs a marker table with a pseudo gating strategy for each cell type}
\usage{
GateMeClass_annotate(exp_matrix,
    marker_table,
    GMM_parameterization = "V",
    reject_option = F,
    sampling = 0.1,
    RSS = T,
    k = 20,				
    verbose = T,
    narrow_marker_table = F,
    seed = 1)

\arguments{
  \item{exp_matrix}{
  Expression matrix, class = matrix
  }
  \item{marker_table}{
  Manually curated or extracted marker table, class = data.frame
  }
  \item{RSS}{
  RSS (Ranked Set Sampling). It is particularly advised in combination with GMM_parameterization = "V" to have a better resolution of the marker distribution, class = logical
  }
  \item{GMM_parameterization}{
  GMM (Gaussian-Mixture-Model) parameter: "V" (Varying) or "E" (Equal), class = character
  }
  \item{reject_option}{
  This parameter tries to detect cell types not defined in the marker table using MNN algorithm, class = logical
  }
  \item{sampling}{
  GMM (Gaussian-Mixture-Model) parameter: "V" (Varying) or "E" (Equal), class = character
  }
  \item{k}{
  k parameter of k-NN (k-Nearest-Neighbour) used to refine uncertain labels to the most similar already annotated, class = numeric
  }
  \item{narrow_marker_table}{
  format of marker table. TRUE for using a marker table with two columns, "Cell" with the name of cell type and "Gate" with the gating strategy (e.g., CD3+CD4+). FALSE for using a column for each marker (as in the above example), class = logical
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
