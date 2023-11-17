\name{GateMeClass_annotate}
\alias{GateMeClass_annotate}
\title{Annotation of a cytometry dataset using a marker table}
\description{The annotation module of GateMeClass inputs a dataset in the form of a matrix and a marker table and returns a vector of cell annotations}
\usage{
GateMeClass_annotate(exp_matrix,
    marker_table,
    GMM_parameterization = "V",
    reject_option = FALSE,
    sampling = 0.1,
    RSS = TRUE,
    k = 20,				
    verbose = TRUE,
    narrow_marker_table = FALSE,
    seed = 1)
}
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

\value{
    The function returns a list with the following elements:
    \item{labels}{
        The labels returned by GateMeClass, class = character
    }    
    \item{marker_table}{
        The marker table used by GateMeClass, class = data.frame
    }    
    \item{cell_signatures}{
        The marker signatures of each cell attributed by GateMeClass, class = data.frame
    }
}
\examples{
# For examples refer to https://github.com/simo1c/GateMeClass
}
