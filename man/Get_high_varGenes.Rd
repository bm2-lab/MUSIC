\name{Get_high_varGenes}
\alias{Get_high_varGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify highly dispersion differentially expressed genes.
}
\description{
By calculating the difference of dispersion between case and control, MASCOT identifys highly dispersion differentially expressed genes for subsequent analysis.
}
\usage{
Get_high_varGenes(expression_profile, perturb_information, x.low.cutoff = 0.01, y.cutoff = 0, num.bin = 30, plot = FALSE, plot_path = "./get_high_var_genes.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile}{
  A dataframe showing the expression profile after all the filterings.
}
  \item{perturb_information}{
  A character vector showing the perturbation of each sample after all the filterings.
}
  \item{x.low.cutoff}{
  Bottom cutoff on x-axis for identifying variable genes
}
  \item{y.cutoff}{
  Bottom cutoff on y-axis for identifying variable genes
}
  \item{num.bin}{
  Number of bins to use in the scaled analysis (default is 30)
}
  \item{plot}{
  FALSE by default. If TRUE, plot the graphs.
}
  \item{plot_path}{
  The save path of the graph you plot. It works only when the parameter "plot" is TRUE.
}
}
\value{
A dataframe showing the expression profile only for the selected highly dispersion differentially expressed genes.
}
\author{
Bin Duan
}
