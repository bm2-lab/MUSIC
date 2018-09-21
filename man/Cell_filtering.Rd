\name{Cell_filtering}
\alias{Cell_filtering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filtering perturbed cells with invalid edits.
}
\description{
gRNA targets Cas9 to a specific gene locus, but only in 70\%-80\% will generate true loss-of-function of the targeted gene(Sternberg and Doudna, 2015). Therefore, to estimate the ranking of impact of different perturbation, it is necessary to filter cells with invalid edits.
}
\usage{
Cell_filtering(expression_profile, perturb_information, cpu_num = 2, cell_num_threshold = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile}{
  A dataframe showing the expression profile after performing the function of "Cell_qc" and "Data_imputation".
}
  \item{perturb_information}{
  A character vector showing the information of sample after performing the function of "Cell_qc()" and "Data_imputation()".}
  \item{cpu_num}{
  The cpu number for parallel computation. 
}
  \item{cell_num_threshold}{
  A cutoff, the minimal perturbed cell number for each perturbation.
}
}
\value{
  \item{expression_profile }{The expression profile after performing these filtering steps.}
  \item{perturb_information }{The information (perturbation names and sample names) of cells retained after performing these filtering steps.}
  \item{perturb_information_abandon }{The information (perturbation names and sample names) abandoned after performing these filtering steps.}
  \item{filter_record }{The summary of filtering by these steps.}
}
\references{
1. Sternberg, S.H. & Doudna, J.A. Expanding the Biologist's Toolkit with CRISPR-Cas9. Mol Cell 58, 568-574 (2015).
2. Lappalainen, T. et al. Transcriptome and genome sequencing uncovers functional variation in humans. Nature 501, 506-511 (2013).
}
\author{
Bin Duan
}
