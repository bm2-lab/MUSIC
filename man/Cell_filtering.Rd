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
Cell_filtering(expression_profile, perturb_information, cpu_num = 4, cell_num_threshold = 30, umi = 0.01, pvalue = 0.05, vargene_min_num = 5, filtered_rate = 0.9, plot = FALSE, plot_path = "./invalid_rate.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile}{
  A dataframe showing the expression profile after performing the function of "Cell_qc()" and "Data_imputation()".
}
  \item{perturb_information}{
  A character vector showing the information of sample after performing the function of "Cell_qc()" and "Data_imputation()".
}
  \item{cpu_num}{
  The cpu number for parallel computation. The default is 4. Parallel computation is strongly recommeneded to use because this step may take long time without parallel computation.
}
  \item{cell_num_threshold}{
  A cutoff, the minimal perturbed cell number for each perturbation. The default is 30. 
}
  \item{umi}{
  The cutoff of average umi to select the differentially expressed genes. The default is 0.01.
}
  \item{pvalue}{
  The p value to select the differentially expressed genes. The default is 0.05.
}
  \item{vargene_min_num}{
  The minimal number of differentially expressed genes. The default is 5. For a perturbation, if the number of differentially expressed genes are less than 5, this perturbation will be filtered directory. 
}
  \item{filtered_rate}{
  The default is 0.9. For a specific perturbation, if the influenced cells filtered amount to 90\% or higher among all, then such a perturbation was filtered.
}
  \item{plot}{
  FALSE by default. If TRUE, plot the graph to show the ratio of filtered cells for each perturbation.
}
  \item{plot_path}{
  The path of the graph you plot. It works only when the parameter "plot" is TRUE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{expression_profile }{The expression profile after performing these filtering steps.}
  \item{perturb_information }{The information (perturbation names and sample names) of cells retained after performing these filtering steps.}
  \item{perturb_information_abandon }{The information (perturbation names and sample names) abandoned after performing these filtering steps.}
  \item{filter_record }{The summary of filtering by these steps.}
  \item{zero_rate }{The proportion of zero expression value in all cells for each perturbation.}
}
\references{
1. Sternberg, S.H. & Doudna, J.A. Expanding the Biologist's Toolkit with CRISPR-Cas9. Mol Cell 58, 568-574 (2015).
2. Lappalainen, T. et al. Transcriptome and genome sequencing uncovers functional variation in humans. Nature 501, 506-511 (2013).
}
\author{
Bin Duan
}
