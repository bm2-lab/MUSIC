\name{Get_high_var_genes}
\alias{Get_high_var_genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify highly dispersion differentially expressed genes.
}
\description{
By calculating the difference of dispersion between case and control, MASCOT identifys highly dispersion differentially expressed genes for subsequent analysis.
}
\usage{
Get_high_var_genes(expression_profile, sample_info_gene, x.low.cutoff = 0.0125, x.high.cutoff = 5, y.cutoff = 1, do.spike = FALSE, num.bin = 30, plot = FALSE, plot_path = "~/get_high_var_genes.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile}{
  A dataframe showing the expression profile after all the filterings.
}
  \item{sample_info_gene}{
  A character vector showing the knockout gene of each sample after all the filterings.
}
  \item{x.low.cutoff}{
  Bottom cutoff on x-axis for identifying variable genes
}
  \item{x.high.cutoff}{
  Top cutoff on x-axis for identifying variable genes
}
  \item{y.cutoff}{
  Bottom cutoff on y-axis for identifying variable genes
}
  \item{do.spike}{
  FALSE by default. If TRUE, all genes starting with ^ERCC will be used to normalize the whole data. For most singe cell CRISPR screening data, generally, there are no spike-in genes.
}
  \item{num.bin}{
  Total number of bins to use in the scaled analysis (default is 30)
}
  \item{plot}{
  FALSE by default. If TURE, plot the graphs.
}
  \item{plot_path}{
  The path of graph you plot. It works only when "plot" is TURE.
}
}
\value{
A dataframe showing the expression profile only for the selected highly dispersion differentially expressed genes.
}
\author{
Bin Duan
}
