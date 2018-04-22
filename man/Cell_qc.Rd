\name{Cell_qc}
\alias{Cell_qc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quality control for single-cell CRISPR screening data.
}
\description{
MASCOT evaluates the cell quality by three factors, i.e. the number of genes detected (default 500), the number of unique molecular identifier induced (default 1500), and the percentage of mitochondrial genes detected (default 10\% among all the detected genes). Only cells with the first two factors above the thresholds and the third factor below the threshold are retained.
}
\usage{
Cell_qc(expression_profile, sample_info_gene, species = "Hs", gene_low = 500, gene_high = 10000, mito_high = 0.1, umi_low = 1500, umi_high = Inf, plot = FALSE, plot_path = "~/quality_control.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile}{
  A dataframe showing the expression profile.
}
  \item{sample_info_gene}{
  A character vector showing the knockout gene of each sample. The vector's column name is the sample name.
}
  \item{species}{
  The species the cells belong to. Currently species="Hs" for homo sapiens or species="Mm" for mus musculus are available.
}
  \item{gene_low}{
  The minimun gene number for cells. Cells with genes below this threshold are filtered(default 500).
}
  \item{gene_high}{
  The maximum gene number for cells. cells with genes above this threshold are filtered(default 10000).
}
  \item{mito_high}{
  The maximum percentage of mitochondrial genes detected(default 10\% among all the detected genes).
}
  \item{umi_low}{
  The minimum number of unique molecular identifier induced (default 1500).
}
  \item{umi_high}{
  The maximun number of unique molecular identifier induced (default Inf)
}
  \item{plot}{
  FALSE by default. If TURE, plot the graph.
}
  \item{plot_path}{
  The path of the graph you plot. It works only when the parameter "plot" is TRUE.
}
}
\value{
\item{expression_profile_qc}{A dataframe showing expression profile after quality control.}
\item{sample_info_gene_qc}{A character vector showing the knockout gene of each sample after quality control.}
}
\references{
http://satijalab.org/seurat/
}
\author{
Bin Duan
}
