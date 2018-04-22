\name{Cell_filtering}
\alias{Cell_filtering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filtering cells and calculating the overall knockout efficiency.
}
\description{
Filtering cells with three aspects. (1) Filtering out cells with an extremely large proportions of zero knockout expression in the control cells. (2) Reducing false positive rate of gene knockout by considering sgRNA efficiency. (3) filtering out knockout cells without sufficient cell number to capture the corresponding perturbation phenotype.
On the other hand, this function calculates the overall knockout efficiency.
}
\usage{
Cell_filtering(expression_profile, sample_info_gene, sample_info_sgRNA, nonzero = 0.01, grna_cell_num = 10, fold_change = 0.5, cellNum = 30, plot = FALSE, plot_path_zeroRatio = "~/nonzero_ratio.pdf", plot_path_sgRNAefficiency = "~/sgRNA_efficiency.pdf", plot_path_capturePhenotype = "~/phenotype_capture.pdf", plot_path_KO_efficiency = "~/KO_efficiency.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile}{
  A dataframe showing the expression profile after performing cell quality control.
}
  \item{sample_info_gene}{
  A character vector showing the knockout gene of each sample. The vector's column name is the sample name.
}
  \item{sample_info_sgRNA}{
  A character vector showing the sgRNA detected in each sample. The vector's column name is the sample name.
}
  \item{nonzero}{
  Its value is between zero and one. It means how much percentage the nonzero value of a knockout gene in control cells should have. The default is 0.01.
}
  \item{grna_cell_num}{
  The minimum number of cells detected with same sgRNA to calculate sgRNA efficiency. The default is 10.
}
  \item{fold_change}{
  The fold change is set default 0.5 between case and control.
}
  \item{cellNum}{
  The minimum cell number for corresponding phenotype capture. The default is 30.
}
  \item{plot}{
  FALSE by default. If TRUE, plot the graph.
}
  \item{plot_path_zeroRatio}{
  The path of graph showing nonzero ratio of each knockout in control cells. It works only when "plot" is TRUE.
}
  \item{plot_path_sgRNAefficiency}{
  The path of graph showing each sgRNA efficiecy. It works only when "plot" is TRUE.
}
  \item{plot_path_capturePhenotype}{
  The path of graph showing each knockout's cell number. It works only when "plot" is TRUE.
}
  \item{plot_path_KO_efficiency}{
  The path of graph showing each knockout' perturbation efficiency. It works only when "plot" is TRUE.
}
}
\value{
\item{expression_profile_qc_zr_se}{The expression profile after performing filtering steps: "quality control", "zeroRatio" and "sgRNA efficiency".}
\item{sample_info_gene_qc_zr_se}{A character vector showing the knockout gene of each sample after performing filtering steps: "quality control", "zeroRatio" and "sgRNA efficiency".}
\item{expression_profile_qc_zr_se_pc}{The expression profile after performing filtering steps: "quality control", "zeroRatio", "sgRNA efficiency" and "phenotype capture".}
\item{sample_info_gene_qc_zr_se_pc}{A character vector showing the knockout gene of each sample after performing filtering steps: "quality control", "zeroRatio", "sgRNA efficiency" and phenotype capture.}
\item{KO_efficiency}{A numeric vector showing the perturbation efficiency of each knockout}
\item{grna_efficiency}{A numeric vector showing the efficiency of each sgRNA}
\item{nonzeroRatio}{A dataframe showing the ratio of nonzero knockout value in control cells for each knockout}
}
\author{
Bin Duan
}
