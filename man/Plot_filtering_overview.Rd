\name{Plot_filtering_overview}
\alias{Plot_filtering_overview}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot cells filtered by each filtering steps.
}
\description{
Plotting a graph to show cells filtered by all four filtering steps, including "quality control", "zeroRatio", "sgRNA efficiency" and "phenotype capture".
}
\usage{
Plot_filtering_overview(sample_info_gene_original = c(), sample_info_gene_qc = c(), nonzeroRatio = c(), sample_info_gene_qc_zr_se = c(), sample_info_gene_qc_zr_se_pc = c(), plot_path = "~/overview_of_cell_filterings.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample_info_gene_original}{
  A character vector showing the knockout gene of each sample without any filtering steps.
}
  \item{sample_info_gene_qc}{
  A character vector showing the knockout gene of each sample after filtering step "quality control". It obtains from the result of "singleCellCRISPRscreen_qc"
}
  \item{nonzeroRatio}{
  A dataframe showing the ratio of nonzero knockout value in control cells for each knockout. It obtains from the result of "cellFilteringAndKOefficiencyCalculating".
}
  \item{sample_info_gene_qc_zr_se}{
  A character vector showing the knockout gene of each sample after filtering steps "quality control", "zeroRatio" and "sgRNA efficiency". It obtains from the result of "cellFilteringAndKOefficiencyCalculating".
}
  \item{sample_info_gene_qc_zr_se_pc}{
  A character vector showing the knockout gene of each sample after filtering steps "quality control", "zeroRatio", "sgRNA efficiency" and "phenotype capture". It obtains from the result of "cellFilteringAndKOefficiencyCalculating".
}
  \item{plot_path}{
  The path of the graph.
}
}
\value{
A dataframe showing the cell number filtered by corresponding filtering step.
}
\author{
Bin Duan
}
