\name{Get_offtarget}
\alias{Get_offtarget}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluating the offtarget effect of sgRNAs.
}
\description{
sgRNA offtarget effect may exist for these novel type of data due to the application of CRISPR knockout screening technique. By our analysis for such data, offtargets were not found, so you may skip this step if you want to save time or you don't have sgRNA sequence . The pipeline works fine even though this step isn't performed.
}
\usage{
Get_offtarget(offTarget_results, expression_profile, sample_info_gene, sample_info_sgRNA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{offTarget_results}{
  The result calculated by another R package CRISPRseek.
}
  \item{expression_profile}{
  A dataframe showing the expression profile after all the filterings.
}
  \item{sample_info_gene}{
  A character vector showing the knockout gene of each sample after all the filterings.
}
  \item{sample_info_sgRNA}{
  A character vector showing the sgRNA detected in each sample after all the filterings.
}
}
\details{
MASCOT integrates sgRNA sequence information with its corresponding knockout gene expression to determine whether this sgRNA has induced the off-target effect.
1). CRISPRseek is performed to predict possible off-targets based on the sgRNA sequence information.
2). The correlations of expression values between the corresponding knockout gene and the possible off-targets are calculated for the case and control respectively.
3). If a significantly increasing of correlations between case and control is detected, the possible off-target effect for this knockout is reported in MASCOT.
}
\value{
A hash table showing the knockout genes and their corresponding offtargets obtained.
}
\references{
Zhu LJ, Holmes BR, Aronin N, Brodsky MH. CRISPRseek: a bioconductor package to identify target-specific guide RNAs for CRISPR-Cas9 genome-editing systems. PloS one. 2014;9(9):e108424 doi: 10.1371/journal.pone.0108424
}
\author{
Bin Duan
}
