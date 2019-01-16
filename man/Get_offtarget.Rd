\name{Get_offtarget}
\alias{Get_offtarget}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluating the offtarget effect of sgRNAs.
}
\description{
sgRNA offtarget effect may exist for these novel type of data due to the application of CRISPR screening technique. For CRISPRi technique, MUSIC won’t consider this step, since CRISPRi knockdown is specific with minimal off-target effects. MUSIC integrates sgRNA sequence information with its corresponding perturbation expression to determine whether the sgRNA has induced an off-target effect. By our analysis for such data, offtargets were not found, so you may skip this step if you want to save time or you don't have sgRNA sequence . The pipeline works fine even though this step isn't performed.
}
\usage{
Get_offtarget(offTarget_results, expression_profile, perturb_information, sgRNA_information)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{offTarget_results}{
  The result calculated by another R package CRISPRseek.
}
  \item{expression_profile}{
  A dataframe showing the expression profile after all the filterings.
}
  \item{perturb_information}{
  A character vector showing the perturbation of each sample after all the filterings.
}
  \item{sgRNA_information}{
  A character vector showing the sgRNA detected in each sample after all the filterings.
}
}
\value{
A hash table showing the perturbations and their corresponding offtargets detected.
}
\references{
Zhu LJ, Holmes BR, Aronin N, Brodsky MH. CRISPRseek: a bioconductor package to identify target-specific guide RNAs for CRISPR-Cas9 genome-editing systems. PloS one. 2014;9(9):e108424 doi: 10.1371/journal.pone.0108424
}
\author{
Bin Duan
}
