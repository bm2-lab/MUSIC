\name{Get_topics}
\alias{Get_topics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performing topic model.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Get_topics(expression_profile_varGene, perturb_information, topic_number = c(4:6), seed_num = 2018, burnin = 1000, thin = 100, iter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile_varGene}{
  A dataframe showing the expression profile only for the selected highly dispersion differentially expressed genes.
}
  \item{perturb_information}{
  A character vector showing the perturbation of each sample after all the filterings.
}
  \item{topic_number}{
  The range of topic number.
}
  \item{seed_num}{
  Object of class "integer"; used to set the seed for Gibbs sampling. Default 2018.
}
  \item{burnin}{
  Object of class "integer"; number of omitted Gibbs iterations at beginning, by default 1000.
}
  \item{thin}{
  Object of class "integer"; number of omitted in-between Gibbs iterations, by default equals iter.
}
  \item{iter}{
  Object of class "integer"; number of Gibbs iterations, by default equals 1000.
}
}
\value{
A list of "LDA"" class with the topic number you choose.
}
\references{
Blei D.M., Ng A.Y., Jordan M.I. (2003). Latent Dirichlet Allocation. Journal of Machine Learning Research, 3, 993–1022.
}
\author{
Bin Daun
}
