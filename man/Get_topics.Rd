\name{Get_topics}
\alias{Get_topics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performing topic model to single-cell CRISPR screening data.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Get_topics(expression_profile, perturb_information, topic_number = c(4:6), seed_num = 2018, burnin = 0, thin = 500, iter = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile}{
  A dataframe showing the expression profile only for the selected highly dispersion differentially expressed genes.
}
  \item{perturb_information}{
  A character vector showing the perturbation of each sample after all the filterings.
}
  \item{topic_number}{
  The range of topic number. The default is 4 to 6. In most cases, 4 is better.
}
  \item{seed_num}{
  Object of class "integer"; used to set the seed for Gibbs sampling. Default 2018.
}
  \item{burnin}{
  Object of class "integer"; number of omitted Gibbs iterations at beginning, by default 0.
}
  \item{thin}{
  Object of class "integer"; number of omitted in-between Gibbs iterations, by default equals iter.
}
  \item{iter}{
  Object of class "integer"; number of Gibbs iterations, by default equals 500.
}
}
\value{
  \item{models }{A list of "LDA"" class with the topic number you choose.}
  \item{perturb_information }{A character vector showing the perturbation of each sample after all the filterings.}
}
\references{
Blei D.M., Ng A.Y., Jordan M.I. (2003). Latent Dirichlet Allocation. Journal of Machine Learning Research, 3, 993–1022.
}
\author{
Bin Daun
}
