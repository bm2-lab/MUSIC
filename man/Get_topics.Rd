\name{Get_topics}
\alias{Get_topics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtaining the result by performing topic model.
}
\description{
Obtaining the result by performing topic model with selection of the optimal topic number automatically.
}
\usage{
Get_topics(expression_profile_var_gene, sample_info_gene, topic_number_min = 3, topic_number_max = 8, alpha = 0.5, seed_num = 2017, burnin = 1000, thin = 100, iter = 1000, plot = FALSE, plot_path = "~/select_topic_number.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile_var_gene}{
  A dataframe showing the expression profile only for the selected highly dispersion differentially expressed genes.
}
  \item{sample_info_gene}{
  A character vector showing the knockout gene of each sample after all the filterings.
}
  \item{topic_number_min}{
  The minimum candidate topic number.
}
  \item{topic_number_max}{
  The maximum candidate topic number.
}
  \item{alpha}{
  alpha is the weight of specificity score. Its value is between 0 and 1. The default is 0.5.
}
  \item{plot}{
  FALSE by default. If TRUE, plot the graph.
}
  \item{plot_path}{
  The path of the graph you plot. It works only when the parameter "plot" is TRUE.
}
  \item{seed_num}{
  Object of class "integer"; used to set the seed for Gibbs sampling. Default 2017.
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
Object of class "LDA" with the optimal topic number.
}
\references{
Blei D.M., Ng A.Y., Jordan M.I. (2003). Latent Dirichlet Allocation. Journal of Machine Learning Research, 3, 993–1022.
}
\author{
Bin Duan
}
