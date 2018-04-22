\name{Rank_overall}
\alias{Rank_overall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtaining the overall perturbation effect ranking list.
}
\description{
MASCOT prioritizes the gene knockout effect either as an overall perturbation effect on the cells, or in a functional topic-specific way. Here, the function calcultes and obtains the overall perturbation effect ranking list.
}
\usage{
Rank_overall(distri_diff, offTarget_hash = hash())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distri_diff}{
  A dataframe showing the topic distribution difference between case and control which can obtain from "Get_distribution_diff()".
}
  \item{offTarget_hash}{
  A hash table showing offtarget information which can obtain from "Get_Offtarget()". The default is a null hash.
}
}
\details{
MASCOT obtains the overall perturbation effect ranking list by calculating the distribution difference of the perturbation between case and control as the Jensen-Shannon divergence (JSD). Then the overall knockout efficiency and background difference are calculated and applied to rectify the distribution difference.
}
\value{
\item{rank_overall_result_summary}{The overall perturbation effect ranking list with the main information.}
\item{rank_overall_result_detail}{The overall perturbation effect ranking list with the detailed information.}
}
\references{
Hinrich Schütze; Christopher D. Manning (1999). Foundations of Statistical Natural Language Processing. Cambridge, Mass: MIT Press. p.304.
}
\author{
Bin Duan
}
