\name{Rank_topic_specific}
\alias{Rank_topic_specific}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtaining the topic-specific ranking list.
}
\description{
MASCOT prioritizes the gene knockout effect either as an overall perturbation impact on the cells, or in a functional topic-specific way. Here, the function calcultes and obtains the topic-specific ranking list.
}
\usage{
Rank_topic_specific(rank_overall_result_detail, alpha = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rank_overall_result_detail}{
  The overall perturbation effect ranking list with the detailed information. It can obtain from "Rank_overall()".
}
  \item{alpha}{
  alpha is the weight of impact degree. Its value is between 0 and 1. The default is 0.5 .
}
}
\details{
For a specific topic, MASCOT prioritizes the perturbation impact by calculating this specific topic probability difference between case and control with the rectifications of overall knockout efficiency and background difference. Intuitively, the ranking list is obtained by evaluating the perturbation effect on this specific topic, where the perturbation should influence this topic as much as possible while keeping other topics as unaffected as possible. 
}
\value{
\item{rank_topic_specific_result_detail}{The topic-specific ranking list with detailed information.}
\item{rank_topic_specific_result_summary}{The topic-specific ranking list with the main information.}
}
\author{
Bin Duan
}
