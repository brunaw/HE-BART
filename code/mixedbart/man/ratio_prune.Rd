% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prune_ratio.R
\name{ratio_prune}
\alias{ratio_prune}
\title{Final ratio for a prune step.}
\usage{
ratio_prune(
  old_tree,
  tree,
  current_node,
  pars,
  p,
  var_in_prune,
  nodes_to_prune,
  results_f
)
}
\arguments{
\item{old_tree}{The previous tree.}

\item{tree}{The current tree.}

\item{current_node}{The current pruned node.}

\item{pars}{The full list of parameters.}

\item{p}{The number of available predictors}

\item{var_in_prune}{The variable that was split in the node chosen to
be pruned.}

\item{nodes_to_prune}{The nodes to prune.}

\item{results_f}{The current results to find the number of second
generation internal nodes.}
}
\value{
The final ratio for the candidate tree.
}
\description{
The final ratio is to be used as the acceptance
criteria in the MCMC of the b-cart model.
}
\author{
Bruna Wundervald, \email{brunadaviesw@gmail.com}.
}
