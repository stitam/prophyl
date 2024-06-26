% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trip_to_tip.R
\name{trip_to_tip}
\alias{trip_to_tip}
\title{Trip to tip analysis}
\usage{
trip_to_tip(tree_tbl, target, focus = NULL, min_sc_prob = 0)
}
\arguments{
\item{tree_tbl}{tibble}
}
\description{
Trip to tip analysis counts the number of ancestral state changes between a
common ancestor and each tip which descends from this common ancestor. The 
function currently counts the number of state changes from root.
}
\examples{
\dontrun{
res <- trip_to_tip(tree_tbl, "k_serotype", min_sc_prob = 0.5)
}
}
