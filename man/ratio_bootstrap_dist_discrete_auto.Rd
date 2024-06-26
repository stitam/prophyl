% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ratio_bootstrap_dist_discrete_auto.R
\name{ratio_bootstrap_dist_discrete_auto}
\alias{ratio_bootstrap_dist_discrete_auto}
\title{Calculate relative risks for a geographic matrix}
\usage{
ratio_bootstrap_dist_discrete_auto(
  x,
  y,
  geo_mat,
  time_mat2,
  MRCA_mat2,
  geo_mat_ref,
  time_mat2_ref,
  MRCA_mat2_ref,
  Pmax,
  Pmin
)
}
\arguments{
\item{x}{integer; vector of tree tip indices. If tree is not bootstrapped, 
a sequence from 1 until the number of tips within the tree. If tree is
bootstrapped, a bootstrap sample of this sequence (equal length with
replacement).}

\item{y}{integer; vector of tree tip indices. If tree is not bootstrapped, 
a sequence from 1 until the number of tips within the tree. If tree is
bootstrapped, a bootstrap sample of this sequence (equal length with
replacement).}

\item{geo_mat}{matrix; a symmetric geographic matrix, rows and columns are
about tips of the tree, cells contain geographical information, e.g. whether
the isolates for the two tips are from the same country.}

\item{time_mat2}{matrix; a symmetric logical matrix, rows and columns are
about tips of the tree, cells tell whether the difference between the
sampling dates of two isolates are below a certain threshold, e.g. 2 years.
Pairs which pass receive value 1. pairs that don't receive value 0.}

\item{MRCA_mat2}{matrix; a symmetric temporal matrix, rows and columns are
about tips of the tree, cells contain the cophenetic distance of tips, in
years.}

\item{geo_mat_ref}{matrix; symmetric logical matrix, rows and columns are
about tips of a tree, cells indicate whether pair is a geographic reference
pair.}

\item{time_mat2_ref}{matrix; symmetric logical matrix, rows and columns are
about tips of a tree, cells indicate whether the difference between the
sampling dates of two isolates are below a certain threshold, e.g. 2 years.
Pairs which pass receive value 1. pairs that don't receive value 0.}

\item{MRCA_mat2_ref}{matrix; a symmetric temporal matrix, rows and columns are
about tips of the tree, cells contain the cophenetic distance of tips, in
years.}

\item{Pmax}{numeric; vector of max temporal windows}

\item{Pmin}{numeric; vector of min temporal windows}
}
\description{
Calculate relative risks for a geographic matrix
}
\note{
within the pipeline application it seems \code{time_mat2} and
\code{time_mat2_ref} are the same, and similarly, \code{MRCA_mat2} and
\code{MRCA_mat2_ref} are the same. TODO: check and simplify if possible.
}
