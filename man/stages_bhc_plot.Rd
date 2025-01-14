% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/11-custom-functions.R
\name{stages_bhc_plot}
\alias{stages_bhc_plot}
\title{Plotting of staged event tree for every step of the BHC algorithm}
\usage{
stages_bhc_plot(
  object,
  plot_var = NULL,
  score = function(x) {     return(-BIC(x)) },
  max_iter = Inf,
  scope = NULL,
  ignore = object$name_unobserved,
  trace = 0
)
}
\arguments{
\item{object}{an object of class \code{sevt} with fitted probabilities and
data, as returned by \code{full} or \code{sevt_fit}.}

\item{plot_var}{variable for plotting the barplots.}

\item{score}{the score function to be maximized.}

\item{max_iter}{the maximum number of iterations per variable.}

\item{scope}{names of variables that should be considered for the optimization.}

\item{ignore}{vector of stages which will be ignored and left untouched,
by default the name of the unobserved stages stored in
\code{object$name_unobserved}.}

\item{trace}{if >0 increasingly amount of info
is printed (via \code{message}).}
}
\value{
The final staged event tree obtained and print a series of plots.
}
\description{
Displaying the plot of a staged event tree together with the
probabilities respective to the \code{plot_var} variable.
}
\details{
The function simultaniously performs Backward Hill Climbing as
implemented in \code{stages_bhc()}. Additionally, at every merge of two stages the
current staged tree together with the probabilities corresponding to the
\code{plot_var} are plotted.
}
