#' fitHeavyTail: Mean and Covariance Matrix Estimation under Heavy Tails
#'
#' Estimation of mean vector and covariance matrix for heavy-tailed
#' distributions. In particular, a multivariate Student t distribution is
#' fitted to the data (possibly containing NAs).
#' Special emphasis is placed on a fast implementation.
#' The package is based on the paper:
#' Rui ZHOU, Junyan Liu, Sandeep Kumar, and Daniel P. Palomar, "Robust factor
#' analysis parameter estimation," Lecture Notes in Computer Science, 2020.
#' <https://arxiv.org/abs/1909.12530>
#'
#' @section Functions:
#' \code{\link{momentsTyler}}, \code{\link{momentsCauchy}}, and \code{\link{momentsStudentt}}
#'
#' @section Help:
#' For a quick help see the README file:
#' \href{https://github.com/dppalomar/fitHeavyTail/blob/master/README.md}{GitHub-README}.
#'
#' @author Daniel P. Palomar and Rui ZHOU
#'
#' @docType package
#' @name fitHeavyTail-package
NULL
