#' fitHeavyTail: Mean and Covariance Matrix Estimation under Heavy Tails
#'
#' Robust estimation methods for the mean vector, scatter matrix,
#' and covariance matrix (if it exists) from data (possibly containing NAs)
#' under multivariate heavy-tailed distributions such as angular Gaussian
#' (via Tyler's method), Cauchy, and Student's t distributions. Additionally,
#' a factor model structure can be specified for the covariance matrix.
#'
#'
#' @section Functions:
#' \code{\link{fit_Tyler}}, \code{\link{fit_Cauchy}}, and \code{\link{fit_mvt}}
#'
#' @section Help:
#' For a quick help see the README file:
#' \href{https://github.com/dppalomar/fitHeavyTail/blob/master/README.md}{GitHub-README}.
#'
#' For more details see the vignette:
#' \href{https://CRAN.R-project.org/package=fitHeavyTail/vignettes/CovarianceEstimationHeavyTail.html}{CRAN-vignette}.
#'
#' @author Daniel P. Palomar and Rui Zhou
#'
#' @references
#' Ying Sun, Prabhu Babu, and Daniel P. Palomar, "Regularized Tyler's Scatter Estimator: Existence, Uniqueness, and Algorithms,"
#' IEEE Trans. on Signal Processing, vol. 62, no. 19, pp. 5143-5156, Oct. 2014. <https://doi.org/10.1109/TSP.2014.2348944>
#'
#' Ying Sun, Prabhu Babu, and Daniel P. Palomar, "Regularized Robust Estimation of Mean and Covariance Matrix Under Heavy-Tailed Distributions,"
#' IEEE Trans. on Signal Processing, vol. 63, no. 12, pp. 3096-3109, June 2015. <https://doi.org/10.1109/TSP.2015.2417513>
#'
#' Chuanhai Liu and Donald B. Rubin, "ML estimation of the t-distribution using EM and its extensions, ECM and ECME,"
#' Statistica Sinica (5), pp. 19-39, 1995.
#'
#' Rui Zhou, Junyan Liu, Sandeep Kumar, and Daniel P. Palomar, "Robust factor analysis parameter estimation,"
#' Lecture Notes in Computer Science (LNCS), 2019. <https://arxiv.org/abs/1909.12530>
#'
#' @docType package
#' @name fitHeavyTail-package
NULL
