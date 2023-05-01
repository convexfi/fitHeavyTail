#' fitHeavyTail: Mean and Covariance Matrix Estimation under Heavy Tails
#'
#' Robust estimation methods for the mean vector, scatter matrix,
#' and covariance matrix (if it exists) from data (possibly containing NAs)
#' under multivariate heavy-tailed distributions such as angular Gaussian
#' (via Tyler's method), Cauchy, and Student's t distributions. Additionally,
#' a factor model structure can be specified for the covariance matrix. The
#' latest revision also includes the multivariate skewed t distribution.
#' The package is based on the papers: Sun, Babu, and Palomar (2014);
#' Sun, Babu, and Palomar (2015); Liu and Rubin (1995);
#' Zhou, Liu, Kumar, and Palomar (2019); Pascal, Ollila, and Palomar (2021).
#'
#'
#' @section Functions:
#' \code{\link{fit_Tyler}}, \code{\link{fit_Cauchy}}, \code{\link{fit_mvt}}, and \code{\link{fit_mvst}}.
#'
#' @section Help:
#' For a quick help see the README file:
#' \href{https://github.com/convexfi/fitHeavyTail/blob/master/README.md}{GitHub-README}.
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
#' Chuanhai Liu, Donald B. Rubin, and Ying Nian Wu, "Parameter Expansion to Accelerate EM: The PX-EM Algorithm,"
#' Biometrika, Vol. 85, No. 4, pp. 755-770, Dec., 1998
#'
#' Rui Zhou, Junyan Liu, Sandeep Kumar, and Daniel P. Palomar, "Robust factor analysis parameter estimation,"
#' Lecture Notes in Computer Science (LNCS), 2019. <https://arxiv.org/abs/1909.12530>
#'
#' Esa Ollila, Daniel P. Palomar, and Frédéric Pascal, "Shrinking the Eigenvalues of M-estimators of Covariance Matrix,"
#' IEEE Trans. on Signal Processing, vol. 69, pp. 256-269, Jan. 2021. <https://doi.org/10.1109/TSP.2020.3043952>
#'
#' Frédéric Pascal, Esa Ollila, and Daniel P. Palomar, "Improved estimation of the degree of freedom parameter of
#' multivariate t-distribution," in Proc. European Signal Processing Conference (EUSIPCO), Dublin, Ireland, Aug. 23-27, 2021.
#' <https://doi.org/10.23919/EUSIPCO54536.2021.9616162>
#'
#' @docType package
#' @name fitHeavyTail-package
NULL
