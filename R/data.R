#' Synthetic heavy-tailed data.
#'
#' Synthetic heavy-tailed data generated with \code{mvtnorm::rmvt(n=T, delta=mu, sigma=(nv-2)/nv*R_cov, df=nv)}
#'
#' @format A list containing the following:
#' \describe{
#'   \item{N}{The dimension}
#'   \item{T}{The number of samples}
#'   \item{mu}{The mean vector}
#'   \item{cov}{The covariance matrix}
#'   \item{nv}{The degrees of freedom}
#'   \item{X}{The TxN data matrix}
#' }
"heavy_data"
