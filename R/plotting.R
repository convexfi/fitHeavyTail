#
#
#
plotConvergence <- function(res_fit) {
  if (is.null(res_fit$iterates_record))
    stop("Fitting result does not contain iteration converge. Make sure to use \"return_iterates = TRUE\" when doing the fitting.")

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Please install package \"ggplot2\"", call. = FALSE)

  if (!requireNamespace("reshape2", quietly = TRUE))
    stop("Please install package \"reshape2\"", call. = FALSE)

  p_all <- list()
  data <- data.frame(iteration = 0:(length(res_fit$iterates_record)-1))
  if (!is.null(res_fit$iterates_record[[1]]$nu)) {
    data$nu <- sapply(res_fit$iterates_record, `[[`, "nu")
    data$nu_div_nu_2 <- fnu(data$nu)
    p_nu <- ggplot2::ggplot(data, ggplot2::aes(x = iteration, y = nu)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::ggtitle("Convergence of nu")
    print(p_nu)
    p_all <- c(p_all, list(p_nu))

    p_nu_div_nu_2 <- ggplot2::ggplot(data, ggplot2::aes(x = iteration, y = nu_div_nu_2)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::ggtitle("Convergence of nu/(nu-2)") + ggplot2::ylab(nu/(nu-2))
    print(p_nu_div_nu_2)
    p_all <- c(p_all, list(p_nu_div_nu_2))
  }
  if (!is.null(res_fit$iterates_record[[1]]$mu)) {
    mu_matrix <- sapply(res_fit$iterates_record, `[[`, "mu")
    rownames(mu_matrix) <- paste0("mu", 1:nrow(mu_matrix))
    data <- cbind(data, as.data.frame(t(mu_matrix)))
    # ggplot(data, aes(x = iteration, y = mu2)) +
    #   geom_line() + geom_point()
    p_mu <- ggplot2::ggplot(reshape2::melt(data, measure.vars = rownames(mu_matrix)), ggplot2::aes(x = iteration, y = value, col = variable)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::ggtitle("Convergence of mu")
    print(p_mu)
    p_all <- c(p_all, list(p_mu))
  }

  if (!is.null(res_fit$iterates_record[[1]]$scatter)) {
    diag_scatter_matrix <- sapply(res_fit$iterates_record, function(x) diag(x$scatter))
    rownames(diag_scatter_matrix) <- paste0("diag_scatter", 1:nrow(diag_scatter_matrix))
    data <- cbind(data, as.data.frame(t(diag_scatter_matrix)))
    p_scatter <- ggplot2::ggplot(reshape2::melt(data, measure.vars = rownames(diag_scatter_matrix)), ggplot2::aes(x = iteration, y = value, col = variable)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::ggtitle("Convergence of scatter matrix")
    print(p_scatter)
    p_all <- c(p_all, list(p_scatter))
  }
  return(invisible(p_all))
}
