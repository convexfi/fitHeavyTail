plotConvergence <- function(res_fit) {
  if (is.null(res_fit$iterations_record))
    stop("Fitting result does not contain iteration converge. Make sure to use \"return_iterates = TRUE\" when doing the fitting.")

  p_all <- list()
  data <- data.frame(iteration = 0:(length(res_fit$iterations_record)-1))
  if (!is.null(res_fit$iterations_record[[1]]$nu)) {
    data$nu <- sapply(res_fit$iterations_record, `[[`, "nu")
    p_nu <- ggplot(data, aes(x = iteration, y = nu)) +
      geom_line() + geom_point() +
      ggtitle("Convergence of nu")
    print(p_nu)
    p_all <- c(p_all, list(p_nu))
  }
  if (!is.null(res_fit$iterations_record[[1]]$mu)) {
    mu_matrix <- sapply(res_fit$iterations_record, `[[`, "mu")
    rownames(mu_matrix) <- paste0("mu", 1:nrow(mu_matrix))
    data <- cbind(data, as.data.frame(t(mu_matrix)))
    # ggplot(data, aes(x = iteration, y = mu2)) +
    #   geom_line() + geom_point()
    p_mu <- ggplot(reshape2::melt(data, measure.vars = rownames(mu_matrix)), aes(x = iteration, y = value, col = variable)) +
      geom_line() + geom_point() +
      ggtitle("Convergence of mu")
    print(p_mu)
    p_all <- c(p_all, list(p_mu))
  }

  if (!is.null(res_fit$iterations_record[[1]]$scale)) {
    diag_scale_matrix <- sapply(res_fit$iterations_record, function(x) diag(x$scale))
    rownames(diag_scale_matrix) <- paste0("diag_scale", 1:nrow(diag_scale_matrix))
    data <- cbind(data, as.data.frame(t(diag_scale_matrix)))
    p_scale <- ggplot(reshape2::melt(data, measure.vars = rownames(diag_scale_matrix)), aes(x = iteration, y = value, col = variable)) +
      geom_line() + geom_point() +
      ggtitle("Convergence of scale matrix")
    print(p_scale)
    p_all <- c(p_all, list(p_scale))
  }
  return(p_all)
}
