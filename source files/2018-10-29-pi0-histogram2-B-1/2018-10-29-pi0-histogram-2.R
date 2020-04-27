library("mvtnorm")
library("copula")

library("doParallel")
library("doRNG")

library("ggplot2")

theme_set(theme_bw(base_size = 20))
theme_update(plot.title = element_text(hjust = 0.5))
update_geom_defaults("point", list(size = 2))
update_geom_defaults("line", list(size = 1.2))


# Plots the histogram of SS and BSS over L tests, respectively, and calculates
# the resulting mean bias and root mean squared error.
main_4 <- function() {
  seed <- seed(820177) #820177

  # default setting
  default <- list(rho = 0.5, n = 100, m = 100, pi0 = 0.5, mu = 1, lambda = 0.5)
  L <- 1000 # number of test per setting
  B <- 1 # bootstrap repetition

  registerDoParallel(cores = max(detectCores() - 1, 1))
  set.seed(seed)
  t <- Sys.time()
  results <- estimate_pi0(default$rho, default$n, default$m, default$pi0,
                          default$mu, default$lambda, L, B)
  attr(results, "rng") <- NULL
  print(Sys.time() - t)
  closeAllConnections()
  pi0_SS <- c(bias = mean(results[, "SS"] - default$pi0),
              RMSE = sqrt(mean((results[, "SS"] - default$pi0)^2)))
  pi0_BSS <- c(bias = mean(results[, "BSS"] - default$pi0),
               RMSE = sqrt(mean((results[, "BSS"] - default$pi0)^2)))
  plot_histograms(results, default$pi0)

  results <- list(pi0_hat = results, pi0 = default$pi0, pi0_hat_SS = pi0_SS,
                  pi0_hat_BSS = pi0_BSS)
  save(results, file = "2018-10-29-pi0-histogram-2.RData")
  results
}


# functions ####

calculate_Pi0 <- function(p, lambda) {
  (1 - ecdf(p)(lambda)) / (1 - lambda)
}

estimate_pi0 <- function(rho_, n_, m_, pi0_, mu_, lambda_, L, B) {
  m0_ <- round(m_ * pi0_)
  mu <- c(numeric(m0_), rep(mu_, m_ - m0_))
  C_ <- normalCopula(rho_, dim = m_)
  kn_ <- ceiling(sqrt(n_))
  foreach(l = 1:L, .export = lsf.str(.GlobalEnv),
          .packages = c("mvtnorm", "copula"), .combine = rbind) %dorng% {
    U <- rCopula(n_, C_)
    X <- matrix(nrow = n_, ncol = m_)
    for (j in 1:m_)
      X[, j] <- qnorm(U[, j], mu[j], sd = 1)

    mu_hat <- apply(X, 2, mean)
    sd_hat <- apply(X, 2, sd)

    BSS <- numeric(B)
    X_star <- matrix(nrow = kn_, ncol = m_)
    for (b in 1:B) {
      for (j in 1:m_)
        X_star[, j] <- rnorm(kn_, mu_hat[j], sd_hat[j])
      T_star <- sqrt(kn_) * apply(X_star, 2, mean)/apply(X_star, 2, sd)
      p <- p_values(T_star, kn_)
      BSS[b] <- calculate_Pi0(p, lambda_)
    }
    p <- p_values(sqrt(n_) * mu_hat/sd_hat, n_)
    c(SS = calculate_Pi0(p, lambda_), BSS = mean(BSS))
  }
}

p_values <- function(T, n) {
  2 * (1 - pt(abs(T), df = n - 1))
}

plot_histograms <- function(results, pi0) {
  gg <- ggplot() +
    ggtitle("Histogram of Schweder-Spjotvoll \n estimates") +
    labs(x = expression(hat(pi)[0] - pi[0]), colour = "") +
    coord_cartesian(xlim = c(-0.5, 0.5)) + # zoom
    geom_histogram(aes(results[, "SS"] - pi0, stat(density)), binwidth = 0.01)
  ggsave("2018-10-29-pi0-histogram-2-SS.png", width = 6.4, height = 6.4)
  print(gg)

  gg <- ggplot() +
    ggtitle("Histogram of \n Bootstrap-Schweder-Spjotvoll \n estimates") +
    labs(x = expression(hat(pi)[0] - pi[0]), colour = "") +
    coord_cartesian(xlim = c(-0.5, 0.5)) + # zoom
    geom_histogram(aes(results[, "BSS"] - pi0, stat(density)), binwidth = 0.01)
  ggsave("2018-10-29-pi0-histogram-2-BSS.png", width = 6.4, height = 6.4)
  print(gg)
}

seed <- function(seed = NULL) {
  # Prints a given random seed or creates a new one.
  if (is.null(seed))
    seed <- as.numeric(paste(floor(runif(6, 0, 10)), collapse = ""))
  print(list(seed = seed))
  seed
}

#