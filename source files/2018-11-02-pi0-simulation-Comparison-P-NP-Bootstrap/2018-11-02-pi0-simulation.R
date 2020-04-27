library("mvtnorm")
library("copula")

library("doParallel")
library("doRNG")


# Simulation of standard (parametric) and non-parametric BSS.
main <- function() {
  seed <- seed(968450) #968450

  # default setting
  default <- list(rho = 0.5, n = 100, m = 100, pi0 = 0.5, mu = 1, lambda = 0.5)

  rho <- 0:19/20 # equi-correlation
  n <- c(10*1:10, 100*2:10) # sample size
  m <- c(10*1:10, 100*2:10) # number of null hypotheses
  pi0 <- 0:20/20 # true proportion of true null hypotheses
  mu <- 1:19/10 # mean / effect size
  lambda <- 1:19/20 # tuning parameter
  L <- 500 # number of test per setting
  B <- 1000 # bootstrap repetition

  compare_biasAndRMSE(rho, n, m, pi0, mu, lambda, default, L, B, seed)
}


# functions ####

calculate_Pi0 <- function(p, lambda) {
  (1 - ecdf(p)(lambda)) / (1 - lambda)
}

compare_biasAndRMSE <- function(rho, n, m, pi0, mu, lambda, default, L, B, seed) {
  registerDoParallel(cores = max(detectCores() - 1, 1))
  # load("2018-11-02-pi0-simulation.RData")
  results <- list()
  # results for equi-correlation rho ####
  set.seed(seed)
  l <- length(rho)
  results$rho <- matrix(nrow = 2*l, ncol = 4)
  colnames(results$rho) <- c("type", "rho", "bias", "RMSE")
  results$rho[1:l, 1] <- "BSS"
  results$rho[1:l, 2] <- rho
  results$rho[(l + 1):(2*l), 1] <- "BSS_nonparametric"
  results$rho[(l + 1):(2*l), 2] <- rho
  i <- 1
  for (rho_ in rho) {
    print(paste("rho = ", rho_, sep = ""))
    t <- Sys.time()
    pi0_hat <- estimate_pi0(rho_, default$n, default$m, default$pi0, default$mu,
                            default$lambda, L, B)
    results$rho[i, 3] <- mean(pi0_hat[, 1] - default$pi0)
    results$rho[i, 4] <- sqrt(mean((pi0_hat[, 1] - default$pi0)^2))
    results$rho[l + i, 3] <- mean(pi0_hat[, 2] - default$pi0)
    results$rho[l + i, 4] <- sqrt(mean((pi0_hat[, 2] - default$pi0)^2))
    i <- i + 1
    print(Sys.time() - t)
  }
  save(results,
       file = paste(Sys.Date(), "-pi0-simulation.RData", sep = ""))
  # results for sample size n ####
  set.seed(seed)
  l <- length(n)
  results$n <- matrix(nrow = 2*l, ncol = 4)
  colnames(results$n) <- c("type", "n", "bias", "RMSE")
  results$n[1:l, 1] <- "BSS"
  results$n[1:l, 2] <- n
  results$n[(l + 1):(2*l), 1] <- "BSS_nonparametric"
  results$n[(l + 1):(2*l), 2] <- n
  i <- 1
  for (n_ in n) {
    print(paste("n = ", n_, sep = ""))
    t <- Sys.time()
    pi0_hat <- estimate_pi0(default$rho, n_, default$m, default$pi0, default$mu,
                            default$lambda, L, B)
    results$n[i, 3] <- mean(pi0_hat[, 1] - default$pi0)
    results$n[i, 4] <- sqrt(mean((pi0_hat[, 1] - default$pi0)^2))
    results$n[l + i, 3] <- mean(pi0_hat[, 2] - default$pi0)
    results$n[l + i, 4] <- sqrt(mean((pi0_hat[, 2] - default$pi0)^2))
    i <- i + 1
    print(Sys.time() - t)
  }
  save(results,
       file = paste(Sys.Date(), "-pi0-simulation.RData", sep = ""))
  # results for number of null hypotheses m ####
  set.seed(seed)
  l <- length(m)
  results$m <- matrix(nrow = 2*l, ncol = 4)
  colnames(results$m) <- c("type", "m", "bias", "RMSE")
  results$m[1:l, 1] <- "BSS"
  results$m[1:l, 2] <- m
  results$m[(l + 1):(2*l), 1] <- "BSS_nonparametric"
  results$m[(l + 1):(2*l), 2] <- m
  i <- 1
  for (m_ in m) {
    print(paste("m = ", m_, sep = ""))
    t <- Sys.time()
    pi0_hat <- estimate_pi0(default$rho, default$n, m_, default$pi0, default$mu,
                            default$lambda, L, B)
    results$m[i, 3] <- mean(pi0_hat[, 1] - default$pi0)
    results$m[i, 4] <- sqrt(mean((pi0_hat[, 1] - default$pi0)^2))
    results$m[l + i, 3] <- mean(pi0_hat[, 2] - default$pi0)
    results$m[l + i, 4] <- sqrt(mean((pi0_hat[, 2] - default$pi0)^2))
    i <- i + 1
    print(Sys.time() - t)
  }
  save(results,
       file = paste(Sys.Date(), "-pi0-simulation.RData", sep = ""))
  # results for true proportion pi0 ####
  set.seed(seed)
  l <- length(pi0)
  results$pi0 <- matrix(nrow = 2*l, ncol = 4)
  colnames(results$pi0) <- c("type", "pi0", "bias", "RMSE")
  results$pi0[1:l, 1] <- "BSS"
  results$pi0[1:l, 2] <- pi0
  results$pi0[(l + 1):(2*l), 1] <- "BSS_nonparametric"
  results$pi0[(l + 1):(2*l), 2] <- pi0
  i <- 1
  for (pi0_ in pi0) {
    print(paste("pi0 = ", pi0_, sep = ""))
    t <- Sys.time()
    pi0_hat <- estimate_pi0(default$rho, default$n, default$m, pi0_, default$mu,
                            default$lambda, L, B)
    results$pi0[i, 3] <- mean(pi0_hat[, 1] - pi0_)
    results$pi0[i, 4] <- sqrt(mean((pi0_hat[, 1] - pi0_)^2))
    results$pi0[l + i, 3] <- mean(pi0_hat[, 2] - pi0_)
    results$pi0[l + i, 4] <- sqrt(mean((pi0_hat[, 2] - pi0_)^2))
    i <- i + 1
    print(Sys.time() - t)
  }
  save(results,
       file = paste(Sys.Date(), "-pi0-simulation.RData", sep = ""))
  # results for effect size mu ####
  set.seed(seed)
  l <- length(mu)
  results$mu <- matrix(nrow = 2*l, ncol = 4)
  colnames(results$mu) <- c("type", "mu", "bias", "RMSE")
  results$mu[1:l, 1] <- "BSS"
  results$mu[1:l, 2] <- mu
  results$mu[(l + 1):(2*l), 1] <- "BSS_nonparametric"
  results$mu[(l + 1):(2*l), 2] <- mu
  i <- 1
  for (mu_ in mu) {
    print(paste("mu = ", mu_, sep = ""))
    t <- Sys.time()
    pi0_hat <- estimate_pi0(default$rho, default$n, default$m, default$pi0, mu_,
                            default$lambda, L, B)
    results$mu[i, 3] <- mean(pi0_hat[, 1] - default$pi0)
    results$mu[i, 4] <- sqrt(mean((pi0_hat[, 1] - default$pi0)^2))
    results$mu[l + i, 3] <- mean(pi0_hat[, 2] - default$pi0)
    results$mu[l + i, 4] <- sqrt(mean((pi0_hat[, 2] - default$pi0)^2))
    i <- i + 1
    print(Sys.time() - t)
  }
  save(results,
       file = paste(Sys.Date(), "-pi0-simulation.RData", sep = ""))
  # results for tuning parameter lambda ####
  set.seed(seed)
  l <- length(lambda)
  results$lambda <- matrix(nrow = 2*l, ncol = 4)
  colnames(results$lambda) <- c("type", "lambda", "bias", "RMSE")
  results$lambda[1:l, 1] <- "BSS"
  results$lambda[1:l, 2] <- lambda
  results$lambda[(l + 1):(2*l), 1] <- "BSS_nonparametric"
  results$lambda[(l + 1):(2*l), 2] <- lambda
  i <- 1
  for (lambda_ in lambda) {
    print(paste("lambda = ", lambda_, sep = ""))
    t <- Sys.time()
    pi0_hat <- estimate_pi0(default$rho, default$n, default$m, default$pi0, default$mu,
                            lambda_, L, B)
    results$lambda[i, 3] <- mean(pi0_hat[, 1] - default$pi0)
    results$lambda[i, 4] <- sqrt(mean((pi0_hat[, 1] - default$pi0)^2))
    results$lambda[l + i, 3] <- mean(pi0_hat[, 2] - default$pi0)
    results$lambda[l + i, 4] <- sqrt(mean((pi0_hat[, 2] - default$pi0)^2))
    i <- i + 1
    print(Sys.time() - t)
  }
  save(results,
       file = paste(Sys.Date(), "-pi0-simulation.RData", sep = ""))
  #
  closeAllConnections()
  results
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

    BSS_par <- numeric(B)
    BSS_nonpar <- numeric(B)
    X_star <- matrix(nrow = kn_, ncol = m_)
    for (b in 1:B) {
      for (j in 1:m_)
        X_star[, j] <- rnorm(kn_, mu_hat[j], sd_hat[j])
      T_star <- sqrt(kn_) * apply(X_star, 2, mean)/apply(X_star, 2, sd)
      p <- p_values(T_star, kn_)
      BSS_par[b] <- calculate_Pi0(p, lambda_)
    }
    for (b in 1:B) {
      for (j in 1:m_)
        X_star[, j] <- X[sample(n_, kn_, replace = TRUE), j]
      T_star <- sqrt(kn_) * apply(X_star, 2, mean)/apply(X_star, 2, sd)
      p <- p_values(T_star, kn_)
      BSS_nonpar[b] <- calculate_Pi0(p, lambda_)
    }
    c(BSS = mean(BSS_par),
      BSS_np = mean(BSS_nonpar))
  }
}

p_values <- function(T, n) {
  2 * (1 - pt(abs(T), df = n - 1))
}

seed <- function(seed = NULL) {
  # Prints a given random seed or creates a new one.
  if (is.null(seed))
    seed <- as.numeric(paste(floor(runif(6, 0, 10)), collapse = ""))
  print(list(seed = seed))
  seed
}

#