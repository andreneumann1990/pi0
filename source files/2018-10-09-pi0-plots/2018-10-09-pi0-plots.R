library(ggplot2) #ggplot


main_2 <- function() {
  # default setting:
  # rho = 0.5, n = 100, m = 100, pi0 = 0.5, mu = 1, lambda = 0.5

  # load data
  load("2018-10-08-pi0-simulation.RData")

  for (i in 1:length(results))
    plotData(results[[i]])
}


# functions ####

theme_set(theme_bw(base_size = 20))
theme_update(plot.title = element_text(hjust = 0.5))
update_geom_defaults("point", list(size = 2))
update_geom_defaults("line", list(size = 1.2))

plotData <- function(data) {
  var <- colnames(data)[2]
  data <- as.data.frame(data)
  for (i in 2:4)
    data[, i] <- as.numeric(as.character(data[, i])) # convert factor to num
  ylab <- "bias"
  gg <- ggplot(data, aes_string(x = var, y = ylab, group = "type", colour = "type")) +
    labs(x = parse(text = var), y = ylab, colour = "") +
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    geom_line()
  ggsave(paste("2018-10-09-pi0-plots-", var, "-", ylab, ".png", sep = ""), width = 6.4, height = 6.4)
  print(gg)
  ylab <- "RMSE"
  gg <- ggplot(data, aes_string(x = var, y = ylab, group = "type", colour = "type")) +
    labs(x = parse(text = var), y = ylab, colour = "") +
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    geom_line()
  ggsave(paste("2018-10-09-pi0-plots-", var, "-", ylab, ".png", sep = ""), width = 6.4, height = 6.4)
  print(gg)
}


####