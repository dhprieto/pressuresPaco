## ----setup, echo = TRUE, message = FALSE-----------------------------------
library("sensobol")
library("tidyverse")
library("data.table")
library("ggplot2")

theme_AP <- function() {
  theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA),
          strip.background = element_rect(fill = "white"), 
          legend.position = "top")
}

## ----binned_mean, cache=TRUE, echo = TRUE, fig.height=2, fig.cap = "Scatterplot of $y$ against $x_i$, $i=1,2,3$. The red dots show the mean $y$ value in each bin (we have set the number of bins arbitrarely at 30), and $N=2^{10}$. The model is the polynomial function in \\cite{Becker2014}, where $y=3 x_1 ^ 2 + 2 x_1 x_2 - 2 x_3$, $x_i \\sim \\mathcal{U}(0,1)$."----

poli <- function(X1, X2, X3) 3 * X1^2 + 2 * X1 * X2 - 2 * X3
poli_fun <- function(X) return(mapply(poli, X[, 1], X[, 2], X[, 3]))

set.seed(2)
N <- 2^10
params <- paste("$x_", 1:3, "$", sep = "")
mat <- sobol_matrices(N = N, params = params)
Y <- poli_fun(mat)

data.table(cbind(mat, Y)) %>%
  .[1:N] %>%
  melt(., measure.vars = params) %>%
  ggplot(., aes(value, Y)) + geom_point(size = 0.2) + stat_summary_bin(fun = "mean",
  geom = "point", colour = "red", size = 0.7) + facet_wrap(~variable) + labs(x = "$x$",
  y = "$y$") + theme_AP()


## ----ishi_plot, cache=TRUE, dependson="binned_mean", echo = TRUE, fig.height=2, fig.cap = "Scatterplot of $y$ against $x_i$, $i=1,2,3$. The red dots show the mean $y$ value in each bin (we have set the number of bins arbitrarely at 30), and $N=2^{10}$. The model is the \\cite{Ishigami1990} function."----

# PLOT ISHIGAMI FUNCTION TO ILLUSTRATE Ti ----------------------------------------

Y <- ishigami_Fun(mat)

data.table(cbind(mat, Y)) %>%
  .[1:N] %>%
  melt(., measure.vars = params) %>%
  ggplot(., aes(value, Y)) + geom_point(size = 0.2) + stat_summary_bin(fun = "mean",
  geom = "point", colour = "red", size = 0.7) + facet_wrap(~variable) + labs(x = "$x$",
  y = "$y$") + theme_AP()


## ----visualization_matrices, echo = TRUE, dependson=c("theme", "settings"), fig=TRUE, fig.height=2, cache=TRUE, fig.cap = "Sampling methods. Each dot is a sampling point. $N=2^{10}$."----

N <- 2^10
params <- paste("X", 1:3, sep = "")
type <- c("QRN", "LHS", "R")
set.seed(2)
prove <- lapply(type, function(type) sobol_matrices(N = N, params = params, type = type))
names(prove) <- type
lapply(prove, data.table) %>%
  lapply(., function(x) x[1:N]) %>%
  rbindlist(., idcol = "Method") %>%
  ggplot(., aes(X1, X2)) + geom_point(size = 0.2) + facet_wrap(~Method) + labs(x = "$x_1$",
  y = "$x_2$") + theme_AP()


## ----settings_sobolg, cache=TRUE-------------------------------------------
N <- 2^10
k <- 8
params <- paste("$x_", 1:k, "$", sep = "")
R <- 10^3
type <- "norm"
conf <- 0.95


## ----matrix_sobolg, cache=TRUE, dependson='settings_sobolg'----------------
set.seed(2)
mat <- sobol_matrices(N = N, params = params)


## ----model_sobolg, cache=TRUE, dependson='matrix_sobolg'-------------------
y <- sobol_Fun(mat)


## ----unc_sobolg, cache=TRUE, dependson="model_sobolg", fig.height=2, fig.cap="Empirical distribution of the Sobol' G model output.", message = FALSE----
plot_uncertainty(Y = y, N = N) + labs(y = "Counts", x = "$y$")


## ----scatter_sobolg, cache=TRUE, dependson=c("model_sobolg", "matrix_sobolg", "settings_sobolg"), fig.height=5, fig.cap="Scatter plots of model inputs against the model output for the Sobol' G function."----
plot_scatter(data = mat, N = N, Y = y, params = params)


## ----multiscatter_sobolg, cache=TRUE, dependson=c("model_sobolg", "matrix_sobolg", "settings_sobolg"), fig.height=5, fig.cap="Scatter plot matrix of pairs of model inputs for the Sobol' G function. The topmost and bottommost label facets refer to the $x$ and the $y$ axis respectively."----

plot_multiscatter(data = mat, N = N, Y = y, params = paste("$x_", 1:4, "$", sep = ""))


## ----indices_sobolg, cache=TRUE, dependson='model_sobolg'------------------

set.seed(2)
ind <- sobol_indices(Y = y, N = N, params = params, boot = TRUE, R = R, type = type,
  conf = conf)


## ----print_sobolg, cache=TRUE, dependson='indices_sobolg'------------------

cols <- colnames(ind$results)[1:5]
ind$results[, `:=`((cols), round(.SD, 3)), .SDcols = (cols)]
ind


## ----dummy_sobolg, cache=TRUE, dependson='indices_sobolg'------------------
set.seed(2)
ind.dummy <- sobol_dummy(Y = y, N = N, params = params, boot = TRUE, R = R)


## ----plot_indices_sobolg, cache=TRUE, dependson="indices_sobolg", fig.height=2.5, fig.cap="Sobol' indices of the Sobol' G function."----
plot(ind, dummy = ind.dummy)


## ----dynamics_population, cache=TRUE, echo = TRUE, fig.cap="Dynamics of the logistic population growth model for $N_0=3$, $r = 0.6$ and $K=100$.", fig.height=2.5----
X <- 3
r <- 0.6
K <- 100
for (i in 1:20) X[i + 1] <- X[i] + r * X[i] * (1 - X[i]/K)

dt <- data.table(X)[, `:=`(t, 0:20)]

ggplot(dt, aes(t, X)) + geom_line() + labs(x = "$t$", y = "$N$") + theme_bw() + theme(panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent",
    color = NA), legend.key = element_rect(fill = "transparent", color = NA))


## ----settings_population, cache=TRUE---------------------------------------
N <- 2^13
params <- c("$r$", "$K$", "$N_0$")
matrices <- c("A", "B", "AB", "BA")
first <- total <- "azzini"
order <- "second"
R <- 10^3
type <- "percent"
conf <- 0.95


## ----model_population, cache=TRUE------------------------------------------
population_growth <- function(r, K, X0) {
  X <- X0
  for (i in 0:20) X <- X + r * X * (1 - X/K)
  return(X)
}


## ----model_population_mapply, cache=TRUE, dependson='model:population'-----
population_growth_run <- function(dt) {
  return(mapply(population_growth, dt[, 1], dt[, 2], dt[, 3]))
}


## ----matrix_population, cache=TRUE, dependson='settings_population'--------
set.seed(2)
mat <- sobol_matrices(matrices = matrices, N = N, params = params, order = order,
  type = "LHS")


## ----transform_matrix_population, cache=TRUE, dependson="matrix_population"----
mat[, "$r$"] <- qnorm(mat[, "$r$"], 1.7, 0.3)
mat[, "$K$"] <- qnorm(mat[, "$K$"], 40, 1)
mat[, "$N_0$"] <- qunif(mat[, "$N_0$"], 10, 50)


## ----run_model_population, cache=TRUE, dependson=c("transform_matrix_population", "model_population_mapply", "model_population")----
y <- population_growth_run(mat)


## ----unc_population, cache=TRUE, dependson="run_model_population", fig.height=2, fig.cap="Empirical distribution of the logistic population growth model output.", message=FALSE----
plot_uncertainty(Y = y, N = N) + labs(y = "Counts", x = "$y$")


## ----quantiles, cache=TRUE, dependson='run_model_population'---------------
quantile(y, probs = c(0.01, 0.025, 0.5, 0.975, 0.99, 1))


## ----scatter_population, cache=TRUE, dependson=c("model_population", "matrix_population", "settings_population"), fig.height=3, fig.cap="Hexbin plot of model inputs against the model output for the population growth model."----
plot_scatter(data = mat, N = N, Y = y, params = params, method = "bin")


## ----multiscatter_population, cache=TRUE, dependson=c("model_population", "matrix_population", "settings_population"), fig.height=3, fig.cap="Scatterplot matrix of pairs of model inputs for the population growth model."----
plot_multiscatter(data = mat, N = N, Y = y, params = params, smpl = 2^11)


## ----indices_population, cache=TRUE, dependson=c("run_model_population", "settings_population")----
set.seed(2)
ind <- sobol_indices(matrices = matrices, Y = y, N = N, params = params, first = first,
  total = total, order = order, boot = TRUE, R = R, parallel = "no", type = type,
  conf = conf)


## ----print_population, cache=TRUE, dependson='indices_population'----------
cols <- colnames(ind$results)[1:5]
ind$results[, (cols):= round(.SD, 3), .SDcols = (cols)]
ind


## ----dummy_population, cache=TRUE, dependson=c("run_model_population", "indices_population")----
set.seed(2)
ind.dummy <- sobol_dummy(Y = y, N = N, params = params, boot = TRUE, R = R)


## ----plot_indices_population, cache=TRUE, dependson="indices_population", fig.height=3, fig.cap = "First and total-order Sobol' indices of the population growth model."----
plot(ind, dummy = ind.dummy)


## ----plot_indices2_population, cache=TRUE, dependson="indices_population", fig.height=2, fig.cap = "Second-order Sobol' indices."----
plot(ind, order = "second")


## ----settings_budworm, cache=TRUE------------------------------------------
N <- 2^9
params <- c("r_b", "K", "beta", "alpha", "r_s", "K_s", "K_e", "r_e", "P", "T_e")
order <- "first"
R <- 10^3
type <- "norm"
conf <- 0.95
times <- seq(0, 150, 1)
timeOutput <- seq(25, 150, 25)


## ----model_budworm, cache=TRUE---------------------------------------------
budworm_fun <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dB <- r_b * B * (1 - B/(K * S) * (T_e^2 + E^2)/E^2) - beta * B^2/((alpha^S)^2 +
      B^2)
    dS <- r_s * S * (1 - (S * K_e)/(E * K_s))
    dE <- r_e * E * (1 - E/K_e) - P * (B/S) * E^2/(T_e^2 + E^2)
    list(c(dB, dS, dE))
  })
}


## ----matrix_budworm, cache=TRUE, dependson='settings_budworm'--------------
set.seed(2)
mat <- sobol_matrices(N = N, params = params, order = order)


## ----transform_matrix_budworm, cache=TRUE, dependson='matrix_budworm'------
mat[, "r_b"] <- qunif(mat[, "r_b"], 1.52, 1.6)
mat[, "K"] <- qunif(mat[, "K"], 100, 355)
mat[, "beta"] <- qunif(mat[, "beta"], 20000, 43200)
mat[, "alpha"] <- qunif(mat[, "alpha"], 1, 2)
mat[, "r_s"] <- qunif(mat[, "r_s"], 0.095, 0.15)
mat[, "K_s"] <- qunif(mat[, "K_s"], 24000, 25440)
mat[, "K_e"] <- qunif(mat[, "K_e"], 1, 1.2)
mat[, "r_e"] <- qunif(mat[, "r_e"], 0.92, 1)
mat[, "P"] <- qunif(mat[, "P"], 0.0015, 0.00195)
mat[, "T_e"] <- qunif(mat[, "T_e"], 0.7, 0.9)


## ----plot_dynamics_budworm, cache=TRUE, echo = TRUE, fig=TRUE, dependson="transform_matrix_budworm", fig.height=2, fig.cap = "Dynamics of the spruce budworm and forest model. The vertical, dashed lines mark the times at which we will conduct the sensitivity analysis. Initial state values: $B=1,S=0.07,E=1$. The parameter values are the mean values of the distributions shown in Table 6."----
y.diff <- data.table(deSolve::ode(y = c(B = 0.1, S = 0.07, E = 1), times = seq(0,
  200, 1), func = budworm_fun, parms = colMeans(mat)))

melt(y.diff, measure.vars = c("B", "S", "E")) %>%
  ggplot(., aes(time, value)) + geom_line(size = 1) + geom_vline(xintercept = timeOutput,
  lty = 2) + labs(x = expression(italic(t)), y = "Value") + facet_wrap(~variable,
  scales = "free_y") + theme_AP()


## ----parallel_budworm, cache=TRUE------------------------------------------
library("foreach")
library("parallel")
library("doParallel")


## ----run_model_budworm, cache=TRUE, dependson=c("transform_matrix_budworm", "model_budworm"), message = FALSE----
n.cores <- makeCluster(floor(detectCores() * 0.75))
registerDoParallel(n.cores)
y <- foreach(i = 1:nrow(mat), .combine = "rbind", .packages = "sensobol") %dopar%
  {
    sobol_ode(d = mat[i, ], times = times, timeOutput = timeOutput, state = c(B = 0.1,
      S = 7, E = 1), func = budworm_fun)
  }
stopCluster(n.cores)


## ----arrange_output_budworm, cache=TRUE, dependson='run_model_budworm'-----
full.dt <- data.table(y)
print(full.dt)


## ----melt_budworm, cache=TRUE, dependson='arrange_output_budworm'----------
indices.dt <- melt(full.dt, measure.vars = c("B", "S", "E"))
print(indices.dt)


## ----new_params_budworm, cache=TRUE, echo = TRUE----------------------------
params <- c("$r_B$", "$K$", "$\\beta$", "$\\alpha$", "$r_S$", "$K_S$", "$K_E$", "$r_E$",
  "$P$", "$T_E$")


## ----sobol_budworm, cache=TRUE, dependson=c('melt_budworm',
## 'new_params_budworm')----
ncpus <- floor(detectCores() * 0.75)
set.seed(2)
indices <- indices.dt[, sobol_indices(Y = value, N = N, params = params, order = order,
  boot = TRUE, first = "jansen", R = R, parallel = "multicore", ncpus = ncpus)$results,
  .(variable, time)]


## ----dummy_sensobol, cache=TRUE, dependson='sobol_budworm'-----------------
indices.dummy <- indices.dt[, sobol_dummy(Y = value, N = N, params = params), .(variable,
  time)]


## ----plot_sobol_budworm_t, cache=TRUE, dependson='sobol_budworm',fig.height=7, fig.cap='Evolution of Sobol' indices through time in the spruce budworm and forest model. The dashed, horizontal blue line shows the $T_i$ of the dummy parameter.'----
ggplot(indices, aes(time, original, fill = sensitivity, color = sensitivity, group = sensitivity)) +
  geom_line() +
  geom_ribbon(aes(ymin = indices[sensitivity %in% c("Si", "Ti")]$low.ci,
    ymax = indices[sensitivity %in% c("Si", "Ti")]$high.ci,
    color = sensitivity), alpha = 0.1, linetype = 0) +
  geom_hline(data = indices.dummy[, parameters:= NULL][sensitivity == "Ti"],
    aes(yintercept = original, color = sensitivity, group = time),
    lty = 2, size = 0.1) +
  guides(linetype = FALSE, color = FALSE) +
  facet_grid(parameters ~ variable) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  labs(x = expression(italic(t)), y = "Sobol' indices") +
  theme_AP() + theme(legend.position = "top")


## ----benchmark_packages, cache=TRUE----------------------------------------
library("microbenchmark")


## ----benchmark_settings, cache=TRUE----------------------------------------
N <- 2^11
parameters <- c("N", "k")
R <- 10^2


## ----benchmark_matrix, cache=TRUE, dependson='benchmark_settings'----------
set.seed(2)
dt <- sobol_matrices(matrices = "A", N = N, params = parameters)
dt[, 1] <- floor(qunif(dt[, 1], 10, 10^2 + 1))
dt[, 2] <- floor(qunif(dt[, 2], 3, 100))


## ----benchmark_model, cache=TRUE, dependson='benchmark_matrix'-------------
n.cores <- makeCluster(floor(detectCores() * 0.75))
registerDoParallel(n.cores)
set.seed(2)
y <- foreach(i = 1:nrow(dt), .packages = c("sensobol", "sensitivity")
  ) %dopar% {
    params <- paste("x", 1:dt[i, "k"], sep = "")
    N <- dt[i, "N"]
    out <- microbenchmark::microbenchmark(
       "sensobol" = {
         params <- paste("X", 1:length(params), sep = "")
         mat <- sensobol::sobol_matrices(N = N, params = params, type = "R")
         y <- sensobol::metafunction(mat)
         ind <- sensobol::sobol_indices(Y = y, N = N, params = params,
           first = "jansen", total = "jansen", boot = TRUE, R = R)$results},
      "sensitivity" = {
        X1 <- data.frame(matrix(runif(length(params) * N), nrow = N))
        X2 <- data.frame(matrix(runif(length(params) * N), nrow = N))
        x <- sensitivity::soboljansen(model = sensobol::metafunction, 
          X1, X2, nboot = R)},
      times = 1)
  }
stopCluster(n.cores)


## ----benchmark_arrange, cache=TRUE, dependson='benchmark_model'------------
out <- rbindlist(y)[, time := time / 1e+06]


## ----plot_benchmark, cache=TRUE, dependson="benchmark", echo = TRUE, fig=TRUE, fig.height=2, fig.cap="Benchmark of the sensitivity and sensobol packages. The comparison has been done with the Jansen estimators."----
ggplot(out, aes(time, expr)) + geom_violin() + labs(x = "Time (Milliseconds)", y = "") +
  theme_AP()


## ----benchmark_time, cache=TRUE, dependson='benchmark_arrange'-------------
out[, median(time), expr]


## ----vars_settings, cache=TRUE---------------------------------------------
star.centers <- 100
h <- 0.1
params <- paste("X", 1:8, sep = "")


## ----vars_matrix, cache=TRUE, dependson = 'vars_settings'------------------
set.seed(2)
mat <- vars_matrices(star.centers = star.centers, h = h, params = params)


## ----vars_model, cache=TRUE, dependson='vars_matrix'-----------------------
y <- sobol_Fun(mat)


## ----vars_to, cache=TRUE, dependson='vars_model'---------------------------
ind <- vars_to(Y = y, star.centers = star.centers, params = params, h = h)
ind

