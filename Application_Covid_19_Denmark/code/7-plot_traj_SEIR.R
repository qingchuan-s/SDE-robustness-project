setwd('~/Application_Covid_19_Denmark/estimates')


###############################################################################
# plot SEIR traj with fixed random seeds
# we selected 50 trajectories with smallest square mean error 
# the code for select the seeds is below

df_full <- read_excel("Data Denmark.xlsx")
n_pop <- 5.83e6

start_date <- 225
end_date <- 404

df_infected <-df_full[225:404,]

observations <- df_infected$approx_infectious_cases/n_pop
n_time <- length(observations)

N <-  n_pop  
i0 <- df_infected$approx_infectious_cases[1]/N
r0 <- sum(df_full$New_cases[1:(start_date-8)])/N
initials <- c(s= 1 - 1.5*i0 - r0, e= 0.5*i0, i = i0)

data1 <- data.frame(df_infected$Date_reported)
colnames(data1) <- c('Time')
data1$observations <- df_infected$approx_infectious_cases/N
data1$Time <- as.Date(data1$Time)



simulate_sir_sde <- function(N, h, par, start, W1, W2) {
  alpha = par[1]; gamma = par[2]; beta = par[3]
  sigma1.sq <- par[4]^2; sigma2.sq <- par[5]^2; sigma3.sq <- par[6]^2
  S <- numeric(N); E <- numeric(N); I <- numeric(N)
  S[1] <- start[1]; E[1] <- start[2]; I[1] <- start[3]
  for (i in 2:N) {
    S[i] <- S[i-1] + h * (-alpha * S[i-1] * I[i-1]) + sqrt(h * sigma1.sq) * W1[i]
    E[i] <- E[i-1] + h * (alpha * S[i-1] * I[i-1] - gamma * E[i-1]) + sqrt(h * sigma2.sq) * W2[i]
    I[i] <- I[i-1] + h * (gamma * E[i-1] - beta * I[i-1]) + sqrt(h * sigma3.sq) * W2[i]
    if(I[i] <= 0)I[i] = 0
  }
  return(matrix(c(S, E, I), ncol = 3))
}

########## ode
SEIR_fun <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    ds    <- -alpha  * s * i
    de <- alpha  * s * i- gamma * e
    di <- gamma * e - beta * i
    
    return(list(c(ds, de, di)))
  })
}

pars.ode <- c(9.984568e-01, 8.952011e-01, 8.709919e-01, 4.669617e-05)
pars <- list(alpha = pars.ode[1], gamma = pars.ode[2], beta = pars.ode[3])
yini  <- c(s= 1-1.5*pars.ode[4]-r0, e = 0.5 * pars.ode[4], i = pars.ode[4])
times <- seq(0, length(observations)-1, by = 0.5)
out   <- ode(yini, times, SEIR_fun, pars)
df_ode_seir <- data.frame(out)


df_ode_seir_matched <- df_ode_seir %>% filter(time %in% seq(0, length(observations)-1, by = 1))
df_ode_seir_matched$time <- data1$Time
df_ode_seir_matched$legend_label <- "Deterministic SEIR model"

#==============================================================================#
best_seeds <- read.csv('best_seeds SEIR.csv')[,1]

pars_opt <- c(9.736867e-01, 5.006845e-01, 8.390785e-01, 1.478588e-05, 1e-5, 4.817349e-04)

param_list <- list(
  list(par = pars_opt, label = "Optimal a")
)


sim_all <- lapply(param_list, function(pinfo) {
  sim_dfs <- lapply(seq_along(best_seeds), function(i) {
    s <- best_seeds[i]
    set.seed(s); W1 <- rnorm(n_time); W2 <- rnorm(n_time)
    I_sim <- simulate_sir_sde(n_time, 1, pinfo$par, initials, W1, W2)[,3]
    data.frame(
      Time   = data1$Time,
      I      = I_sim,
      seed   = i,
      label  = pinfo$label
    )
  })
  bind_rows(sim_dfs)
}) %>% bind_rows()


summary_sim <- sim_all %>%
  group_by(label, Time) %>%
  summarise(
    min_I    = min(I),
    max_I    = max(I),
    median_I = median(I),
    .groups  = "drop"
  )


summary_sim$legend_label <- "Stochastic SEIR model - range"
data1$legend_label <- "Observations"
df_ode_seir_matched$legend_label <- "Deterministic SEIR model"
summary_sim$legend_label_median <- "Stochastic SEIR model - median"

p2 <- ggplot() +
  geom_ribbon(data = summary_sim, aes(x = Time, ymin = min_I , ymax = max_I, fill = legend_label), alpha = 0.3) +
  geom_line(data = summary_sim, aes(x = Time, y = median_I , color = legend_label_median), size = 1) +
  geom_line(data = data1, aes(x = Time, y = observations, color = legend_label), size = 1) +
  geom_line(data = df_ode_seir_matched, aes(x = time, y = i, color = legend_label), size = 1) +
  
  scale_color_manual(
    name = "Optimal parameter estimates",
    values = c(
      "Observations" = "black",
      "Deterministic SEIR model" = "blue",
      "Stochastic SEIR model - median" = "red"
    )
  ) +
  scale_fill_manual(
    name = "",
    values = c("Stochastic SEIR model - range" = "darkgray")
  ) +
  labs(
    y = "(Proportion of) Infectious"
  ) +
  scale_x_date(date_labels = "%Y-%m-%d") +
  theme_bw() +
  theme(
    axis.text.y=element_blank(), axis.title.y=element_blank(),
    axis.text.x=element_blank(), axis.title.x=element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    legend.position = c(0.20, 0.85),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.background = element_rect(
      fill = "transparent", size = 0.5, linetype = "solid", colour = "transparent"
    ),
    strip.text = element_text(size = 15)
  )+ guides(
    color = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )+ylim(c(0, 0.0075))
p2<- ggdraw(p2) +
  draw_label("B", x = 0.025, y = 0.97, fontface = 'bold', size = 30, hjust = 0, vjust = 1)


################################################################################
# ============================================================================ #
################################################################################
# The code to generate seeds of the trajectories

seeds_all <- sample.int(1e8, 1e4)

pars_opt <- c(9.736867e-01, 5.006845e-01, 8.390785e-01, 1.478588e-05, 1e-5, 4.817349e-04)


errs1 <- numeric(length(seeds_all))

for(i in seq_along(seeds_all)) {
  s <- seeds_all[i]
  set.seed(s); W1 <- rnorm(n_time); W2 <- rnorm(n_time)
  I1 <- simulate_sir_sde(n_time,1, pars_opt, initials, W1, W2)[,3]
  errs1[i] <- sum((I1 - observations)^2)  
  
}

e1_min <- min(errs1); e1_max <- max(errs1)

score <- numeric(length(seeds_all))
for(i in seq_along(seeds_all)) {
  
  n1 <- (errs1[i] - e1_min) / (e1_max - e1_min)
  score[i] <- n1
}

best_seeds <- seeds_all[order(score)[1:50] ]

#write.csv(best_seeds, 'best_seeds SEIR.csv', row.names = F)

