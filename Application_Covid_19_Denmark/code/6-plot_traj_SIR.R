setwd('~/Application_Covid_19_Denmark/estimates')

library(cowplot)
library(latex2exp)
###############################################################################
# plot SIR traj with fixed random seeds
# p1: optimal parameter estimates
# p3: alpha = 0.4
# p3: alpha = 0.9
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
initials <- c(s= 1 - i0 - r0,  i = i0)

data1 <- data.frame(df_infected$Date_reported)
colnames(data1) <- c('Time')
data1$observations <- df_infected$approx_infectious_cases/N
data1$Time <- as.Date(data1$Time)


simulate_sir_sde <- function(N, h, par, start, W1, W2) {
  alpha = par[1]; beta = par[2]
  sigma1.sq <- par[3]^2; sigma2.sq <- par[4]^2
  S <- numeric(N); I <- numeric(N)
  S[1] <- start[1]; I[1] <- start[2]
  for (i in 2:N) {
    S[i] <- S[i-1] + h * (-alpha * S[i-1] * I[i-1]) + sqrt(h * sigma1.sq) * W1[i]
    I[i] <- I[i-1] + h * (alpha * S[i-1] * I[i-1] - beta * I[i-1]) + sqrt(h * sigma2.sq) * W2[i]
    if(I[i] <= 0)I[i] = 0
  }
  return(matrix(c(S, I), ncol = 2))
}

########## ode
SIR_fun <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    ds    <- -alpha  * s * i
    di <- alpha  * s * i - beta * i
    
    return(list(c(ds, di)))
  })
}

pars.ode <- c(6.967358e-01, 6.324725e-01, 3.183146e-05)
times <- seq(0, length(observations)-1, by = 0.5)
yini  <- c(s= 1-pars.ode[3]-r0, i = pars.ode[3])

pars1 <- list(alpha = pars.ode[1], beta = pars.ode[2])
out1   <- ode(yini, times, SIR_fun, pars1)
df_ode_sir1 <- data.frame(out1)

df_ode_sir_matched1<- df_ode_sir1 %>% filter(time %in% seq(0, length(observations)-1, by = 1))
df_ode_sir_matched1$time <- data1$Time
df_ode_sir_matched1$legend_label <- "Deterministic SIR model"

#====
pars3 <- list(alpha =0.4, beta = 0.4/(pars.ode[1]/pars.ode[2]))

times <- seq(0, length(observations)-1, by = 0.5)
out3   <- ode(yini, times, SIR_fun, pars3)
df_ode_sir3 <- data.frame(out3)

df_ode_sir_matched3 <- df_ode_sir3 %>% filter(time %in% seq(0, length(observations)-1, by = 1))
df_ode_sir_matched3$time <- data1$Time
df_ode_sir_matched3$legend_label <- "Deterministic SIR model"

#====
pars4 <- list(alpha =0.9, beta = 0.9/(pars.ode[1]/pars.ode[2]))

times <- seq(0, length(observations)-1, by = 0.5)
out4   <- ode(yini, times, SIR_fun, pars4)
df_ode_sir4 <- data.frame(out4)

df_ode_sir_matched4 <- df_ode_sir4 %>% filter(time %in% seq(0, length(observations)-1, by = 1))
df_ode_sir_matched4$time <- data1$Time
df_ode_sir_matched4$legend_label <- "Deterministic SIR model"

#==============================================================================#
best_seeds <- read.csv('best_seeds SIR.csv')[,1]

pars_opt <- c(7.365011e-01, 6.401112e-01, 3.047673e-02, 5.575764e-05)
pars1 <- c(0.4,0.4/(pars_opt[1]/pars_opt[2]), 2.150762e-02, 7.845494e-05)
pars2 <- c(0.9,0.9/(pars_opt[1]/pars_opt[2]), 3.258088e-02, 6.404162e-05)


param_list <- list(
  list(par = pars_opt, label = "Optimal a"),
  list(par = pars1,    label = "a = 0.4"),
  list(par = pars2,    label = "a = 0.9")
)


sim_all <- lapply(param_list, function(pinfo) {
  sim_dfs <- lapply(seq_along(best_seeds), function(i) {
    s <- best_seeds[i]
    set.seed(s); W1 <- rnorm(n_time); W2 <- rnorm(n_time)
    I_sim <- simulate_sir_sde(n_time, 1, pinfo$par, initials, W1, W2)[,2]
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


summary_sim$legend_label <- "Stochastic SIR model - range"
data1$legend_label <- "Observations"
df_ode_sir_matched$legend_label <- "Deterministic SIR model"
summary_sim$legend_label_median <- "Stochastic SIR model - median"

summary_sim1 <- summary_sim %>% filter(label == "Optimal a")
summary_sim3 <- summary_sim %>% filter(label == "a = 0.4")
summary_sim4 <- summary_sim %>% filter(label == "a = 0.9")

p1 <- ggplot() +
  geom_ribbon(data = summary_sim1, aes(x = Time, ymin = min_I , ymax = max_I, fill = legend_label), alpha = 0.3) +
  geom_line(data = summary_sim1, aes(x = Time, y = median_I , color = legend_label_median), size = 1) +
  geom_line(data = data1, aes(x = Time, y = observations, color = legend_label), size = 1) +
  geom_line(data = df_ode_sir_matched1, aes(x = time, y = i, color = legend_label), size = 1) +
  
  scale_color_manual(
    name = "Optimal parameter estimates",
    values = c(
      "Observations" = "black",
      "Deterministic SIR model" = "blue",
      "Stochastic SIR model - median" = "red"
    )
  ) +
  scale_fill_manual(
    name = "",
    values = c("Stochastic SIR model - range" = "darkgray")
  ) +
  labs(
    y = "(Proportion of) Infectious"
  ) +
  scale_x_date(date_labels = "%Y-%m-%d") +
  theme_bw() +
  theme(
    #axis.text.y=element_blank(), axis.title.y=element_blank(),
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
p1 <- ggdraw(p1) +
  draw_label("A", x = 0.025, y = 0.97, fontface = 'bold', size = 30, hjust = 0, vjust = 1)

#************************#
p3 <- ggplot() +
  geom_ribbon(data = summary_sim3, aes(x = Time, ymin = min_I , ymax = max_I, fill = legend_label), alpha = 0.3) +
  geom_line(data = summary_sim3, aes(x = Time, y = median_I , color = legend_label_median), size = 1) +
  geom_line(data = data1, aes(x = Time, y = observations, color = legend_label), size = 1) +
  geom_line(data = df_ode_sir_matched3, aes(x = time, y = i, color = legend_label), size = 1) +
  
  scale_color_manual(
    name = TeX("$\\alpha = 0.4$"),
    values = c(
      "Observations" = "black",
      "Deterministic SIR model" = "blue",
      "Stochastic SIR model - median" = "red"
    )
  ) +
  scale_fill_manual(
    name = "",
    values = c("Stochastic SIR model - range" = "darkgray")
  ) +
  labs(
    y = "(Proportion of) Infectious"
  ) +
  scale_x_date(date_labels = "%Y-%m-%d") +
  theme_bw() +
  theme(
    #axis.text.y=element_blank(), axis.title.y=element_blank(),
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
p3 <-  ggdraw(p3) +
  draw_label("C", x = 0.025, y = 0.97, fontface = 'bold', size = 30, hjust = 0, vjust = 1)


#************************#
p4 <- ggplot() +
  geom_ribbon(data = summary_sim4, aes(x = Time, ymin = min_I , ymax = max_I, fill = legend_label), alpha = 0.3) +
  geom_line(data = summary_sim4, aes(x = Time, y = median_I , color = legend_label_median), size = 1) +
  geom_line(data = data1, aes(x = Time, y = observations, color = legend_label), size = 1) +
  geom_line(data = df_ode_sir_matched4, aes(x = time, y = i, color = legend_label), size = 1) +
  
  scale_color_manual(
    name = TeX("$\\alpha = 0.9$"),
    values = c(
      "Observations" = "black",
      "Deterministic SIR model" = "blue",
      "Stochastic SIR model - median" = "red"
    )
  ) +
  scale_fill_manual(
    name = "",
    values = c("Stochastic SIR model - range" = "darkgray")
  ) +
  labs(
    y = "(Proportion of) Infectious"
  ) +
  scale_x_date(date_labels = "%Y-%m-%d") +
  theme_bw() +
  theme(
    axis.text.y=element_blank(), axis.title.y=element_blank(),
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
p4 <- ggdraw(p4) +
  draw_label("D", x = 0.025, y = 0.97, fontface = 'bold', size = 30, hjust = 0, vjust = 1)



(p1 + p2) / (p3 + p4)

################################################################################
# ============================================================================ #
################################################################################
# The code to generate seeds of the trajectories

seeds_all <- sample.int(1e8, 1e5)

pars_opt <- c(7.365011e-01, 6.401112e-01, 3.047673e-02, 5.575764e-05)
pars1 <- c(0.4,0.4/( pars_opt[1] /  pars_opt[2]), 2.150762e-02, 7.845494e-05)
pars2 <- c(0.9,0.9/( pars_opt[1] /  pars_opt[2]), 3.258088e-02, 6.404162e-05)

# 1) Precompute the distribution of e1,e2,e3 so you know their ranges:
errs1 <- numeric(length(seeds_all))

for(i in seq_along(seeds_all)) {
  s <- seeds_all[i]
  set.seed(s); W1 <- rnorm(n_time); W2 <- rnorm(n_time)
  I1 <- simulate_sir_sde(n_time,1, pars_opt, initials, W1, W2)[,2]

  errs1[i] <- sum((I1 - observations)^2)
}

e1_min <- min(errs1); e1_max <- max(errs1)

score <- numeric(length(seeds_all))
for(i in seq_along(seeds_all)) {
  
  n1 <- (errs1[i] - e1_min) / (e1_max - e1_min)
  score[i] <- n1
}

best_seeds <- seeds_all[order(score)[1:2000] ]

data1 <- data.frame(df_infected$Date_reported)
colnames(data1) <- c('Time')
data1$observations <- df_infected$approx_infectious_cases/N

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
param_list <- list(
  list(par = pars_opt, label = "Optimal a"),
  list(par = pars1,    label = "a = 0.4"),
  list(par = pars2,    label = "a = 0.9")
)


errs2 <- numeric(length(best_seeds))
errs3 <- numeric(length(best_seeds))

for(i in 1:length(best_seeds)) {
  s <- best_seeds[i]
  set.seed(s); W1 <- rnorm(n_time); W2 <- rnorm(n_time)
  I2 <- simulate_sir_sde(n_time,1, pars1,    initials, W1, W2)[,2]
  set.seed(s); W1 <- rnorm(n_time); W2 <- rnorm(n_time)
  I3 <- simulate_sir_sde(n_time,1, pars2,    initials, W1, W2)[,2]
  
  errs2[i] <- sum((I2 - observations)^2)
  errs3[i] <- sum((I3 - observations)^2)
}

e2_min <- min(errs2); e2_max <- max(errs2)
e3_min <- min(errs3); e3_max <- max(errs3)

score <- numeric(length(best_seeds))
for(i in 1:length(best_seeds)) {

  n2 <- (errs2[i] - e2_min) / (e2_max - e2_min)
  n3 <- (errs3[i] - e3_min) / (e3_max - e3_min)
  
  score[i] <- n2 + n3
}


idx <- c(1:50) 

best_seeds <- best_seeds[order(score)[idx]]

#write.csv(best_seeds, 'best_seeds SIR.csv', row.names = F)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
param_list <- list(
  list(par = pars_opt, label = "Optimal a"),
  list(par = pars1,    label = "a = 0.4"),
  list(par = pars2,    label = "a = 0.9")
)

sim_all <- lapply(param_list, function(pinfo) {
  sim_dfs <- lapply(seq_along(best_seeds), function(i) {
    s <- best_seeds[i]
    set.seed(s); W1 <- rnorm(n_time); W2 <- rnorm(n_time)
    I_sim <- simulate_sir_sde(n_time, 1, pinfo$par, initials, W1, W2)[,2]
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


ggplot() +

  geom_ribbon(data = summary_sim,
              aes(x = Time, ymin = min_I, ymax = max_I, fill = "SDE range"),
              alpha = 0.3) +
  geom_line(data = summary_sim,
            aes(x = Time, y = median_I, color = "SDE median"),
            size = 1) +
  geom_line(data = data1,
            aes(x = Time, y = observations, color = "Observations"),
            size = 1) +
  facet_wrap(~ label, ncol = 2, scales = "fixed") +
  scale_color_manual(name = "",
                     values = c("Observations" = "black",
                                "SDE median"   = "red")) +
  scale_fill_manual(name = "",
                    values = c("SDE range" = "gray")) +
  labs(x = "Time", y = "Proportion Infectious") +
  theme_bw() +
  theme(
    strip.text      = element_text(size = 14),
    axis.title      = element_text(size = 14),
    axis.text       = element_text(size = 12),
    legend.position = "bottom",
    legend.text     = element_text(size = 12)
  ) 

