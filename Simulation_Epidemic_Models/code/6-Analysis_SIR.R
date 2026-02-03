setwd('Simulation_Epidemic_Models/estimates')

mat_ode_lse <- read.csv('est_ode_lse_xptb.csv')
mat_sde_lse <- read.csv('est_sde_lse_xptb.csv')
mat_ode_mle <- read.csv('est_ode_mle_xptb.csv')
mat_sde_mle <- read.csv('est_sde_mle_xptb.csv')

N_observation  = 1000

#### Full data set ####
df_p1.full <- data.frame(Estimates = mat_ode_lse[1:N_observation,1])
df_p1.full[(N_observation+1):(4*N_observation),] <- c(mat_ode_mle[1:N_observation,1], 
                                                      mat_sde_lse[1:N_observation,1], 
                                                      mat_sde_mle[1:N_observation,1])
df_p1.full$Data_model <- c(rep("Deterministic SIR model", (2*N_observation)),rep("Stochastic SIR model", (2*N_observation)))
df_p1.full$Fit_model <- c(rep("Deterministic SIR model", N_observation),rep("Stochastic SIR model", N_observation),
                          rep("Deterministic SIR model", N_observation),rep("Stochastic SIR model", N_observation))
df_p1.full$par <- "alpha"

df_p2.full <- data.frame(Estimates =mat_ode_lse[1:N_observation,2])
df_p2.full[(N_observation+1):(4*N_observation),] <- c(mat_ode_mle[1:N_observation,2], 
                                                      mat_sde_lse[1:N_observation,2], 
                                                      mat_sde_mle[1:N_observation,2])
df_p2.full$Data_model <- c(rep("Deterministic SIR model", (2*N_observation)),rep("Stochastic SIR model", (2*N_observation)))
df_p2.full$Fit_model <- c(rep("Deterministic SIR model", N_observation),rep("Stochastic SIR model", N_observation),
                          rep("Deterministic SIR model", N_observation),rep("Stochastic SIR model", N_observation))
df_p2.full$par <- "beta"

df_p1.full1 <- df_p1.full%>%filter(Data_model=="Deterministic SIR model")
df_p1.full2 <- df_p1.full%>%filter(Data_model=="Stochastic SIR model")
df_p2.full1 <- df_p2.full%>%filter(Data_model=="Deterministic SIR model")
df_p2.full2 <- df_p2.full%>%filter(Data_model=="Stochastic SIR model")

df_R.full1 <- df_p1.full1 %>%
  mutate(Estimates = Estimates / df_p2.full1$Estimates,
         par = "R")

df_R.full2 <- df_p1.full2 %>%
  mutate(Estimates = Estimates / df_p2.full2$Estimates,
         par = "R")

labels_tex <- c(TeX("Deterministic SIR model"), TeX("Stochastic SIR model"))

p1 <- ggplot(data = df_p1.full1 ,aes(x = Estimates, color=Fit_model)) +
  geom_density(fill="grey",alpha=0.2,linewidth=0.8)+
  geom_vline(aes(xintercept= 0.5),color="black", linetype="dashed", linewidth=0.5)+
  scale_color_manual(name = "Fitted model", labels = labels_tex, values = c(4,7))+
  labs(title="",x= "", y = "")+  theme(legend.position = "none")+
  facet_grid( Data_model~. ) + xlim(c(0.4,0.75)) + ylim(c(0,100)) +
  theme(
    axis.text.x = element_blank(), 
    strip.text.y = element_blank(),
    strip.background = element_blank() 
  ) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )

p2 <- ggplot(data = df_p1.full2 ,aes(x = Estimates, color=Fit_model))+
  geom_density(fill="grey",alpha=0.2,linewidth=0.8)+
  geom_vline(aes(xintercept= 0.5),color="black", linetype="dashed", linewidth=0.5)+
  scale_color_manual(name = "Fitted model", labels = labels_tex, values = c(4,7))+
  labs(title=" ",x= TeX("$\\hat{\\alpha}$"), y = "")+ 
  theme(legend.position = "none", strip.text.y = element_blank(),
        strip.background = element_blank() )+
  facet_grid( Data_model~. ) + xlim(c(0.4,0.75)) + ylim(c(0,45)) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )


p3 <- ggplot(data = df_p2.full1 ,aes(x = Estimates, color=Fit_model))+
  geom_density(fill="grey",alpha=0.2,linewidth=0.8)+
  geom_vline(aes(xintercept= 0.3),color="black", linetype="dashed", linewidth=0.5)+
  scale_color_manual(name = "Fitted model", labels = labels_tex, values = c(4,7))+
  labs(title=" ",x= "", y = "")+  theme(legend.position = "none")+
  facet_grid( Data_model ~ .) + xlim(c(0.22,0.48))+ ylim(c(0,100))+
  theme(
    axis.text.x = element_blank(),  
    axis.text.y = element_blank(),
    strip.text.y = element_blank(),
    strip.background = element_blank() 
  ) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )
p4 <- ggplot(data = df_p2.full2 ,aes(x = Estimates, color=Fit_model))+
  geom_density(fill="grey",alpha=0.2,linewidth=0.8)+
  geom_vline(aes(xintercept= 0.3),color="black", linetype="dashed", linewidth=0.5)+
  scale_color_manual(name = "Fitted model", labels = labels_tex, values = c(4,7))+
  labs(title=" ",x= TeX("$\\hat{\\beta}$"), y = " ")+  theme(legend.position = "none")+
  facet_grid( Data_model ~ .) + xlim(c(0.22,0.48)) + ylim(c(0,45))+
  theme(
    #axis.text.x = element_blank(),  # remove x-axis numbers
    axis.text.y = element_blank(), strip.text.y = element_blank(),
    strip.background = element_blank() 
  ) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )

p5 <- ggplot(data = df_R.full1 ,aes(x = Estimates, color=Fit_model))+
  geom_density(fill="grey",alpha=0.2,linewidth=0.8)+
  geom_vline(aes(xintercept= 5/3),color="black", linetype="dashed", linewidth=0.5)+
  scale_color_manual(name = "Fitted model", labels = labels_tex, values = c(4,7))+
  labs(title=" ",x= "", y = "")+  #theme(legend.position = "none")+
  facet_grid( Data_model ~ .) + xlim(c(1.54,1.75)) + ylim(c(0,100))+
  theme(
    axis.text.x = element_blank(),  # remove x-axis numbers
    axis.text.y = element_blank()
  ) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )
p6 <- ggplot(data = df_R.full2 ,aes(x = Estimates, color=Fit_model))+
  geom_density(fill="grey",alpha=0.2,linewidth=0.8)+
  geom_vline(aes(xintercept= 5/3),color="black", linetype="dashed", linewidth=0.5)+
  scale_color_manual(name = "Fitted model", labels = labels_tex, values = c(4,7))+
  labs(title=" ",x= TeX("$\\hat{R}_0$"), y = " ")+  theme(legend.position = "none")+
  facet_grid( Data_model ~ .) + xlim(c(1.54,1.75)) + ylim(c(0,45))+
  theme(
    #axis.text.x = element_blank(),  # remove x-axis numbers
    axis.text.y = element_blank()
  ) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )


p_lab <- 
  ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = "Data-generating model", angle = -90) +
  coord_cartesian(clip = "off")+
  theme_void()

p <- (p1 + p3 +p5 + p2 + p4 +  p6) + plot_layout(guides = "collect")& theme(legend.position = "top") 

(p|p_lab ) +
  plot_layout(widths = c(1,0.05))

