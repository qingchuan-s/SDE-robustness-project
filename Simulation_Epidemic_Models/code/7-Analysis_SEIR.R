setwd('Simulation_Epidemic_Models/estimates')

mat_ode_lse <- read.csv('est_ode_lse_seir.csv')
mat_sde_lse <- read.csv('est_sde_lse_seir.csv')
mat_ode_mle <- read.csv('est_ode_mle_seir.csv')
mat_sde_mle <- read.csv('est_sde_mle_seir.csv')

N_observation = 1000

df_p1 <- data.frame(matrix(c(mat_ode_lse[1:N_observation,1],
                             mat_ode_mle[1:N_observation,1], 
                             mat_sde_lse[1:N_observation,1], 
                             mat_sde_mle[1:N_observation,1]),ncol = 1))
df_p2 <- data.frame(matrix(c(mat_ode_lse[1:N_observation,2],
                             mat_ode_mle[1:N_observation,2], 
                             mat_sde_lse[1:N_observation,2], 
                             mat_sde_mle[1:N_observation,2]),ncol = 1))
df_R <- data.frame(matrix(c(mat_ode_lse[1:N_observation,1]/mat_ode_lse[1:N_observation,2],
                            mat_ode_mle[1:N_observation,1]/mat_ode_mle[1:N_observation,2], 
                            mat_sde_lse[1:N_observation,1]/mat_sde_lse[1:N_observation,2], 
                            mat_sde_mle[1:N_observation,1]/mat_sde_mle[1:N_observation,2]),ncol = 1))

colnames(df_p1) <-  colnames(df_p2) <-  colnames(df_R) <-  c("Estimates") 
df_p1$`Data-generating model` <- df_p2$`Data-generating model` <- df_R$`Data-generating model` <- 
  c(rep('deterministic SEIR model', 2*N_observation),
    rep('stochastic SEIR model', 2*N_observation))
df_p1$`Fitted model` <- df_p2$`Fitted model` <-  df_R$`Fitted model` <-
  c(rep('deterministic SIR model', N_observation),
    rep('stochastic SIR model', N_observation),
    rep('deterministic SIR model', N_observation),
    rep('stochastic SIR model', N_observation))

#==============================================================================#
# Parameter with original data set
df_p1.ode <- df_p1[(df_p1$`Data-generating model` == 'deterministic SEIR model'),] 
df_p2.ode <- df_p2[(df_p2$`Data-generating model` == 'deterministic SEIR model'),]
df_R.ode <- df_R[(df_R$`Data-generating model` == 'deterministic SEIR model'),] 


df_p1.sde <- df_p1[(df_p1$`Data-generating model` == 'stochastic SEIR model'),]
df_p2.sde <- df_p2[(df_p2$`Data-generating model` == 'stochastic SEIR model'),]
df_R.sde <- df_R[(df_R$`Data-generating model` == 'stochastic SEIR model'),]

p1 <- ggplot(data = df_p1.ode)+
  geom_density(aes(x= Estimates, color = `Fitted model`),fill='grey',alpha = 0.1,adjust = 3, linewidth =0.8)+
  geom_vline(xintercept = 0.5/1.3, linetype = 2, linewidth = 0.5)+
  scale_color_manual("Fitted model",labels = c('deterministic SIR model', 'stochastic SIR model'),
                     values = c(4,7)) +
  facet_grid(`Data-generating model`~ .)+
  xlim(0.28,0.58)+ylim(0,130)+
  labs(x = NULL, y=NULL) + 
  theme(legend.position = "none",
        strip.text.y = element_blank(),
        axis.text.x=element_blank())+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )

p2 <- ggplot(data = df_p2.ode)+
  geom_density(aes(x= Estimates, color = `Fitted model`),fill='grey',alpha = 0.1,adjust = 3, linewidth =0.8)+
  geom_vline(xintercept = 0.3/1.3, linetype = 2, linewidth = 0.5)+
  scale_color_manual("Fitted model",labels = c('deterministic SIR model', 'stochastic SIR model'),
                     values = c(4,7)) +
  facet_grid(`Data-generating model`~ .)+
  xlim(0.15,0.36)+ylim(0,130)+
  labs(x = NULL, y=NULL) + theme(legend.position = "none",strip.text.y = element_blank(),
                                 axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )

p3 <- ggplot(data = df_p1.sde)+
  geom_density(aes(x= Estimates, color = `Fitted model`),fill='grey',alpha = 0.1,adjust = 3, linewidth =0.8)+
  geom_vline(xintercept = 0.5/1.3, linetype = 2, linewidth = 0.5)+
  scale_color_manual("Fitted model",labels = c('deterministic SIR model', 'stochastic SIR model'),
                     values = c(4,7)) +
  facet_grid(`Data-generating model`~ .)+
  xlim(0.28,0.58)+ylim(0,40)+
  labs(x = NULL, y=NULL) + 
  theme(legend.position = "none",
        strip.text.y = element_blank())+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )
p4 <- ggplot(data = df_p2.sde)+
  geom_density(aes(x= Estimates, color = `Fitted model`),fill='grey',alpha = 0.1,adjust = 3, linewidth =0.8)+
  geom_vline(xintercept = 0.3/1.3, linetype = 2, linewidth = 0.5)+
  scale_color_manual("Fitted model",labels = c('deterministic SIR model', 'stochastic SIR model'),
                     values = c(4,7)) +
  facet_grid(`Data-generating model`~ .)+
  xlim(0.15,0.36)+ylim(0,40)+
  labs(x = NULL, y=NULL) + 
  theme(legend.position = "none",
        strip.text.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )

p.R1 <- ggplot(data = df_R.ode )+
  geom_density(aes(x= Estimates, color = `Fitted model`), fill= 'grey',alpha = 0.1,adjust = 3, linewidth =0.8)+
  geom_vline(xintercept = 5/3, linetype = 2, linewidth = 0.5)+
  facet_grid(`Data-generating model`~ .)+
  scale_color_manual("Fitted model",labels = c('deterministic SIR model', 'stochastic SIR model'),
                     values = c(4,7)) +
  xlim(1.42,1.83)+ylim(0,130)+
  labs(x = NULL, y=NULL) + 
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.text.y=element_blank())+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )

p.R2 <-ggplot(data = df_R.sde )+
  geom_density(aes(x= Estimates, color = `Fitted model`), fill= 'grey',alpha = 0.1,adjust = 3, linewidth =0.8)+
  geom_vline(xintercept = 5/3, linetype = 2, linewidth = 0.5)+
  facet_grid( `Data-generating model`~ .)+
  scale_color_manual("Fitted model",labels = c('deterministic SIR model', 'stochastic SIR model'),
                     values = c(4,7)) +
  xlim(1.42,1.83)+ylim(0,40)+
  labs(x = NULL, y=NULL) + theme(axis.text.y=element_blank())+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )

p_lab1 <- 
  ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = TeX("$\\hat{\\alpha}$ "), angle = 0) +
  coord_cartesian(clip = "off")+
  theme_void()
p_lab2 <- 
  ggplot() + 
  annotate(geom = "text", x = 2, y = 1, label = TeX("$\\hat{\\beta}$"), angle = 0) +
  coord_cartesian(clip = "off")+
  theme_void()

p_lab3 <-   ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = TeX("$\\hat{R}_0$ "), angle = 0) +
  coord_cartesian(clip = "off")+
  theme_void()

p_lab4 <- 
  ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = "Data-generating Model", angle = -90) +
  coord_cartesian(clip = "off")+
  theme_void()


design <- "
BBBBBCCCCCFFFFFA
DDDDDEEEEEGGGGGA
HHHHHIIIIIJJJJJ#
"
p_lab4 + p1+p2+ p3+p4+p.R1+p.R2+p_lab1 + p_lab2+ p_lab3+
  plot_layout(heights = c( 1,1,0.06),design = design,guides = "collect")& theme(legend.position = "top")


