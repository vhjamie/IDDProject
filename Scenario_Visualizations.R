library(ggplot2)
library(tidyverse)
library(stringr)
library(cowplot)
# Scenario 1: Status quo (2021 parameters)
scenario1 <- run_simp_coupled_HIV_function(beta = 0.0008,
                                           c_GS = 950,
                                           c_GG = 54,
                                           c_SG = 67.6,
                                           sigma_GS = 0.36,
                                           sigma_GG = 0.27,
                                           sigma_SG = 0.18,
                                           alpha = 0.009532, 
                                           theta = 0.00001, 
                                           delta = 0.007, 
                                           mu1 = 0.2149,
                                           mu2 = 0.04228,
                                           tau_prep_g = 0.0000214,
                                           tau_prep_s = 0.00583,
                                           phi_prep_g = 0.015, # Need to test for sensitivity of a range of retention rates
                                           phi_prep_s = 0.015, # 
                                           tau_art_g = 0.3605,
                                           phi_art_g = 0.015,
                                           tau_art_s = 0.3605,
                                           phi_art_s = 0.015,
                                           initial.state = initial.state,
                                           max.time = 5,
                                           freq.dependent = TRUE)

# Incidence data for scenario 1
scenario1 <- add_column(scenario1, "incidence_s" = c(rep.int(0, nrow(scenario1))))
for (i in (2:6)) {
  beta = 0.0008
  c_GS = 950
  sigma_GS = 0.36
  scenario1$incidence_s[i] = (0.5*beta*c_GS*sigma_GS*scenario1$I_g[i-1]/(scenario1$I_s[i-1]+scenario1$A_s[i-1]+scenario1$S_s[i-1]+scenario1$PrEP_s[i-1]))*scenario1$S_s[i-1]
}

scenario1 <- add_column(scenario1, "incidence_g" = c(rep.int(0, nrow(scenario1))))
for (i in (2:6)) {
  beta = 0.0008
  c_GG = 54
  c_SG = 67.6
  sigma_GG = 0.27
  sigma_SG = 0.18
  scenario1$incidence_g[i] = (beta*c_GG*sigma_GG*scenario1$I_g[i-1] + 0.5*beta*c_SG*sigma_SG*scenario1$I_s[i-1])/(scenario1$I_g[i-1]+scenario1$A_g[i-1]+scenario1$S_g[i-1]+scenario1$PrEP_g[i-1]) * scenario1$S_g[i-1]
}
scenario1$D <- scenario1$D_g + scenario1$D_s

# Scenario 2: ART Scale up 
## 2a: Targeted intervention (30% interpolation)
scenario2a <- run_simp_coupled_HIV_function(beta = 0.0008,
                                           c_GS = 950,
                                           c_GG = 54,
                                           c_SG = 67.6,
                                           sigma_GS = 0.36,
                                           sigma_GG = 0.27,
                                           sigma_SG = 0.18,
                                           alpha = 0.009532, 
                                           theta = 0.00001, 
                                           delta = 0.007, 
                                           mu1 = 0.2149,
                                           mu2 = 0.04228,
                                           tau_prep_g = 0.0000214,
                                           tau_prep_s = 0.00583,
                                           phi_prep_g = 0.015, 
                                           phi_prep_s = 0.015, 
                                           tau_art_g = 0.3605,
                                           phi_art_g = 0.015,
                                           tau_art_s = 0.46865,
                                           phi_art_s = 0.015,
                                           initial.state = initial.state,
                                           max.time = 5,
                                           freq.dependent = TRUE)

scenario2a <- add_column(scenario2a, "incidence_s" = c(rep.int(0, nrow(scenario2a))))
for (i in (2:6)) {
  beta = 0.0008
  c_GS = 950
  sigma_GS = 0.36
  scenario2a$incidence_s[i] = (0.5*beta*c_GS*sigma_GS*scenario2a$I_g[i-1]/(scenario2a$I_s[i-1]+scenario2a$A_s[i-1]+scenario2a$S_s[i-1]+scenario2a$PrEP_s[i-1]))*scenario2a$S_s[i-1]
}

scenario2a <- add_column(scenario2a, "incidence_g" = c(rep.int(0, nrow(scenario2a))))
for (i in (2:6)) {
  beta = 0.0008
  c_GG = 54
  c_SG = 67.6
  sigma_GG = 0.27
  sigma_SG = 0.18
  scenario2a$incidence_g[i] = (beta*c_GG*sigma_GG*scenario2a$I_g[i-1] + 0.5*beta*c_SG*sigma_SG*scenario2a$I_s[i-1])/(scenario2a$I_g[i-1]+scenario2a$A_g[i-1]+scenario2a$S_g[i-1]+scenario2a$PrEP_g[i-1]) * scenario2a$S_g[i-1]
}

## 2b: Targeted intervention (50% scale up)
scenario2b <- run_simp_coupled_HIV_function(beta = 0.0008,
                                            c_GS = 950,
                                            c_GG = 54,
                                            c_SG = 67.6,
                                            sigma_GS = 0.36,
                                            sigma_GG = 0.27,
                                            sigma_SG = 0.18,
                                            alpha = 0.009532, 
                                            theta = 0.00001, 
                                            delta = 0.007, 
                                            mu1 = 0.2149,
                                            mu2 = 0.04228,
                                            tau_prep_g = 0.0000214,
                                            tau_prep_s = 0.00583,
                                            phi_prep_g = 0.015, 
                                            phi_prep_s = 0.015, 
                                            tau_art_g = 0.3605,
                                            phi_art_g = 0.015,
                                            tau_art_s = 0.54075,
                                            phi_art_s = 0.015,
                                            initial.state = initial.state,
                                            max.time = 5,
                                            freq.dependent = TRUE)

scenario2b <- add_column(scenario2b, "incidence_s" = c(rep.int(0, nrow(scenario2b))))
for (i in (2:6)) {
  beta = 0.0008
  c_GS = 950
  sigma_GS = 0.36
  scenario2b$incidence_s[i] = (0.5*beta*c_GS*sigma_GS*scenario2b$I_g[i-1]/(scenario2b$I_s[i-1]+scenario2b$A_s[i-1]+scenario2b$S_s[i-1]+scenario2b$PrEP_s[i-1]))*scenario2b$S_s[i-1]
}

scenario2b <- add_column(scenario2b, "incidence_g" = c(rep.int(0, nrow(scenario2b))))
for (i in (2:6)) {
  beta = 0.0008
  c_GG = 54
  c_SG = 67.6
  sigma_GG = 0.27
  sigma_SG = 0.18
  scenario2b$incidence_g[i] = (beta*c_GG*sigma_GG*scenario2b$I_g[i-1] + 0.5*beta*c_SG*sigma_SG*scenario2b$I_s[i-1])/(scenario2b$I_g[i-1]+scenario2b$A_g[i-1]+scenario2b$S_g[i-1]+scenario2b$PrEP_g[i-1]) * scenario2b$S_g[i-1]
}

##2c: Targeted intervention (70% scale up)
scenario2c <- run_simp_coupled_HIV_function(beta = 0.0008,
                                            c_GS = 950,
                                            c_GG = 54,
                                            c_SG = 67.6,
                                            sigma_GS = 0.36,
                                            sigma_GG = 0.27,
                                            sigma_SG = 0.18,
                                            alpha = 0.009532, 
                                            theta = 0.00001, 
                                            delta = 0.007, 
                                            mu1 = 0.2149,
                                            mu2 = 0.04228,
                                            tau_prep_g = 0.0000214,
                                            tau_prep_s = 0.00583,
                                            phi_prep_g = 0.015, 
                                            phi_prep_s = 0.015, 
                                            tau_art_g = 0.3605,
                                            phi_art_g = 0.015,
                                            tau_art_s = 0.61285,
                                            phi_art_s = 0.015,
                                            initial.state = initial.state,
                                            max.time = 5,
                                            freq.dependent = TRUE)

scenario2c <- add_column(scenario2c, "incidence_s" = c(rep.int(0, nrow(scenario2c))))
for (i in (2:6)) {
  beta = 0.0008
  c_GS = 950
  sigma_GS = 0.36
  scenario2c$incidence_s[i] = (0.5*beta*c_GS*sigma_GS*scenario2c$I_g[i-1]/(scenario2c$I_s[i-1]+scenario2c$A_s[i-1]+scenario2c$S_s[i-1]+scenario2c$PrEP_s[i-1]))*scenario2c$S_s[i-1]
}

scenario2c <- add_column(scenario2c, "incidence_g" = c(rep.int(0, nrow(scenario2c))))
for (i in (2:6)) {
  beta = 0.0008
  c_GG = 54
  c_SG = 67.6
  sigma_GG = 0.27
  sigma_SG = 0.18
  scenario2c$incidence_g[i] = (beta*c_GG*sigma_GG*scenario2c$I_g[i-1] + 0.5*beta*c_SG*sigma_SG*scenario2c$I_s[i-1])/(scenario2c$I_g[i-1]+scenario2c$A_g[i-1]+scenario2c$S_g[i-1]+scenario2c$PrEP_g[i-1]) * scenario2c$S_g[i-1]
}

##2d: General intervention (10% scale up for the entire population)
scenario2d <- run_simp_coupled_HIV_function(beta = 0.0008,
                                            c_GS = 950,
                                            c_GG = 54,
                                            c_SG = 67.6,
                                            sigma_GS = 0.36,
                                            sigma_GG = 0.27,
                                            sigma_SG = 0.18,
                                            alpha = 0.009532, 
                                            theta = 0.00001, 
                                            delta = 0.007, 
                                            mu1 = 0.2149,
                                            mu2 = 0.04228,
                                            tau_prep_g = 0.0000214,
                                            tau_prep_s = 0.00583,
                                            phi_prep_g = 0.015, 
                                            phi_prep_s = 0.015, 
                                            tau_art_g = 0.39655,
                                            phi_art_g = 0.015,
                                            tau_art_s = 0.39655,
                                            phi_art_s = 0.015,
                                            initial.state = initial.state,
                                            max.time = 5,
                                            freq.dependent = TRUE)

scenario2d <- add_column(scenario2d, "incidence_s" = c(rep.int(0, nrow(scenario2d))))
for (i in (2:6)) {
  beta = 0.0008
  c_GS = 950
  sigma_GS = 0.36
  scenario2d$incidence_s[i] = (0.5*beta*c_GS*sigma_GS*scenario2d$I_g[i-1]/(scenario2d$I_s[i-1]+scenario2d$A_s[i-1]+scenario2d$S_s[i-1]+scenario2d$PrEP_s[i-1]))*scenario2d$S_s[i-1]
}

scenario2d <- add_column(scenario2d, "incidence_g" = c(rep.int(0, nrow(scenario2d))))
for (i in (2:6)) {
  beta = 0.0008
  c_GG = 54
  c_SG = 67.6
  sigma_GG = 0.27
  sigma_SG = 0.18
  scenario2d$incidence_g[i] = (beta*c_GG*sigma_GG*scenario2d$I_g[i-1] + 0.5*beta*c_SG*sigma_SG*scenario2d$I_s[i-1])/(scenario2d$I_g[i-1]+scenario2d$A_g[i-1]+scenario2d$S_g[i-1]+scenario2d$PrEP_g[i-1]) * scenario2d$S_g[i-1]
}

# Scenario 3: Prediction for policy change reducing KP-led PrEP
scenario3 <- run_simp_coupled_HIV_function(beta = 0.0008,
                                            c_GS = 950,
                                            c_GG = 54,
                                            c_SG = 67.6,
                                            sigma_GS = 0.36,
                                            sigma_GG = 0.27,
                                            sigma_SG = 0.18,
                                            alpha = 0.009532, 
                                            theta = 0.00001, 
                                            delta = 0.007, 
                                            mu1 = 0.2149,
                                            mu2 = 0.04228,
                                            tau_prep_g = 0.0000214,
                                            tau_prep_s = 0.0000214,
                                            phi_prep_g = 0.015, 
                                            phi_prep_s = 0.025, 
                                            tau_art_g = 0.39655,
                                            phi_art_g = 0.015,
                                            tau_art_s = 0.39655,
                                            phi_art_s = 0.015,
                                            initial.state = initial.state,
                                            max.time = 5,
                                            freq.dependent = TRUE)

scenario3 <- add_column(scenario3, "incidence_s" = c(rep.int(0, nrow(scenario3))))
for (i in (2:6)) {
  beta = 0.0008
  c_GS = 950
  sigma_GS = 0.36
  scenario3$incidence_s[i] = (0.5*beta*c_GS*sigma_GS*scenario3$I_g[i-1]/(scenario3$I_s[i-1]+scenario3$A_s[i-1]+scenario3$S_s[i-1]+scenario3$PrEP_s[i-1]))*scenario3$S_s[i-1]
}

scenario3 <- add_column(scenario3, "incidence_g" = c(rep.int(0, nrow(scenario3))))
for (i in (2:6)) {
  beta = 0.0008
  c_GG = 54
  c_SG = 67.6
  sigma_GG = 0.27
  sigma_SG = 0.18
  scenario3$incidence_g[i] = (beta*c_GG*sigma_GG*scenario3$I_g[i-1] + 0.5*beta*c_SG*sigma_SG*scenario3$I_s[i-1])/(scenario3$I_g[i-1]+scenario3$A_g[i-1]+scenario3$S_g[i-1]+scenario3$PrEP_g[i-1]) * scenario3$S_g[i-1]
}
scenario3$D <- scenario3$D_g + scenario3$D_s
## Visualization
## Deaths averted
deathsAverted_df_g <- data_frame("scenario1" = scenario1$D_g, 
                                   "scenario2a" = scenario2a$D_g,
                                   "scenario2b" = scenario2b$D_g,
                                   "scenario2c" = scenario2c$D_g,
                                   "scenario2d" = scenario2d$D_g)
deathsAverted_df_g$scenario2a_averted = (deathsAverted_df_g$scenario1 - deathsAverted_df_g$scenario2a)
deathsAverted_df_g$scenario2b_averted = (deathsAverted_df_g$scenario1 - deathsAverted_df_g$scenario2b)
deathsAverted_df_g$scenario2c_averted = (deathsAverted_df_g$scenario1 - deathsAverted_df_g$scenario2c)
deathsAverted_df_g$scenario2d_averted = (deathsAverted_df_g$scenario1 - deathsAverted_df_g$scenario2d)

deathsAverted_df_s <- data_frame("scenario1" = scenario1$D_s, 
                                 "scenario2a" = scenario2a$D_s,
                                 "scenario2b" = scenario2b$D_s,
                                 "scenario2c" = scenario2c$D_s,
                                 "scenario2d" = scenario2d$D_s)
deathsAverted_df_s$scenario2a_averted = (deathsAverted_df_s$scenario1 - deathsAverted_df_s$scenario2a)
deathsAverted_df_s$scenario2b_averted = (deathsAverted_df_s$scenario1 - deathsAverted_df_s$scenario2b)
deathsAverted_df_s$scenario2c_averted = (deathsAverted_df_s$scenario1 - deathsAverted_df_s$scenario2c)
deathsAverted_df_s$scenario2d_averted = (deathsAverted_df_s$scenario1 - deathsAverted_df_s$scenario2d)

deathsAverted_df_all <- data_frame("time" = c(2021:2026),
                                   "scenario2a" = deathsAverted_df_g$scenario2a_averted + deathsAverted_df_s$scenario2a_averted,
                                   "scenario2b" = deathsAverted_df_g$scenario2b_averted + deathsAverted_df_s$scenario2b_averted,
                                   "scenario2c" = deathsAverted_df_g$scenario2c_averted + deathsAverted_df_s$scenario2c_averted,
                                   "scenario2d" = deathsAverted_df_g$scenario2d_averted + deathsAverted_df_s$scenario2d_averted)

# Plotting deaths Averted Scenario 2a-2d
ggplot(deathsAverted_df_all, aes(x= time)) +
  geom_line(aes(y = scenario2a, colour = "Targeted Intervention - 30%")) +
  geom_line(aes(y = scenario2b, colour = "Targeted Intervention - 50%")) +
  geom_line(aes(y = scenario2c, colour = "Targeted Intervention - 70%")) +
  geom_line(aes(y = scenario2d, colour = "Generalized intervention")) +
  labs(title = "Projected deaths averted by different ART intervention strategies in Thailand",
       x = "Time (years)",
       y = "Deaths averted") + 
  scale_color_manual(name = "Legend",
                     values = c("Targeted Intervention - 30%"="blue",
                                   "Targeted Intervention - 50%"="red",
                                   "Targeted Intervention - 70%"="green",
                                   "Generalized intervention"="purple")) 


# Plotting comparison _ Status quo our model vs. SPECTRUM AIM Model
incidence_df <- data_frame("time" = c(2022:2026),
                           "Our_model" = scenario1$incidence_g[2:6] + scenario1$incidence_s[2:6],
                           "AIMmodel" = c(7208.816347785235,	7150.787396832486,	7095.767046566996,	7046.995582902279, 7000.304214385396))

ggplot(incidence_df, aes(x=time)) +
  geom_line(aes(y = Our_model, colour = "Our model projection")) +
  geom_line(aes(y = AIMmodel, colour = "SPECTRUM AIM Model")) +
  labs(title = str_wrap("Comparison of projected HIV incidence by our model and the SPECTRUM AIM Model for Thailand"),
       x = "Time (years)",
       y = "Incidence (cases)") + 
  scale_color_manual(name = "Legend",
                     values = c("Our model projection"="blue",
                                "SPECTRUM AIM Model"="red"))

# Visualize effect of the policy
policyChange_df <- data.frame("time" = c(2022:2026),
                              "scenario1_D" = scenario1$D[2:6],
                              "scenario3_D" = scenario3$D[2:6],
                              #"scenario1_PrEP" = scenario1$PrEP_g[2:6] + scenario1$PrEP_s[2:6],
                              "scenario1_PrEP" =  scenario1$PrEP_s[2:6],
                              "scenario3_PrEP" = scenario3$PrEP_s[2:6] 
                              #"scenario3_PrEP" = scenario3$PrEP_g[2:6] + scenario3$PrEP_s[2:6],
                              "scenario1_S_s" = scenario1$S_s[2:6],
                              "scenario3_S_s" = scenario3$S_s[2:6],
                              "scenario1_inc" = scenario1$incidence_g[2:6] + scenario1$incidence_s[2:6],
                              "scenario3_inc" = scenario3$incidence_g[2:6] + scenario3$incidence_s[2:6])
# PrEP users plot
prep_users_plot <- ggplot(policyChange_df, aes(x=time)) +
  geom_line(aes(y = scenario1_PrEP, colour = "Status quo")) +
  geom_line(aes(y = scenario3_PrEP, colour = "New policies")) +
  labs(title = str_wrap("Comparison of projected PrEP users"),
       x = "Time (years)",
       y = "PrEP users") + 
  scale_color_manual(name = "Legend",
                     values = c("Status quo"="blue",
                                "New policies"="red"))
susceptibility_plot <- ggplot(policyChange_df, aes(x=time)) +
  geom_line(aes(y = scenario1_S_s, colour = "Status quo")) +
  geom_line(aes(y = scenario3_S_s, colour = "New policies")) +
  labs(title = str_wrap("Comparison of the projected susceptibles"),
       x = "Time (years)",
       y = "Susceptibles") + 
  scale_color_manual(name = "Legend",
                     values = c("Status quo"="blue",
                                "New policies"="red"))
incidence_policy_plot <- ggplot(policyChange_df, aes(x=time)) +
  geom_line(aes(y = scenario1_inc, colour = "Status quo")) +
  geom_line(aes(y = scenario3_inc, colour = "New policies")) +
  labs(title = str_wrap("Comparison projected incidence of HIV"),
       x = "Time (years)",
       y = "Incidence (cases)") + 
  scale_color_manual(name = "Legend",
                     values = c("Status quo"="blue",
                                "New policies"="red"))
incidence_policy_plot
legend <- get_legend(
  # create some space to the left of the legend
  incidence_policy_plot + theme(legend.box.margin = margin(0, 0, 0, 12))
)
plot_grid(prep_users_plot + theme(legend.position="none"),
          susceptibility_plot + theme(legend.position="none"),
          incidence_policy_plot + theme(legend.position="none"),
          legend,
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
