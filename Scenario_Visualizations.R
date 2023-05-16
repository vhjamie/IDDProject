library(ggplot2)
library(tidyverse)
library(stringr)
library(cowplot)
library(lsoda)

simplified_coupled_HIV_func <- function(t, y, param){
  ## parameters
  # population
  alpha = param['alpha'] # birth rate - more specifically rate of people being introduced to the susceptible population once they turn 18
  theta = param['theta']# rate of starting sex work
  delta = param['delta'] # rate of non-AIDS-related mortality
  # transmission
  beta = param['beta'] 
  c_GS = param['c_GS'] # contact rate G>S
  c_GG = param['c_GG'] # contact rate within general population
  c_SG = param['c_SG'] # contact rate S>G
  sigma_GS = param['sigma_GS'] # parameter for type of contact risk %>% taking into account interventions like condom use etc.
  sigma_GG = param['sigma_GG'] # parameter for type of contact risk %>% taking into account interventions like condom use etc.
  sigma_SG = param['sigma_SG'] # parameter for type of contact risk %>% taking into account interventions like condom use etc.
  betaGS = param['beta'] * param['c_GS'] * param['sigma_GS']
  betaGG = param['beta'] * param['c_GG'] * param['sigma_GG']
  betaSG = param['beta'] * param['c_SG'] * param['sigma_SG']
  # death rates
  mu1 = param['mu1'] # AIDS-related mortality rate _ without ART
  mu2 = param['mu2'] # AIDS-related mortality rate _ with ART
  # prep
  tau_prep_g = param['tau_prep_g']# rate of initiating prep in the general population
  tau_prep_s = param['tau_prep_s']# rate of initiating prep in the SWs population
  phi_prep_g = param['phi_prep_g']# rate of prep discontinuation - general population
  phi_prep_s = param['phi_prep_s']# rate of prep discontinuation - SWs population
  # art
  tau_art_g = param['tau_art_g']# rate of initiating ART
  phi_art_g = param['phi_art_g']# rate of ART discontinuation
  tau_art_s = param['tau_art_s']# rate of initiating ART
  phi_art_s = param['phi_art_s']# rate of ART discontinuation
  ## defining states
  S_g = y['S_g']; I_g = y['I_g']; A_g = y['A_g'] 
  S_s = y['S_s']; I_s = y['I_s']; A_s = y['A_s']
  PrEP_g = y['PrEP_g']; PrEP_s = y['PrEP_s']; D = y['D']
  N = y['S_g'] + y['I_g'] + y['A_g'] + y['S_s'] + y['I_s'] + 
    y['A_s'] + y['PrEP_g'] + y['PrEP_s']
  
  ## Equations for the general population
  dS_g = -(beta*c_GG*sigma_GG*I_g + 0.5*beta*c_SG*sigma_SG*I_s)*S_g - tau_prep_g*S_g - delta*S_g - theta*S_g + alpha*N + phi_prep_g*PrEP_g
  dS_s = -(0.5*betaGS*I_g)*S_s - tau_prep_s*S_s - delta*S_s + theta*S_g + phi_prep_s*PrEP_s
  dPrEP_g = -(phi_prep_g)*PrEP_g + tau_prep_g*S_g - delta*PrEP_g
  dPrEP_s = -(phi_prep_s)*PrEP_s + tau_prep_s*S_s - delta*PrEP_s
  dI_g = -delta*I_g - mu1*I_g - tau_art_g*I_g + (beta*c_GG*sigma_GG*I_g + 0.5*beta*c_SG*sigma_SG*I_s)*S_g + phi_art_g*A_g - theta*I_g
  dI_s = -delta*I_s - mu1*I_s - tau_art_s*I_s + (0.5*beta*c_GS*sigma_GS*I_g)*S_s + phi_art_s*A_s + theta*I_g
  dA_g = -delta*A_g - mu2*A_g - phi_art_g*A_g + tau_art_g*I_g
  dA_s = -delta*A_s - mu2*A_s - phi_art_s*A_s + tau_art_s*I_s
  dD_g = mu1*(I_g) + mu2*(A_g)
  dD_s = mu1*(I_s) + mu2*(A_s)
  ## return system of differential equation
  return(list(c(dS_g, dI_g, dA_g, dPrEP_g,  dS_s, dI_s, dA_s, dPrEP_s, dD_g, dD_s)))
}

## Assume frequency dependent because R0 does not change based on population size

run_simp_coupled_HIV_function <- function(beta, c_GS, c_GG, c_SG, sigma_GS, sigma_GG, sigma_SG, alpha, theta, delta, mu1, mu2, tau_prep_g, tau_prep_s, phi_prep_g, phi_prep_s, tau_art_g, tau_art_s, phi_art_g, phi_art_s, initial.state, max.time, freq.dependent) {
  # Determine beta divisor for each sub population
  beta_divisor_g <- ifelse(freq.dependent == TRUE, initial.state["S_g"] + initial.state["I_g"] + initial.state["A_g"] + initial.state["PrEP_g"], 1)
  beta_divisor_s <- ifelse(freq.dependent == TRUE, initial.state["S_s"] + initial.state["I_s"] + initial.state["A_s"] + initial.state["PrEP_s"], 1)
  
  # create param vector to solve the system of equations
  param <- c(beta = beta,
             c_GS = c_GS,
             c_GG = c_GG,
             c_SG = c_SG,
             sigma_GS = sigma_GS/beta_divisor_s,
             sigma_GG = sigma_GG/beta_divisor_g,
             sigma_SG = sigma_SG/beta_divisor_g,
             alpha = alpha, 
             theta = theta, 
             delta = delta, 
             mu1 = mu1, #AIDS_related death rate without ART
             mu2 = mu2, #AIDS_related death rate with ART
             tau_prep_g = tau_prep_g,
             tau_prep_s = tau_prep_s,
             phi_prep_g = phi_prep_g,
             phi_prep_s = phi_prep_s,
             tau_art_g = tau_art_g,
             tau_art_s = tau_art_s,
             phi_art_g = phi_art_g,
             phi_art_s = phi_art_s)
  times <- seq(0, max.time, 1)
  
  #solve the system of equations
  simplified_coupled_HIV_output <- deSolve::lsoda(initial.state, times, simplified_coupled_HIV_func, param)
  
  return(as.data.frame(simplified_coupled_HIV_output))
}
initial.state <- c("S_g" = 70000000, "I_g" = 56000, "A_g" = 344000, "PrEP_g" = 12768,
                   "S_s" = 480000, "I_s" = 16800, "A_s" = 103200, "PrEP_s" = 38304, 
                   "D_g" = 0, "D_s" = 0)

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
