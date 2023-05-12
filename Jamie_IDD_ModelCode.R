
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
  dS_g = -(betaGG*I_g + betaSG*I_s)*S_g - tau_prep_g*S_g - delta*S_g - theta*S_g + alpha*N + phi_prep_g*PrEP_g
  dS_s = -(betaGS*I_g)*S_s - tau_prep_s*S_s - delta*S_s + theta*S_g + phi_prep_s*PrEP_s
  dPrEP_g = -(phi_prep_g)*PrEP_g + tau_prep_g*S_g - delta*PrEP_g
  dPrEP_s = -(phi_prep_s)*PrEP_s + tau_prep_s*S_s - delta*PrEP_s
  dI_g = -delta*I_g - mu1*I_g - tau_art_g*I_g + (beta*c_GG*sigma_GG*I_g + 0.5*beta*c_SG*sigma_SG*I_s)*S_g + phi_art_g*A_g - theta*I_g
  dI_s = -delta*I_s - mu1*I_s - tau_art_s*I_s + (beta*c_GS*sigma_GS*I_g)*S_s + phi_art_s*A_s + theta*I_g
  dA_g = -delta*A_g - mu2*A_g - phi_art_g*A_g + tau_art_g*I_g
  dA_s = -delta*A_s - mu2*A_s - phi_art_s*A_s + tau_art_s*I_s
  dD = mu1*(I_g + I_s) + mu2*(A_g + A_s)
  ## return system of differential equation
  return(list(c(dS_g, dI_g, dA_g, dPrEP_g,  dS_s, dI_s, dA_s, dPrEP_s, dD)))
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
                   "S_s" = 480000, "I_s" = 16800, "A_s" = 103200, "PrEP_s" = 38304, "D" = 0)
# What actually happens - as is - data from before the policy was applied
scenario1 <- run_simp_coupled_HIV_function(beta = 0.0004,
                                           c_GS = 17.7,
                                           c_GG = 1,
                                           c_SG = 1.25,
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
                                           max.time = 10,
                                           freq.dependent = TRUE)

# Scenario with lower ART uptake and retention - Projection scenario for policy change
scenario2 <- run_simp_coupled_HIV_function(beta = 0.0004,
                                           c_GS = 17.7,
                                           c_GG = 1,
                                           c_SG = 1.25,
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
                                           tau_art_g = 0.712,
                                           phi_art_g = 0.015,
                                           tau_art_s = 0.712,
                                           phi_art_s = 0.015,
                                           initial.state = initial.state,
                                           max.time = 10,
                                           freq.dependent = TRUE)

# 
scenario1 <- run_simp_coupled_HIV_function(beta = 0.0004,
                                           c_GS = 17.7,
                                           c_GG = 1,
                                           c_SG = 1.25,
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
                                           tau_art = 0.712,
                                           phi_art = 0.015,
                                           initial.state = initial.state,
                                           max.time = 10,
                                           freq.dependent = TRUE)

## Things to analyze:
# Scenarios should have one that has differential intervention between two populations
# Deaths averted between scenarios
# Prevalence in the two populations
# Incidence 