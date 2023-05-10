alpha = 


simplified_coupled_HIV_func <- function(param, y, t){
  ## parameters
  # population
    alpha = param['alpha'] # birth rate - more specifically rate of people being introduced to the susceptible population once they turn 18
    theta = param['theta']# rate of starting sex work
    delta = # rate of non-AIDS-related mortality
  # transmission
    betaSG = param['betaSG']
    betaGG = param['betaGG']
    betaSS = param['betaSS']
  
  # death rates
    mu = param['mu'] # AIDS-related mortality rate
  # prep
    tau_prep_g = param['tau_prep_g']# rate of initiating prep in the general population
    tau_prep_s = param['tau_prep_s']# rate of initiating prep in the SWs population
    phi_prep_g = param['phi_prep_g']# rate of prep discontinuation - general population
    phi_prep_s = param['phi_prep_g']# rate of prep discontinuation - SWs population
  # art
    tau_art = param['tau_art']# rate of initiating ART
    phi_art = param['phi_art']# rate of ART discontinuation
## defining states
    S_g = y['S_g'], I_g = y['I_g'], A_g = y['A_g'], 
    S_s = y['S_s'], I_s = y['I_s'], A_s = y['A_s'],
    PrEP_g = y['PrEP_g'], PrEP_s = y['PrEP_s'], D = y['D']

## Equations for the general population
  dS_g = -(betaGG*I_g + betaSG*I_s)*S_g - tau_prep_g*S_g - delta*S_g - theta*S_g + alpha*N + phi_prep_g*PrEP
  dS_s = -(betaSG*I_g + betaSS*I_s)*S_s - tau_prep_s*S_s - delta*S_s + theta*S_g + phi_prep_s*PREP
  dPrEP = -(phi_prep_g + phi_prep_s)*PrEP + tau_prep_g*S_g + tau_prep_s*S_s
  dI_g = -delta*I_g - mu*I_g - tau_art*I_g + (betaGG*I_g + betaSG*I_s)*S_g + phi_art*A_g 
  dI_s = -delta*I_s - mu*I_s - tau_art*I_s + (betaSG*I_g + betaSS*I_s)*S_s + phi_art*A_s 
  dA_g = -delta*A_g - mu*A_g - phi_art*A_g + tau_art*I_g
  dA_s = -delta*A_s - mu*A_s - phi_art*A_s + tau_art*I_s
  dD = mu*(I_g + I_s + A_g + A_s)
  ## return system of differential equation
  return(list(c(dS_g, dS_s, dPrEP, dI_g, dI_s, dA_g, dA_s, dD)))
}

## Assume frequency dependent because R0 does not change based on population size
  
run_simp_coupled_HIV_function <- function(param, initial.state, max.time, freq.dependent) {
  beta_divisor_1 <- ifelse(freq.dependent == TRUE, initial.state["S1"] + initial.state["I1"] + initial.state["R1"], 1)
  beta_divisor_2 <- ifelse(freq.dependent == TRUE, initial.state["S2"] + initial.state["I2"] + initial.state["R2"], 1)
}