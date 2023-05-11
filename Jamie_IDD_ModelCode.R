
simplified_coupled_HIV_func <- function(param, y, t){
  ## parameters
  # population
    alpha = param['alpha'] # birth rate - more specifically rate of people being introduced to the susceptible population once they turn 18
    theta = param['theta']# rate of starting sex work
    delta = param['delta'] # rate of non-AIDS-related mortality
  # transmission
    betaSG = param['betaSG']
    betaGG = param['betaGG']
    betaSS = param['betaSS']
  
  # death rates
    mu1 = param['mu1'] # AIDS-related mortality rate _ without ART
    mu2 = param['mu2'] # AIDS-related mortality rate _ with ART
  # prep
    tau_prep_g = param['tau_prep_g']# rate of initiating prep in the general population
    tau_prep_s = param['tau_prep_s']# rate of initiating prep in the SWs population
    phi_prep_g = param['phi_prep_g']# rate of prep discontinuation - general population
    phi_prep_s = param['phi_prep_g']# rate of prep discontinuation - SWs population
  # art
    tau_art = param['tau_art']# rate of initiating ART
    phi_art = param['phi_art']# rate of ART discontinuation
## defining states
    S_g = y['S_g']; I_g = y['I_g']; A_g = y['A_g'] 
    S_s = y['S_s']; I_s = y['I_s']; A_s = y['A_s']
    PrEP_g = y['PrEP_g']; PrEP_s = y['PrEP_s']; D = y['D']
    N = S_g + I_g + A_g + S_s + I_s + A_s + PrEP_g + PrEP_s
## Equations for the general population
  dS_g = -(betaGG*I_g + betaSG*I_s)*S_g - tau_prep_g*S_g - delta*S_g - theta*S_g + alpha*N + phi_prep_g*PrEP_g
  dS_s = -(betaSG*I_g + betaSS*I_s)*S_s - tau_prep_s*S_s - delta*S_s + theta*S_g + phi_prep_s*PrEP_s
  dPrEP_g = -(phi_prep_g)*PrEP_g + tau_prep_g*S_g
  dPrEP_s = -(phi_prep_s)*PrEP_s + tau_prep_s*S_s
  dI_g = -delta*I_g - mu1*I_g - tau_art*I_g + (betaGG*I_g + betaSG*I_s)*S_g + phi_art*A_g 
  dI_s = -delta*I_s - mu1*I_s - tau_art*I_s + (betaSG*I_g + betaSS*I_s)*S_s + phi_art*A_s 
  dA_g = -delta*A_g - mu2*A_g - phi_art*A_g + tau_art*I_g
  dA_s = -delta*A_s - mu2*A_s - phi_art*A_s + tau_art*I_s
  dD = mu1*(I_g + I_s) + mu2*(A_g + A_s)
  ## return system of differential equation
  return(list(c(dS_g, dS_s, dPrEP_g, dPrEP_s, dI_g, dI_s, dA_g, dA_s, dD)))
}

## Assume frequency dependent because R0 does not change based on population size
  
run_simp_coupled_HIV_function <- function(betaSS, betaSG, betaGG, alpha, theta, delta, mu1, mu2, tau_prep_g, tau_prep_s, phi_prep_g, phi_prep_s, tau_art, phi_art, initial.state, max.time, freq.dependent) {
  # Determine beta divisor for each sub population
  beta_divisor_g <- ifelse(freq.dependent == TRUE, initial.state["S_g"] + initial.state["I_g"] + initial.state["A_g"] + initial.state["PrEP_g"], 1)
  beta_divisor_s <- ifelse(freq.dependent == TRUE, initial.state["S_s"] + initial.state["I_s"] + initial.state["A_s"] + initial.state["PrEP_s"], 1)
  beta_divisor_mixing <- ifelse(freq.dependent == TRUE, initial.state["S_s"] + initial.state["I_s"] + initial.state["A_s"] + initial.state["PrEP_s"] + initial.state["S_g"] + initial.state["I_g"] + initial.state["A_g"] + initial.state["PrEP_g"], 1)
  
  # create param vector to solve the system of equations
  param <- c(betaSS = betaSS/beta_divisor_s, 
             betaSG = betaSG/beta_divisor_mixing, 
             betaGG = betaGG/beta_divisor_g, 
             alpha = alpha, 
             theta = theta, 
             delta = delta, 
             mu1 = mu1, #AIDS_related death rate without ART
             mu2 = mu2, #AIDS_related death rate with ART
             tau_prep_g = tau_prep_g,
             tau_prep_s = tau_prep_s,
             phi_prep_g = phi_prep_g,
             phi_prep_s = phi_prep_s,
             tau_art = tau_art,
             phi_art = phi_art)
  times <- seq(0, max.time, 1)
  
  #solve the system of equations
  simplified_coupled_HIV_output <- deSolve::lsoda(initial.state, times, simplified_coupled_HIV_func, param)
  
  return(as.data.frame(simplified_coupled_HIV_output))
}
initial.state <- c("S_g" = 71000000, "I_g" = 1, "A_g" = 0, "PrEP_g" = 0,
                   "S_s" = 600000, "I_g" = 1, "A_s" = 0, "PrEP_s" = 0, "D" = 0)

scenario1 <- run_simp_coupled_HIV_function(betaSS = 0.002900, 
                                           betaSG = 0.00020, 
                                           betaGG = 0.0005,
                                           alpha = 0.009532, 
                                           theta = 0.00001, 
                                           delta = 0.007, 
                                           mu1 = 0.2149,
                                           mu2 = 0.04228,
                                           tau_prep_g = 0.0000214,
                                           tau_prep_s = 0.00583,
                                           phi_prep_g = 0.73,
                                           phi_prep_s = 0.73,
                                           tau_art = 0.712,
                                           phi_art = 0.247,
                                           initial.state,
                                           max.time = 10,
                                           freq.dependent = TRUE)
