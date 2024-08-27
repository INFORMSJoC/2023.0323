## Set parameters
X_t = 0.5 # current state [X_t] or [X_0]
mu = 0.8
delta = 1.1
sigma = 1.3
t_target = 1 # target time [t] for simulating point process

N_sim = 10^4 # number of simulation paths




## Function for calculating ([delta_tilde], [mu_tilde], [sigma_tilde], [kappa_2], [kappa_1], [kappa_0])
f_delta_mu_sigma_tilde_kappas = function(mu, delta, sigma, w=1){
  delta_tilde = sqrt(delta^2 + 2*sigma^2*w)
  mu_tilde = mu * delta^2 / delta_tilde^2
  sigma_tilde = sigma
  
  kappa_2 = (delta_tilde-delta) / (2*sigma_tilde^2)
  kappa_1 = 2*mu*delta*kappa_2 / delta_tilde
  kappa_0 = mu*delta*kappa_1 - sigma^2/2 * (kappa_1^2-2*kappa_2)
  
  return(list(delta_tilde=delta_tilde, mu_tilde=mu_tilde, sigma_tilde=sigma_tilde,
              kappa_2=kappa_2, kappa_1=kappa_1, kappa_0=kappa_0))
}




## Function for simulating interarrival time [tau] (Algorithm 4.1) and calculate its PDF and CDF for plotting
f_simulate_tau = function(X_t, mu, delta, sigma, N_sim=1, decompose=TRUE){
  # Calculate parameters ([delta_tilde], [mu_tilde], [sigma_tilde], [kappa_2], [kappa_1], [kappa_0])
  Parameters_tilde_kappas = f_delta_mu_sigma_tilde_kappas(mu, delta, sigma)
  delta_tilde = Parameters_tilde_kappas$delta_tilde
  mu_tilde = Parameters_tilde_kappas$mu_tilde
  sigma_tilde = Parameters_tilde_kappas$sigma_tilde
  kappa_2 = Parameters_tilde_kappas$kappa_2
  kappa_1 = Parameters_tilde_kappas$kappa_1
  kappa_0 = Parameters_tilde_kappas$kappa_0
  
  time_start = Sys.time() # record time
  
  # Function for calculating PDF of [V_star]
  f_PDF_V_star = function(s){
    varsigma_tilde = sigma_tilde * sqrt((1-exp(-2*delta_tilde*s)) / (2*delta_tilde))
    
    P_0 = mu_tilde^2*(1-exp(-delta_tilde*s)) / (2*delta^2 * (1-2*varsigma_tilde^2*kappa_2))
    P_0 = exp(P_0 * ((delta_tilde+2*delta) + (delta_tilde-2*delta) * exp(-delta_tilde*s)))
    P_0 = P_0 * exp(-kappa_0*s) / sqrt(1-2*varsigma_tilde^2*kappa_2)
    
    PDF_V_star = mu_tilde^2*(1-exp(-delta_tilde*s))^2 / (4*delta^2*(1-2*varsigma_tilde^2*kappa_2)^2)
    PDF_V_star = PDF_V_star * (delta_tilde*(1+exp(-delta_tilde*s)) + delta*(1-exp(-delta_tilde*s)))^2
    PDF_V_star = P_0 * (varsigma_tilde^2 / (1-2*varsigma_tilde^2*kappa_2) + PDF_V_star)
    
    return(PDF_V_star)
  }
  
  # Function for maximising [R_V_star(lambda, s)] in Algorithm 4.2 and return its maximum
  f_max_R_V_star = function(lambda){
    f_R_V_star = function(s){
      R_V_star = f_PDF_V_star(s) / dexp(s, lambda)
      
      return(R_V_star)
    }
    
    R_V_max = optimize(f_R_V_star, lower = 0.0001, upper =20, maximum = TRUE)$objective
    
    return(R_V_max)
  }
  
  # Calculate [lambda_star] and [C_V_Star] in Algorithm 4.2
  lambda_star = optimize(f_max_R_V_star, lower = 0.0001, upper =20, maximum = FALSE)$minimum
  C_V_star = optimize(f_max_R_V_star, lower = 0.0001, upper =20, maximum = FALSE)$objective
  
  # Function for simulating [V_star] (Algorithm 4.2)
  f_simulate_V_star = function(){
    repeat{
      V_hat = rexp(1, lambda_star)
      R_V_hat = f_PDF_V_star(V_hat) / dexp(V_hat, lambda_star)
      AR_V_star = R_V_hat / C_V_star
      
      if (runif(1, 0, 1) < AR_V_star){
        V_star = V_hat 
        break
      }
    }
    
    return(V_star)
  }
  
  # Function for simulating [V_t_star] (Algorithm 4.3)
  f_simulate_V_t_star = function(X_t){
    U_Z = runif(1, 0, 1)
    if (U_Z <= exp(-kappa_1*X_t - kappa_2*X_t^2)){
      V_t_star = Inf
    }
    else{
      Z = log(U_Z) + kappa_1*X_t + kappa_2*X_t^2
      a = kappa_2*X_t*(X_t-2*mu_tilde) - sigma_tilde^2*kappa_2/delta_tilde*Z
      b = 2*mu_tilde*X_t / delta
      c = -Z / (2*delta_tilde*kappa_2)
      V_t_star = -1/delta_tilde * log((-b + sqrt(b^2-4*a*c)) / (2*a))
    }
    
    return(V_t_star)
  }
  
  # Function for calculating PDF and CDF of interarrival time [tau]
  f_PDF_CDF_tau = function(X_t, s){
    nu_tilde = X_t*exp(-delta_tilde*s) + mu_tilde*(1-exp(-delta_tilde*s))
    varsigma_tilde = sigma_tilde * sqrt((1-exp(-2*delta_tilde*s)) / (2*delta_tilde))
    
    P = exp(-kappa_1*X_t - kappa_2*X_t^2 - kappa_0*s) / sqrt(1-2*varsigma_tilde^2*kappa_2)
    
    Q = varsigma_tilde^2*kappa_1^2 + 2*nu_tilde*kappa_1 + 2*nu_tilde^2*kappa_2
    Q = exp(Q / (2*(1-2*varsigma_tilde^2*kappa_2)))
    
    p = -(mu*mu_tilde + varsigma_tilde^2/(1-2*varsigma_tilde^2*kappa_2))
    
    q = (nu_tilde-mu_tilde) * (1-2*varsigma_tilde^2*kappa_2)
    q = q - sigma_tilde^2*exp(-2*delta_tilde*s)/(2*delta_tilde) * (kappa_1+2*nu_tilde*kappa_2)
    q = -q * delta_tilde*(kappa_1+2*nu_tilde*kappa_2) / (1-2*varsigma_tilde^2*kappa_2)^2
    
    PDF_tau = -P*Q * (p+q)
    CDF_tau = 1 - P*Q
    
    return(list(PDF=PDF_tau, CDF=CDF_tau))
  }
  
  # Function for calculating [R_tau(s)] in Algorithm 4.4
  f_R_tau = function(X_t, s){
    R_tau = f_PDF_CDF_tau(X_t, s)$PDF / dexp(s, kappa_0)
    
    return(R_tau)
  }
  
  # Function for maximising [R_tau(s)] in Algorithm 4.4
  f_max_R_tau = function(X_t){
    s_star = optimize(function(s) f_R_tau(X_t, s), lower = 0.0001, upper =20, maximum = TRUE)$mimimun
    C_M = optimize(function(s) f_R_tau(X_t, s), lower = 0.0001, upper =20, maximum = TRUE)$objective
    
    C_0 = (X_t-mu_tilde) - sigma_tilde^2/(2*delta_tilde)*(kappa_1+2*kappa_2*X_t)
    C_0 = C_0 * delta_tilde * (kappa_1+2*kappa_2*X_t)
    C_0 = (C_0 + mu*mu_tilde) / kappa_0
    
    C_inf = (delta_tilde+delta) * (delta_tilde+2*delta)
    C_inf = exp(C_inf * mu_tilde^2*delta_tilde*kappa_2^2/delta^2)
    C_inf = C_inf * sqrt(2*delta_tilde*kappa_2) * exp(-kappa_1*X_t - kappa_2*X_t^2)
    
    C_tau = max(C_0, C_M, C_inf)
    
    return(list(s_star=s_star, C_M=C_M, C_0=C_0, C_inf=C_inf, C_tau=C_tau))
  }
  
  # Function for simulating [tau] directly through its PDF (Algorithm 4.4)
  f_simulate_tau_direct = function(X_t){
    C_tau = f_max_R_tau(X_t)$C_tau
    repeat{
      tau_hat = rexp(1, kappa_0)
      R_tau_hat = f_PDF_CDF_tau(X_t, tau_hat)$PDF / dexp(tau_hat, kappa_0)
      AR_tau = R_tau_hat / C_tau
      
      if (is.nan(AR_tau)){
        next # avoid [AR_tau] to be NAN when [X_t] is too large and [tau_hat] is too small
      }
      else if (runif(1,0,1) < AR_tau){
        tau = tau_hat
        break
      }
    }
    
    return(tau)
  }
  
  # Define parameters ([a_star], [b_star], [c_star])
  a_star = -2*mu_tilde*sigma_tilde^2*kappa_2*X_t / delta
  b_star = X_t * (X_t-2*mu_tilde)
  c_star = mu_tilde*X_t / (delta*kappa_2)
  
  # Define conditions under which Algorithm 4.1 execute
  Condition_1 = kappa_1*X_t + kappa_2*X_t^2 > 0
  Condition_2a = (a_star > 0) & (b_star^2-4*a_star*c_star < 0)
  Condition_2b = (a_star > 0) & (-b_star/(2*a_star) > 1) & (a_star+b_star+c_star > 0)
  Condition_2c = (a_star > 0) & (-b_star/(2*a_star) < 0) & (c_star > 0)
  Condition_2d = (a_star < 0) & (c_star > 0) & (a_star+b_star+c_star > 0)
  
  time_end = Sys.time() # record time
  time_tau = difftime(time_end, time_start, units='sec')
  
  # Generate (a vector of) interarrival time [tau]
  count_decompose = 0
  count_direct = 0 # number of times of two ways to generate [tau]
  time_decompose = 0
  time_direct = 0 # running time of two ways to generate [tau]
  
  tau = rep(0, N_sim)
  for (i in 1:N_sim){
    if (decompose & Condition_1 & (Condition_2a | Condition_2b | Condition_2c | Condition_2d)){
      count_decompose = count_decompose + 1
      time_start = Sys.time() # record time
      
      # Decompose CDF of [tau] into CDFs of [V_star] and [V_t_star]
      V_star = f_simulate_V_star() # generate [V_star] (Algorithm 4.2)
      V_t_star = f_simulate_V_t_star(X_t) # generate [V_t_star] (Algorithm 4.3)
      
      tau[i] = min(V_star, V_t_star)
      
      time_end = Sys.time() # record time
      time_decompose = time_decompose + difftime(time_end, time_start, units='sec')
    }
    else{
      count_direct = count_direct + 1
      time_start = Sys.time() # record time
      
      # Generate [tau] directly through its PDF (Algorithm 4.4)
      tau[i] = f_simulate_tau_direct(X_t)
      
      time_end = Sys.time() # record time
      time_direct = time_direct + difftime(time_end, time_start, units='sec')
    }
  }
  
  time_tau = time_tau + time_decompose + time_direct
  
  # Calculate PDF and CDF of interarrival time [tau] for plotting
  tau_plot = seq(0, 10, 0.01)
  PDF_plot = rep(0, length(tau_plot))
  CDF_plot = rep(0, length(tau_plot))
  for (i in 1:length(tau_plot)){
    PDF_plot[i] = f_PDF_CDF_tau(X_t, tau_plot[i])$PDF
    CDF_plot[i] = f_PDF_CDF_tau(X_t, tau_plot[i])$CDF
  }
  
  # Calculate [R_tau(s)] for plotting
  C_M = f_max_R_tau(X_t)$C_M
  C_0 = f_max_R_tau(X_t)$C_0
  C_inf = f_max_R_tau(X_t)$C_inf
  C_tau = f_max_R_tau(X_t)$C_tau
  
  s_plot = seq(0, 5, 0.01)
  R_tau_plot = rep(0, length(s_plot))
  for (i in 1:length(s_plot)){
    R_tau_plot[i] = f_R_tau(X_t, s_plot[i])
  }
  
  return(list(tau=tau, tau_plot=tau_plot, PDF_plot=PDF_plot, CDF_plot=CDF_plot,
              C_M=C_M, C_0=C_0, C_inf=C_inf, C_tau=C_tau, s_plot=s_plot, R_tau_plot=R_tau_plot,
              count_decompose=count_decompose, count_direct=count_direct,
              time_tau=time_tau, time_decompose=time_decompose, time_direct=time_direct))
}

# Plot sample density and ECDF of simulated interarrival time [tau] and its actual PDF and CDF
tau_simulate = f_simulate_tau(X_t, mu, delta, sigma, N_sim, decompose=FALSE)
tau = tau_simulate$tau
tau_plot = tau_simulate$tau_plot
PDF_tau_plot = tau_simulate$PDF_plot
CDF_tau_plot = tau_simulate$CDF_plot

plot(tau_plot, PDF_tau_plot, type="l", col="red")
lines(density(tau), type='l', col="blue")

plot(tau_plot, CDF_tau_plot, type="l", col="red")
lines(ecdf(tau), col="blue")

# Plot [R_tau(s)] in Algorithm 4.4
C_M = tau_simulate$C_M
C_0 = tau_simulate$C_0
C_inf = tau_simulate$C_inf

s_plot = tau_simulate$s_plot
R_tau_plot = tau_simulate$R_tau_plot

plot(s_plot, R_tau_plot, type="l", col="black")
lines(c(0, 5), c(C_M, C_M), type="l", lty="dotted", col="black")
lines(c(0, 5), c(C_0, C_0), type="l", lty="dotdash", col="black")
lines(c(0, 5), c(C_inf, C_inf), type="l", lty="dashed", col="black")




## Function for simulating pre-event intensity [lambda_pre] and state [X_pre] (Algorithm 4.5)
f_simulate_X_pre = function(X_t, tau, mu, delta, sigma, N_sim=1){
  # Calculate parameters ([delta_tilde], [mu_tilde], [sigma_tilde], [kappa_2], [kappa_1])
  Parameters_tilde_kappas = f_delta_mu_sigma_tilde_kappas(mu, delta, sigma)
  delta_tilde = Parameters_tilde_kappas$delta_tilde
  mu_tilde = Parameters_tilde_kappas$mu_tilde
  sigma_tilde = Parameters_tilde_kappas$sigma_tilde
  kappa_2 = Parameters_tilde_kappas$kappa_2
  kappa_1 = Parameters_tilde_kappas$kappa_1
  
  time_start = Sys.time() # record time
  
  # Define parameters ([nu_tilde], [varsigma_tilde])
  nu_tilde_tau = X_t*exp(-delta_tilde*tau) + mu_tilde*(1-exp(-delta_tilde*tau))
  varsigma_tilde_tau = sigma_tilde * sqrt((1-exp(-2*delta_tilde*tau)) / (2*delta_tilde))
  
  # Define parameters ([nu_tilde], [varsigma_tilde], [A], [B], [C], [p_lambda])
  A = varsigma_tilde_tau^2/2 * (kappa_1 + nu_tilde_tau/varsigma_tilde_tau^2)^2
  B = 2 * varsigma_tilde_tau^2
  C = 1 - 2 * varsigma_tilde_tau^2 * kappa_2
  p_lambda = (A/C) / (A/C + 1/2)
  
  # Define parameters ([nu], [varsigma])
  nu_tau = X_t*exp(-delta*tau) + mu*(1-exp(-delta*tau))
  varsigma_tau = sigma * sqrt((1-exp(-2*delta*tau)) / (2*delta))
  
  # Generate pre-event intensity [lambda_pre]
  B_lambda = rbinom(N_sim, 1, p_lambda)
  J = rpois(N_sim, A/C)
  lambda_pre = rgamma(N_sim, J+3/2+B_lambda, C/B)
  # lambda_pre = rep(0, N_sim)
  # for (i in 1:N_sim){
  #   B_lambda = rbinom(1, 1, p_lambda)
  #   J = rpois(1, A/C)
  #   lambda_pre[i] = rgamma(1, J+3/2+B_lambda, C/B)
  # }
  
  time_end_1 = Sys.time() # record time
  time_lambda_pre = difftime(time_end_1, time_start, units='sec')
  
  # Generate pre-event state [X_pre]
  X = sqrt(lambda_pre)
  p_X = exp(-(X-nu_tau)^2/(2*varsigma_tau^2)) + exp(-(X+nu_tau)^2/(2*varsigma_tau^2))
  p_X = exp(-(X-nu_tau)^2/(2*varsigma_tau^2)) / p_X
  
  B_X = rbinom(N_sim, 1, p_X)
  X_pre = (2*B_X - 1) * X
  
  time_end_2 = Sys.time() # record time
  time_X_pre = difftime(time_end_2, time_start, units='sec')
  
  return(list(lambda_pre=lambda_pre, X_pre=X_pre, time_lambda=time_lambda_pre, time_X=time_X_pre))
}

## Function for calculating Laplace transform of pre-event intensity [lambda_pre]
f_Laplace_lambda_pre = function(u, X_t, tau, mu, delta, sigma){
  # Calculate parameters ([delta_tilde], [mu_tilde], [sigma_tilde], [kappa_2], [kappa_1])
  Parameters_tilde_kappas = f_delta_mu_sigma_tilde_kappas(mu, delta, sigma)
  delta_tilde = Parameters_tilde_kappas$delta_tilde
  mu_tilde = Parameters_tilde_kappas$mu_tilde
  sigma_tilde = Parameters_tilde_kappas$sigma_tilde
  kappa_2 = Parameters_tilde_kappas$kappa_2
  kappa_1 = Parameters_tilde_kappas$kappa_1
  
  # Define parameters ([nu_tilde], [varsigma_tilde])
  nu_tilde_tau = X_t*exp(-delta_tilde*tau) + mu_tilde*(1-exp(-delta_tilde*tau))
  varsigma_tilde_tau = sigma_tilde * sqrt((1-exp(-2*delta_tilde*tau)) / (2*delta_tilde))
  
  # Define parameters ([nu_tilde], [varsigma_tilde], [A], [B], [C], [p_lambda])
  A = varsigma_tilde_tau^2/2 * (kappa_1 + nu_tilde_tau/varsigma_tilde_tau^2)^2
  B = 2 * varsigma_tilde_tau^2
  C = 1 - 2 * varsigma_tilde_tau^2 * kappa_2
  p_lambda = (A/C) / (A/C + 1/2)
  
  # Define Laplace transform of pre-event intensity [lambda_pre]
  f_Psi_derivative = function(u){
    -B/(B*u+C)^(3/2) * (A/(B*u+C) + 1/2) * exp(A/(B*u+C))
  }
  Laplace_lambda_pre = f_Psi_derivative(u) / f_Psi_derivative(0)
  
  return(Laplace_lambda_pre)
}

## Function for simulating pre-event intensity [lambda_pre] by numerical Laplace inversion (Ridout, 2009)
f_simulate_lambda_pre_Ridout = function(X_t, tau, mu, delta, sigma, N_sim=1){
  source("rlaptrans.r")
  
  time_start = Sys.time() # record time
  
  # Generate pre-event intensity [lambda_pre] by Ridout's method
  lambda_pre = rlaptrans(N_sim, function(u) f_Laplace_lambda_pre(u, X_t, tau, mu, delta, sigma))
  
  time_end = Sys.time() # record time
  time_Ridout = difftime(time_end, time_start, units='sec')
  
  return(list(lambda_pre=lambda_pre, time=time_Ridout))
}

## Function for calculating [f(t)], inverse Laplace inversion of [f_hat(s)], by Euler Algorithm (Abate and Whitt, 2006)
f_inverse_Laplace_Euler = function(f_hat, t, M=11){
  k = 1 : (2*M)
  beta = M/3*log(10) + 1i*pi*c(0, k)
  xi = pbinom(2*M-k, size=M, prob=0.5)
  eta = c(1/2, (-1)^k * xi)
  
  f_t_Euler = 10^(M/3) / t * sum(eta * Re(f_hat(beta/t)))
  
  return(f_t_Euler)
}

## Function for calculating PDF, CDF and CCDF of pre-event intensity [lambda_pre] by inverting its Laplace transform
f_PDF_CDF_CCDF_lambda_pre = function(lambda, X_t, tau, mu, delta, sigma){
  f_Laplace_lambda_pre_u = function(u) f_Laplace_lambda_pre(u, X_t, tau, mu, delta, sigma)
  
  # Calculate PDF, CDF and CCDF of [lambda_pre] by inverting its functions of Laplace transform (Abate et al., 2000)
  PDF_lambda_pre = f_inverse_Laplace_Euler(f_Laplace_lambda_pre_u, lambda)
  CDF_lambda_pre = f_inverse_Laplace_Euler(function(u) f_Laplace_lambda_pre_u(u)/u, lambda)
  CCDF_lambda_pre = f_inverse_Laplace_Euler(function(u) (1-f_Laplace_lambda_pre_u(u))/u, lambda) # complementary CDF, i.e., 1-CDF
  
  return(list(PDF=PDF_lambda_pre, CDF=CDF_lambda_pre, CCDF=CCDF_lambda_pre))
}

# Compare running time of our exact simulation and Ridout's simulation method
tau = f_simulate_tau(X_t, mu, delta, sigma)$tau
lambda_pre_simulate_exact = f_simulate_X_pre(X_t, tau, mu, delta, sigma, N_sim)
lambda_pre_simulate_Ridout = f_simulate_lambda_pre_Ridout(X_t, tau, mu, delta, sigma, N_sim)

lambda_pre_simulate_exact$time_lambda # running time of our exact simulation
lambda_pre_simulate_Ridout$time # running time of Ridout's simulation

# Plot sample density and ECDF of pre-event intensity [lambda_pre] and its actual PDF and CDF
lambda_pre_exact = lambda_pre_simulate_exact$lambda_pre
lambda_pre_Ridout = lambda_pre_simulate_Ridout$lambda_pre

lambda_pre_plot = seq(0, 10, 0.01)
PDF_lambda_pre_plot = rep(0, length(lambda_pre_plot))
CDF_lambda_pre_plot = rep(0, length(lambda_pre_plot))
for (i in 1:length(lambda_pre_plot)){
  PDF_CDF_CCDF_lambda_pre_plot = 
    f_PDF_CDF_CCDF_lambda_pre(lambda_pre_plot[i], X_t, tau, mu, delta, sigma)
  PDF_lambda_pre_plot[i] = PDF_CDF_CCDF_lambda_pre_plot$PDF
  if (lambda_pre_plot[i] < 5){
    CDF_lambda_pre_plot[i] = PDF_CDF_CCDF_lambda_pre_plot$CDF
  }
  else{
    # Guarantee accuracy in upper tail of distribution
    CDF_lambda_pre_plot[i] = 1 - PDF_CDF_CCDF_lambda_pre_plot$CCDF
  }
}

plot(lambda_pre_plot, PDF_lambda_pre_plot, type = "l", col="red")
lines(density(lambda_pre_exact), type='l', col="blue")
lines(density(lambda_pre_Ridout), type='l', col="green")

plot(lambda_pre_plot, CDF_lambda_pre_plot, type = "l", col="red")
lines(ecdf(lambda_pre_exact), col="blue")
lines(ecdf(lambda_pre_Ridout), col="green")




## Function for simulating point process [N_t] with quadratic OU intensity (Algorithm 4.6) (with/without self-exciting jump in state [X_t])
f_simulate_N = function(X_0, mu, delta, sigma, t_target, N_sim=1, decompose=TRUE, jump="0", 
                        loss="0.5"){
  # Calculate parameters ([delta_tilde], [mu_tilde], [sigma_tilde], [kappa_2], [kappa_1], [kappa_0])
  Parameters_tilde_kappas = f_delta_mu_sigma_tilde_kappas(mu, delta, sigma)
  delta_tilde = Parameters_tilde_kappas$delta_tilde
  mu_tilde = Parameters_tilde_kappas$mu_tilde
  sigma_tilde = Parameters_tilde_kappas$sigma_tilde
  kappa_2 = Parameters_tilde_kappas$kappa_2
  kappa_1 = Parameters_tilde_kappas$kappa_1
  kappa_0 = Parameters_tilde_kappas$kappa_0
  
  time_start = Sys.time() # record time
  
  # Function for calculating PDF of [V_star]
  f_PDF_V_star = function(s){
    varsigma_tilde = sigma_tilde * sqrt((1-exp(-2*delta_tilde*s)) / (2*delta_tilde))
    
    P_0 = mu_tilde^2*(1-exp(-delta_tilde*s)) / (2*delta^2 * (1-2*varsigma_tilde^2*kappa_2))
    P_0 = exp(P_0 * ((delta_tilde+2*delta) + (delta_tilde-2*delta) * exp(-delta_tilde*s)))
    P_0 = P_0 * exp(-kappa_0*s) / sqrt(1-2*varsigma_tilde^2*kappa_2)
    
    PDF_V_star = mu_tilde^2*(1-exp(-delta_tilde*s))^2 / (4*delta^2*(1-2*varsigma_tilde^2*kappa_2)^2)
    PDF_V_star = PDF_V_star * (delta_tilde*(1+exp(-delta_tilde*s)) + delta*(1-exp(-delta_tilde*s)))^2
    PDF_V_star = P_0 * (varsigma_tilde^2 / (1-2*varsigma_tilde^2*kappa_2) + PDF_V_star)
    
    return(PDF_V_star)
  }
  
  # Function for maximising [R_V_star(lambda, s)] in Algorithm 4.2 and return its maximum
  f_max_R_V_star = function(lambda){
    f_R_V_star = function(s){
      R_V_star = f_PDF_V_star(s) / dexp(s, lambda)
      
      return(R_V_star)
    }
    
    R_V_max = optimize(f_R_V_star, lower = 0.0001, upper =20, maximum = TRUE)$objective
    
    return(R_V_max)
  }
  
  # Calculate [lambda_star] and [C_V_Star] in Algorithm 4.2
  lambda_star = optimize(f_max_R_V_star, lower = 0.0001, upper =20, maximum = FALSE)$minimum
  C_V_star = optimize(f_max_R_V_star, lower = 0.0001, upper =20, maximum = FALSE)$objective
  
  # Function for simulating [V_star] (Algorithm 4.2)
  f_simulate_V_star = function(){
    repeat{
      V_hat = rexp(1, lambda_star)
      R_V_hat = f_PDF_V_star(V_hat) / dexp(V_hat, lambda_star)
      AR_V_star = R_V_hat / C_V_star
      
      if (runif(1, 0, 1) < AR_V_star){
        V_star = V_hat 
        break
      }
    }
    
    return(V_star)
  }
  
  # Function for simulating [V_t_star] (Algorithm 4.3)
  f_simulate_V_t_star = function(X_t){
    U_Z = runif(1, 0, 1)
    if (U_Z <= exp(-kappa_1*X_t - kappa_2*X_t^2)){
      V_t_star = Inf
    }
    else{
      Z = log(U_Z) + kappa_1*X_t + kappa_2*X_t^2
      a = kappa_2*X_t*(X_t-2*mu_tilde) - sigma_tilde^2*kappa_2/delta_tilde*Z
      b = 2*mu_tilde*X_t / delta
      c = -Z / (2*delta_tilde*kappa_2)
      V_t_star = -1/delta_tilde * log((-b + sqrt(b^2-4*a*c)) / (2*a))
    }
    
    return(V_t_star)
  }
  
  # Function for calculating PDF and CDF of interarrival time [tau]
  f_PDF_CDF_tau = function(X_t, s){
    nu_tilde = X_t*exp(-delta_tilde*s) + mu_tilde*(1-exp(-delta_tilde*s))
    varsigma_tilde = sigma_tilde * sqrt((1-exp(-2*delta_tilde*s)) / (2*delta_tilde))
    
    P = exp(-kappa_1*X_t - kappa_2*X_t^2 - kappa_0*s) / sqrt(1-2*varsigma_tilde^2*kappa_2)
    
    Q = varsigma_tilde^2*kappa_1^2 + 2*nu_tilde*kappa_1 + 2*nu_tilde^2*kappa_2
    Q = exp(Q / (2*(1-2*varsigma_tilde^2*kappa_2)))
    
    p = -(mu*mu_tilde + varsigma_tilde^2/(1-2*varsigma_tilde^2*kappa_2))
    
    q = (nu_tilde-mu_tilde) * (1-2*varsigma_tilde^2*kappa_2)
    q = q - sigma_tilde^2*exp(-2*delta_tilde*s)/(2*delta_tilde) * (kappa_1+2*nu_tilde*kappa_2)
    q = -q * delta_tilde*(kappa_1+2*nu_tilde*kappa_2) / (1-2*varsigma_tilde^2*kappa_2)^2
    
    PDF_tau = -P*Q * (p+q)
    CDF_tau = 1 - P*Q
    
    return(list(PDF=PDF_tau, CDF=CDF_tau))
  }
  
  # Function for calculating [R_tau(s)] in Algorithm 4.4
  f_R_tau = function(X_t, s){
    R_tau = f_PDF_CDF_tau(X_t, s)$PDF / dexp(s, kappa_0)
    
    return(R_tau)
  }
  
  # Function for maximising [R_tau(s)] in Algorithm 4.4
  f_max_R_tau = function(X_t){
    s_star = optimize(function(s) f_R_tau(X_t, s), lower = 0.0001, upper =20, maximum = TRUE)$mimimun
    C_M = optimize(function(s) f_R_tau(X_t, s), lower = 0.0001, upper =20, maximum = TRUE)$objective
    
    C_0 = (X_t-mu_tilde) - sigma_tilde^2/(2*delta_tilde)*(kappa_1+2*kappa_2*X_t)
    C_0 = C_0 * delta_tilde * (kappa_1+2*kappa_2*X_t)
    C_0 = (C_0 + mu*mu_tilde) / kappa_0
    
    C_inf = (delta_tilde+delta) * (delta_tilde+2*delta)
    C_inf = exp(C_inf * mu_tilde^2*delta_tilde*kappa_2^2/delta^2)
    C_inf = C_inf * sqrt(2*delta_tilde*kappa_2) * exp(-kappa_1*X_t - kappa_2*X_t^2)
    
    C_tau = max(C_0, C_M, C_inf)
    
    return(list(s_star=s_star, C_M=C_M, C_0=C_0, C_inf=C_inf, C_tau=C_tau))
  }
  
  # Function for simulating [tau] directly through its PDF (Algorithm 4.4)
  f_simulate_tau_direct = function(X_t){
    C_tau = f_max_R_tau(X_t)$C_tau
    repeat{
      tau_hat = rexp(1, kappa_0)
      R_tau_hat = f_PDF_CDF_tau(X_t, tau_hat)$PDF / dexp(tau_hat, kappa_0)
      AR_tau = R_tau_hat / C_tau
      
      if (is.nan(AR_tau)){
        next # avoid [AR_tau] to be NAN when [X_t] is too large and [tau_hat] is too small
      }
      else if (runif(1,0,1) < AR_tau){
        tau = tau_hat
        break
      }
    }
    
    return(tau)
  }
  
  count_decompose = 0 # number of times using Algorithm 4.2 and 4.3
  count_direct = 0 # number of times using Algorithm 4.4
  
  N_t = rep(0, N_sim)
  L_t = rep(0, N_sim) # total loss
  T_event = vector("list", length=N_sim) # all arrival time in every path
  X_event = vector("list", length=N_sim) # state [X_t] at all arrival time in every path
  lambda_event = vector("list", length=N_sim) # intensity [lambda_t] at all arrival time in every path
  for (i in 1:N_sim){
    T_current = 0 # current time
    N_T = 0
    X_T = X_0
    L_T = 0
    
    T_arrival = c() # record all arrival time
    X_arrival = c() # record state [X_t] at all arrival time
    lambda_arrival = c() # record intensity [lambda_t] at all arrival time
    
    repeat{
      ## Generate interarrival time [tau] (Algorithm 4.1)
      # Define parameters ([a_star], [b_star], [c_star])
      a_star = -2*mu_tilde*sigma_tilde^2*kappa_2*X_T / delta
      b_star = X_T * (X_T-2*mu_tilde)
      c_star = mu_tilde*X_T / (delta*kappa_2)
      
      # Define conditions under which Algorithm 4.1 execute
      Condition_1 = kappa_1*X_T + kappa_2*X_T^2 > 0
      Condition_2a = (a_star > 0) & (b_star^2-4*a_star*c_star < 0)
      Condition_2b = (a_star > 0) & (-b_star/(2*a_star) > 1) & (a_star+b_star+c_star > 0)
      Condition_2c = (a_star > 0) & (-b_star/(2*a_star) < 0) & (c_star > 0)
      Condition_2d = (a_star < 0) & (c_star > 0) & (a_star+b_star+c_star > 0)
      
      if (decompose & Condition_1 & (Condition_2a | Condition_2b | Condition_2c | Condition_2d)){
        count_decompose = count_decompose + 1
        
        # Decompose CDF of [tau] into CDFs of [V_star] and [V_t_star]
        V_star = f_simulate_V_star() # generate [V_star] (Algorithm 4.2)
        V_t_star = f_simulate_V_t_star(X_T) # generate [V_t_star] (Algorithm 4.3)
        
        tau = min(V_star, V_t_star)
      }
      else{
        count_direct = count_direct + 1
        
        # Generate [tau] directly through its PDF (Algorithm 4.4)
        tau = f_simulate_tau_direct(X_T)
      }
      
      ## Generate pre-event state position [X_pre] (Algorithm 4.5)
      X_T = f_simulate_X_pre(X_T, tau, mu, delta, sigma)$X_pre
      
      ## Update process state
      T_current = T_current + tau
      if (T_current > t_target){
        break
      }
      else{
        N_T = N_T + 1
        L_T = L_T + eval(parse(text=loss))
        T_arrival = c(T_arrival, T_current)
        X_arrival = c(X_arrival, X_T)
        lambda_arrival = c(lambda_arrival, X_T^2)
      }
      
      ## Add self-exciting jump size
      X_T = X_T + eval(parse(text=jump))
      # X_arrival = c(X_arrival, X_T)
      # lambda_arrival = c(lambda_arrival, X_T^2)
    }
    
    N_t[i] = N_T
    L_t[i] = L_T
    T_event[[i]] = T_arrival
    X_event[[i]] = X_arrival
    lambda_event[[i]] = lambda_arrival
  }
  
  if (length(T_event) < N_sim){
    T_event = c(T_event, vector("list", length=N_sim-length(T_event)))
    X_event = c(X_event, vector("list", length=N_sim-length(T_event)))
    lambda_event = c(lambda_event, vector("list", length=N_sim-length(T_event)))
  }
  
  time_end = Sys.time() # record time
  time_N = difftime(time_end, time_start, units='sec')
  
  return(list(N=N_t, L=L_t, T_event=T_event, X_event=X_event, lambda_event=lambda_event,
              time=time_N, count_decompose=count_decompose, count_direct=count_direct))
}

## Function for calculating expectation of point process [N_t] with quadratic OU intensity
f_E_N = function(X_t, s, mu, delta, sigma, N_t=0){
  E_N = N_t + (mu^2 + sigma^2/(2*delta)) * s
  E_N = E_N + 2*mu/delta * (X_t-mu) * (1-exp(-delta*s))
  E_N = E_N + 1/(2*delta) * ((X_t-mu)^2 - sigma^2/(2*delta)) * (1-exp(-2*delta*s))
  
  return(E_N)
}

# Plot simulated path of OU state [X_t] and intensity [lambda_t]
N_simulate = f_simulate_N(X_t, mu, delta, sigma, t_target=500)
T_event = N_simulate$T_event
X_event = N_simulate$X_event
lambda_event = N_simulate$lambda_event

plot(c(0, T_event[[1]]), c(X_t, X_event[[1]]), type = "l", col="black")
plot(c(0, T_event[[1]]), c(X_t^2, lambda_event[[1]]), type = "l", col="black")

# Plot simulated path of point process [N_t] with quadratic OU intensity
N = N_simulate$N

plot(c(0, T_event[[1]]), 0:N, type = "l", col="black")
hist(T_event[[1]], breaks=50)

# Calculate sample mean error and RMSE of simulated point process [N_t] with quadratic OU intensity
E_N = f_E_N(X_t, t_target, mu, delta, sigma)
N_simulate = f_simulate_N(X_t, mu, delta, sigma, t_target, N_sim=15625, decompose=FALSE)
N = N_simulate$N
N_simulate$time
error_N = mean(N) - E_N
error_percent_N = error_N / E_N
RMSE = sd(N) / sqrt(N_sim)
error_N
error_percent_N
RMSE

# Plot simulated path of OU state [X_t] and intensity [lambda_t]
N_simulate = f_simulate_N(X_t, mu, delta, sigma, t_target=500, jump="rnorm(1, 0, 1)")
T_event = N_simulate$T_event
X_event = N_simulate$X_event
lambda_event = N_simulate$lambda_event

plot(c(0, T_event[[1]]), c(X_t, X_event[[1]]), type = "l", col="black")
plot(c(0, T_event[[1]]), c(X_t^2, lambda_event[[1]]), type = "l", col="black")

# Plot simulated path of point process [N_t] with quadratic OU intensity
N = N_simulate$N

plot(c(0, T_event[[1]]), 0:N, type = "l", col="black")
hist(T_event[[1]], breaks=50)

# Plot ECDF of simulated total loss [L_t]
L_t = f_simulate_N(X_t, mu, delta, sigma, t_target=5, N_sim, jump="rnorm(1, 0, 1)", loss="rexp(1, 2)")$L
plot(ecdf(L_t), xlim=c(0, 20), col="red")

# Calculate estimated PMF of point process [N_t] with quadratic OU intensity
N_simulate = f_simulate_N(X_t, mu, delta, sigma, t_target, N_sim, jump="rnorm(1, 0, 2)")
N = N_simulate$N
EPMF_N = table(N) / N_sim
EPMF_N
