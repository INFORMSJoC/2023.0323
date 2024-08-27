%% Function for calculating PMF of point process [N_(t+s)] with quadratic OU intensity given state [X_t]
function PMF_N = f_PMF_N(X_t, delta, mu, sigma, s, M_max)
% X_t -- (1*1 double) initial value of {X_t} on [t, t+s]
% delta -- (1*1 double)
% mu -- (1*1 double)
% sigma -- (1*1 double)
% s -- (1*1 double)
% M_max -- (1*1 double) largest [m] in [P(N=m|X_t)] to be calculate

    %% Define functions
    syms w

    kappa_2(w) = (-delta + sqrt(delta^2+2*sigma^2*w)) / (2*sigma^2);
    kappa_1(w) = (2*delta*mu*kappa_2(w)) / (delta + 2*sigma^2*kappa_2(w));
    kappa_0(w) = (delta*mu*kappa_1(w)) - 1/2*sigma^2*(kappa_1(w)^2-2*kappa_2(w));

    delta_tilde(w) = sqrt(delta^2 + 2*sigma^2*w);
    mu_tilde(w) = (delta^2*mu) / (delta^2 + 2*sigma^2*w);

    nu_tilde(w) = (X_t*exp(-delta_tilde(w)*s)) + mu_tilde(w)*(1-exp(-delta_tilde(w)*s));
    varsigma_tilde(w) = sigma * sqrt((1-exp(-2*delta_tilde(w)*s)) / (2*delta_tilde(w)));

    A_Lambda(w) = exp(-kappa_1(w)*X_t - kappa_2(w)*X_t^2 - kappa_0(w)*s);
    B_Lambda(w) = 1 - 2*varsigma_tilde(w)^2*kappa_2(w);
    C_Lambda(w) = varsigma_tilde(w)^2/2 ...
                    * (kappa_1(w) + nu_tilde(w)/varsigma_tilde(w)^2)^2;
    D_Lambda(w) = nu_tilde(w)^2 / (2 * varsigma_tilde(w)^2);

    Laplace_Lambda(w) = A_Lambda(w) / sqrt(B_Lambda(w)) ...
                            * exp(C_Lambda(w)/B_Lambda(w) - D_Lambda((w)));

    %% Calculate PMF of N by derivation
    PMF_N = zeros(M_max+1, 1);
    for k_diff = 1 : M_max
        PMF_N(k_diff+1) = (-1)^k_diff / factorial(k_diff) ...
                            * subs(diff(Laplace_Lambda(w), k_diff), 1);
    end
    
end