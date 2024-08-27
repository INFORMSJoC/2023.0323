%% Initialise environment
clc
clear
format long

%% Set parameters
X_0 = [0.5, -0.5, 0.5, -0.5];
mu = repmat(0.8, 1, 4);
delta = repmat(1.1, 1, 4);
sigma = repmat(1.3, 1, 4);
t = [1, 1, 0.5, 0.5];

M_max = 6;  % largest m in P(N=m|X_t) to be calculate

%% Calculate PMFs of point process [N_t] with quadratic OU intensity
PMF_N = zeros(M_max+1, 4);
time_PMF_N = zeros(1, 4);
for k = 1 : 4
    tic
    PMF_N(:, k) = f_PMF_N(X_0(k), delta(k), mu(k), sigma(k), t(k), M_max);
    time_PMF_N(k) = toc;
end

%% Plot PMFs of point process [N_t] with quadratic OU intensity
for k = 1 : 4
    subplot(2, 2, k)
    plot(linspace(0, M_max), spline(0:M_max, PMF_N(:, k), linspace(0, M_max)))
    grid minor
end
