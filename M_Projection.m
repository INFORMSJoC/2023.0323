%% Initialise environment
clc
clear
format long

%% Set parameters
X_0_OU = 0.5; % state [X_t] of OU process
mu_OU = 0;  
delta_OU = 1.1;
sigma_OU = 1.3;
t_target = 0.01; % target time [t] for simulating point process

N_jumps_max = 2;

%% Define parameters
lambda_0 = X_0_OU^2;  
delta = 2*delta_OU;
a = sigma_OU^2 / delta;
sigma = 2*sigma_OU; % parameters ([delta], [a], [sigma])  

kappa = sqrt(delta^2 + 2*sigma^2);  

%% Define functions  
syms z s t
a_s_z = 2*a*delta/(sigma^2)*log(  2*kappa*exp((kappa+delta)/2*s)/(z*sigma^2*(exp(kappa*s)-1)+kappa-delta+(kappa+delta)*exp(kappa*s) )  );   % function a(s,z)
b_s_z = ( z*(kappa+delta+(kappa-delta)*exp(kappa*s))+2*(exp(kappa*s)-1) )/( z*sigma^2*(exp(kappa*s)-1)+kappa-delta+(kappa+delta)*exp(kappa*s) );  % function b(s,z)

%% Simulate point process via projection scheme
N_path = 1000; % number of simulated sample paths
N_T = zeros(N_path,1);  % record point N at time T

time_CPU_start = cputime;
for j_path=1:N_path % loop for N_path sample paths of simulation
    T_jump = zeros(N_jumps_max+1,1);  % default/jump times
    M_z=exp(-z*lambda_0); % initial for function M_T(z) % ezplot(M_z,[0,5])
    
    for i_jump=1:(N_jumps_max+1)
        % Calculate projection h_t
        a_t_z = subs(a_s_z,s,t-T_jump(i_jump));
        b_t_z = subs(b_s_z,s,t-T_jump(i_jump));
        D_t_z = exp(a_t_z)*subs(M_z,z,b_t_z);
        D_t_0 = subs(D_t_z,z,0);
        M_z_t = D_t_z/D_t_0;      % function M_t(z)
        dM_z = diff(M_z_t,z);
        h_t = - subs(dM_z,z,0); % projection function % ezplot(h_t,[0,5])
        
        % Thinning
        t_num = T_jump(i_jump); % initialise
        accept = false;
        while accept==false
            B_bound = double(subs(h_t,t,t_num)); % vpa(subs(h_t,t,t_num));
            Exp_B_bound = -1/B_bound*log(rand);
            t_temp = t_num + Exp_B_bound;
            U_rand = rand;
            if (U_rand*B_bound <= double(subs(h_t,t,t_temp)))  % vpa(subs(h_t,t,t_temp))
                accept = true;
                % t_temp
                T_jump(i_jump+1) = t_temp; % jump occures at t_temp
            else
                % accept = false
                t_num = t_temp;
            end
        end
        if T_jump(i_jump+1)>t_target
             N_T(j_path) = i_jump-1; % point process value N at time T_stop
        break
        end

        % Update M_t(z) function
         M_z_n = subs(M_z_t,t,T_jump(i_jump+1));
         dM_z_n = diff(M_z_n,z);
         M_z = dM_z_n/(subs(dM_z_n,z,0));
        
    end
end
time_looping=cputime-time_CPU_start
 
%% Validation
Expectation_N_con_T=(a*delta)/delta *t_target + 1/delta*( lambda_0 - (a*delta)/delta ) * ( 1-exp(-delta*t_target) ) % xi neq 0

Mean_simulated_N_T=mean(N_T) % mean of simulated point
