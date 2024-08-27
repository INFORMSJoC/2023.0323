%% Initialise environment
clc
clear
format long

%% Set parameters
X_0 = 0.5;
mu = 0.8;
delta = 1.1;
sigma = 1.3;

N_jumps_max = 100000; % set enough number of jumps for each path within time T_stop
T_input = 0.5;

N_m = 4; % index of number of grids

%% Define parameters
delta_tilde = sqrt(delta^2 + 2*sigma^2);
mu_tilde = mu * delta^2 / delta_tilde^2;
sigma_tilde = sigma;

kappa_2 = 1 / (delta_tilde+delta);
kappa_1 = 2*mu*delta*kappa_2 / delta_tilde;
kappa_0 = mu*mu_tilde + sigma_tilde^2*kappa_2;

%% Simulate point process via discretisation scheme
Expectation_N_con_T=zeros(N_m,1);
Mean_simulated_N_T=zeros(N_m,1);
SD_simulated_N_T=zeros(N_m,1);
Time_Loop=zeros(N_m,1);
Bias=zeros(N_m,1);
RMSE=zeros(N_m, 1);
Output=zeros(N_m,6);
for m=1:N_m
    N_grid=250*m; % total number of grids in time space
    h=T_input/N_grid; % length of each grid
    
    N_path= 15625*4^(m-1);   % number of simulated sample paths
    N_T = zeros(N_path,1);  % record point N at time T
    
    time_CPU_start = cputime;
    for j=1:N_path
        % Define Variables
        X_OU_process=zeros(N_grid+1,1);
        X_OU_process(1)=X_0;
        lambda_process= zeros(N_grid+1,1);  % record discretized intensity lambda
        lambda_process(1)=X_0^2;
        T_jump = zeros(N_jumps_max+1,1);  % default/jump arrival times in point process
        lambda_process_n = zeros(N_jumps_max+1,1);  % left-limit of intensity at default/jump arrival times
        
        % Simulate one path
        index_jump=1;  % record the index of arrival times
        for i_jump=1:N_jumps_max
            compensator_hat_change = 0; % compensator change between two successive jumps
            E = -log(rand); % Exp(1)
            
            for i_grid=(index_jump+1):N_grid+1 % time discretization
                 
                X_OU_process(i_grid)=X_OU_process(i_grid-1)-delta*(X_OU_process(i_grid-1)-mu)*h+sigma*sqrt(h)*randn; 
                lambda_process(i_grid) = ( X_OU_process(i_grid))^2;
                compensator_hat_change = compensator_hat_change + h*lambda_process(i_grid); % compensator
                
                if (compensator_hat_change>E)
                    T_jump(i_jump)=i_grid*h; % record the arrival times
                    index_jump = i_grid; % record the index of arrival times
                    lambda_process_n(i_jump) = lambda_process(i_grid);
                    lambda_process(i_grid) = lambda_process_n(i_jump);
                    break
                end
            end
            
            % Terminate loop before time T_stop
            if (i_grid == N_grid+1) % terminate loop execution for one path before time T_stop
                break
            end
            
        end
        
        N_T(j)=i_jump-1; % record N_T for each path
        
    end
    time_looping=cputime-time_CPU_start;
    
    % Validation
    Time_Loop(m)=time_looping;
    Expectation_N_con_T(m)=1/(2*delta)*(2*mu^2*delta*T_input+4*mu*(X_0-mu)*(1-exp(-delta*T_input))+ (X_0-mu)^2*(1-exp(-2*delta*T_input))+sigma^2/(2*delta)*(2*delta*T_input+exp(-2*delta*T_input)-1));
    Mean_simulated_N_T(m)=mean(N_T); % mean of simulated point
    Bias(m)=mean(N_T-Expectation_N_con_T(m));   
    SD_simulated_N_T(m)=std(N_T)/sqrt(N_path); % SD of simulated point
    RMSE(m)=sqrt(Bias(m)^2+SD_simulated_N_T(m)^2);
    Output(m,:)=[Expectation_N_con_T(m)  Mean_simulated_N_T(m) Bias(m) SD_simulated_N_T(m) RMSE(m) Time_Loop(m)]; % output
    
end

Output