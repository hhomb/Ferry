%% ================================================================================================
% Code for: Example III.C.b of Toward a Decision Support System for Energy-Efficient 
%           Ferry Operation on Lake Constance based on Optimal Control 
% Authors: Hannes Homburger, Bastian Jäckl, Stefan Wirtensohn, Christian Stopp,
%          Maximilian T. Fischer, Moritz Diehl, Daniel A. Keim, and Johannes Reuter
% Corresponding Author: Hannes Homburger (hhomburg@htwg-konstanz.de)
% Affiliation: HTWG Konstanz, Institute of System Dynamics
% Conference: IEEE European Control Conference ECC, 2026
%
% Description:
% This MATLAB script implements the Example III.C.b described in:
% "Toward a Decision Support System for Energy-Efficient Ferry
%           Operation on Lake Constance based on Optimal Control"
% Accepted for publication at ECC 2026.
%
% Usage:
% This code is intended to reproduce the sensitivity evaluation of Example III.C.b presented in the paper.
%
% License:
% for academic use only
%
% Disclaimer:
% This code is provided for academic and research purposes only, and 
% comes without any warranty or guarantee of performance.
%
% Last updated: 2026-04-09
% ==================================================================================================

%% Set up CasADi
close all
clear  
clc

addpath('C:\Program Files\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

%% Parameters
% load ferry parameters
ferry_param = ferry_parameters;
n_x = 6; % state dimension  
n_u = 3; % input dimension
F_max = 17000; % [N] max. force

% env parameters for scenario 
xi = env_parameters;

%% simulation parameters
N_sim = 20; % number of internal RK4 steps each with same control
N = 90;     % number of simulation steps with different control

% initial state
x_init = zeros(6,1); % see paper for details
x_init(1) = -100;
x_init(2) = -2000;
x_init(3) = pi/2;
x_init(4) = 1;

% goal state
x_end = zeros(6,1);
x_end(2) = 0;
x_end(3) = pi/2;

%% scaling
s_x = [1; 1; 1; 5; 1; 0.1];
s_u = [10000; 10000; 0.1];

Sx = diag(s_x);
Su = diag(s_u);

%% CasADi OPTIMIZATION SETUP
opti = casadi.Opti();

% Declare variables for CasADi 
x_cas = SX.sym('x_cas', n_x, 1);
u_cas = SX.sym('u_cas', n_u, 1);
h_cas = SX.sym('h_cas', 1, 1);

% Declare system dynamics in CasADi
x_next = F(x_cas, u_cas, xi, ferry_param, h_cas);
for k = 1:N_sim-1
    x_next = F(x_next, u_cas, xi, ferry_param, h_cas);
end
F_cas = Function('F_cas', {x_cas, u_cas, h_cas}, {x_next});

% Declare OCP parameters
x_init_param = opti.parameter(n_x, 1);
x_end_param  = opti.parameter(n_x, 1);
T_param      = opti.parameter(1, 1);

% Derive step size for the internal step
h_param = T_param / (N * N_sim);

% Scale parameters
x_init_s_expr = x_init_param ./ s_x;
x_end_s_expr  = x_end_param ./ s_x;

% Opti decision variables (scaled)
X_s = opti.variable(n_x, N+1);
U_s = opti.variable(n_u, N);

% Physical variables (unscaled -> for dynamics/cost)
X = Sx * X_s;
U = Su * U_s;

% Original cost matrices
R_orig = diag([0.001, 0.001, 10000]);
Q_orig = diag([0, 0, 0, 0.1, 0.1, 10000]);
% Scaled cost matrices
R_s = Su * R_orig * Su; 
Q_s = Sx * Q_orig * Sx; 

% setup OCP
opti.subject_to(X_s(:,1) == x_init_s_expr);

Cost = 0;
for k=1:N   
    % Dynamics constraint parameterized by shrinking h
    opti.subject_to(X_s(:,k+1) == F_cas(X(:,k), U(:,k), h_param) ./ s_x);
    
    % Scaled constraints
    opti.subject_to(U_s(1,k)^2 + U_s(2,k)^2 <= (F_max*R_orig(1,1))^2);
    opti.subject_to(X_s(4,k) >= 0);
    
    % Accumulate cost
    Cost = Cost + h_param * N_sim * ( P(U(:,k),ferry_param) + ... % power
           U_s(:,k)'*R_s*U_s(:,k) + ...                           % regularization
           (x_end_s_expr - X_s(:,k))'*Q_s*(x_end_s_expr - X_s(:,k)) );
end

% terminal constraint
opti.subject_to(X_s(:,end) - x_end_s_expr == 0);
% define cost to be minimized
opti.minimize(Cost);

% Configure solver
p_opts = struct('expand', true);
s_opts = struct('max_iter', 1000, 'print_level', 0, 'sb', 'yes'); 

opti.solver('ipopt', p_opts, s_opts);

%% SHRINKING HORIZON MPC LOOP

% Define initial time conditions
T_initial = 900; 

stoer_idx = [6,12,18,24];       % define indizes for disturbances
factor_array = 0.8:0.02:1.2;    % define coefficients for disturbances

Energy_buffer = zeros(length(factor_array),length(stoer_idx)); % alloc buffer to store energy     

for scenario = 1:length(stoer_idx) % over all indizes
    
    for factor=1:length(factor_array) % over all factors
        
        % init variables
        T_current = T_initial;
        xi_sim = xi;
        % vary environment
        xi_sim(stoer_idx(scenario)) = xi_sim(stoer_idx(scenario))*factor_array(factor);
        
        % init plant
        x_plant = x_init;
        
        % allocate historical data for analysis
        mpc_steps = N - 10; % Stop before T=0 to avoid numerical singularity
        X_history = zeros(n_x, mpc_steps + 1);
        U_history = zeros(n_u, mpc_steps);
        T_history = zeros(1, mpc_steps);
        
        X_history(:,1) = x_plant;
        
        % Set the constant goal state
        opti.set_value(x_end_param, x_end);
        
        %% MPC loop
        disp('Starting Shrinking Horizon MPC...');
        solve_times = zeros(1,mpc_steps);
        
        for i = 1:mpc_steps
            
            % Update Parameters for current iteration
            opti.set_value(x_init_param, x_plant);
            opti.set_value(T_param, T_current);
            
            % Warm Start IPOPT
            if i == 1 % first iterate (cold start)
                % Linear interpolation for the very first solve
                x_init_s_num = x_plant ./ s_x;
                x_end_s_num  = x_end ./ s_x;
                opti.set_initial(X_s, (x_end_s_num - x_init_s_num)*linspace(0,1,N+1) + x_init_s_num);
                opti.set_initial(U_s, [0.1*ones(2,N); zeros(1,N)]);
            else
                % Use previous optimal solution as the initial guess
                opti.set_initial(X_s, sol.value(X_s));
                opti.set_initial(U_s, sol.value(U_s));
            end
            
            % Solve OCP
            try
                tic
                sol = opti.solve();
                solve_times(i) = toc;
                steps_completed = i;
            catch
                break;
            end
          
            
            % Get first control action and apply to plant
            u_applied = Su * sol.value(U_s(:,1));
            U_history(:,i) = u_applied;
            T_history(i) = T_current;
            
            % Simulate the plant forward one control interval
            time_elapsed = 10; %[s]
            number_fine_steps = 100;
            for k=1:number_fine_steps
                x_plant = F(x_plant, u_applied, xi_sim, ferry_param,  time_elapsed/number_fine_steps);
            end
            X_history(:,i+1) = x_plant;
            
            % Shrink the horizon time for the next iteration     
            T_current = T_current - time_elapsed;
        end
        
        disp('MPC Simulation Complete');
        
        %% SOLVER PERFORMANCE
        % get warmstarts
        actual_solve_times = solve_times(2:steps_completed);
        
        avg_time = mean(actual_solve_times);
        max_time = max(actual_solve_times);
        
        disp('--- Solver Performance ---');
        disp(['Average solve time: ' num2str(avg_time) ' seconds']);
        disp(['Maximum solve time: ' num2str(max_time) ' seconds']);
        
        
        %% Get solution and compute closed-loop energy
        Energy = 0;
        P_buf = zeros(1, steps_completed);
        
        % Pre-allocate actual time vectors for accurate plotting
        t_vec_u = zeros(1, steps_completed);
        t_vec_x = zeros(1, steps_completed + 1);
        
        current_time = 0;
        t_vec_x(1) = 0;
        
        for k = 1:steps_completed
            % Calculate power for the applied control action
            P_buf(k) = P(U_history(:,k), ferry_param);
            dt = time_elapsed;
            
            % Integrate power to get energy (Ws)
            Energy = Energy + dt * P_buf(k);
            
            % Build the true time vectors
            t_vec_u(k) = current_time;
            current_time = current_time + dt;
            t_vec_x(k+1) = current_time;
        end
        
        Energy_kWh = Energy / 3600000;
        Energy_buffer(factor,scenario) = Energy_kWh;
    end
    
end

%% Visualization
% Plant states
figure('Name', 'MPC Closed-Loop States', 'Color', 'w');
state_labels = {'X_1', 'X_2', 'X_3', 'X_4', 'X_5', 'X_6'};
for k = 1:6
    subplot(6,1,k)
    plot(t_vec_x, X_history(k, 1:steps_completed+1), 'b-', 'LineWidth', 1.5)
    ylabel(state_labels{k}, 'FontWeight', 'bold');
    grid on;
    if k == 6
        xlabel('Time [s]', 'FontWeight', 'bold');
    end
end

% Applied controls
figure('Name', 'MPC Closed-Loop Controls', 'Color', 'w');
control_labels = {'U_1', 'U_2', 'U_3'};
for k = 1:3
    subplot(3,1,k)
    stairs(t_vec_u, U_history(k, 1:steps_completed), 'r-', 'LineWidth', 1.5)
    ylabel(control_labels{k}, 'FontWeight', 'bold');
    grid on;
    if k == 3
        xlabel('Time [s]', 'FontWeight', 'bold');
    end
end

% Computation times
figure('Name', 'MPC Computation Times', 'Color', 'w');
stairs(2:steps_completed, actual_solve_times, 'k-', 'LineWidth', 1.2)
hold on;
yline(avg_time, 'b--', 'Average', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
yline(max_time, 'r--', 'Max', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
xlabel('MPC Step', 'FontWeight', 'bold');
ylabel('Solver Time [s]', 'FontWeight', 'bold');
title('IPOPT Time per MPC Step');
grid on;

%% Subfunctions

% Explicit RK4 step
function x_next = F(x,u,xi,ferry_param,h)
k1 = f(x,u,xi,ferry_param);
k2 = f(x+k1*h/2,u,xi,ferry_param);
k3 = f(x+k2*h/2,u,xi,ferry_param);
k4 = f(x+k3*h,u,xi,ferry_param);
x_next = x + h/6*(k1+2*k2+2*k3+k4);
end

% 2D rotation matrix
function J = J(psi)
J=[cos(psi),-sin(psi),0;
    sin(psi),cos(psi),0;
    0,0,1];
end

% Quadratic model
function deta_c = quadraticModel(xi,eta)

p = eta(1:2);

Q_x = [xi(1),xi(2);xi(2),xi(3)];
L_x = [xi(4),xi(5)];
c_x = xi(6);
Q_y = [xi(7),xi(8);xi(8),xi(9)];
L_y = [xi(10),xi(11)];
c_y = xi(12);


deta_c = [0.5*p'*Q_x*p + L_x*p + c_x;
    0.5*p'*Q_y*p + L_y*p + c_y];

end

% Power model
function P = P(u,ferry_param)
c_p = ferry_param(10);
P = 2*c_p*((1/4)*(u(1)^2+u(2)^2))^(3/4);
end

% Ferry dynamics (continuous time)
function dx = f(x,u,xi,ferry_param)
m = ferry_param(1);

c_x = ferry_param(2);
A_FW = ferry_param(3);
c_y = ferry_param(4);
A_LW = ferry_param(5);

X_uu = ferry_param(6);
X_u = ferry_param(7);
Y_vv = ferry_param(8);
Y_v = ferry_param(9);

Xa = u(1);
Ya = u(2);
dr = u(3);
psi = x(3);
nu = x(4:6);

tau_w = tau_wind(x,xi);
tau_D = tau_hydro(x,xi);

deta = J(psi) * nu;
dnu  = [m^(-1) .* ([Xa;Ya] + tau_w - tau_D);
    dr/1000];
dx = [deta;dnu];

    function tau_w = tau_wind(x,xi)
        epsilon = 0.1;
        eta = x(1:3);
        nu_ = x(4:6);
        psi_=eta(3);
        xi_wind = xi(1:12);
        deta_wind = quadraticModel(xi_wind,eta);
        nu_r =  nu_ - J(-psi_)*[deta_wind;0];
        u_r = nu_r(1);
        v_r = nu_r(2);
        rho = 1.2041;
        tau_w = 1/2* rho.*[-c_x*u_r*sqrt(u_r^2+epsilon)*A_FW;
            -c_y*v_r*sqrt(v_r^2+epsilon)*A_LW];
    end
    function tau_D = tau_hydro(x,xi)
        epsilon = 0.1;
        eta = x(1:3);
        nu_ = x(4:6);
        psi_=eta(3);
        xi_water = xi(13:24);
        deta_water = quadraticModel(xi_water,eta);
        nu_r=  nu_- J(-psi_)*[deta_water;0];
        u_r = nu_r(1);
        v_r = nu_r(2);
        tau_D = [X_uu*u_r*sqrt(u_r^2+epsilon) + X_u*u_r;
            Y_vv*v_r*sqrt(v_r^2+epsilon) + Y_v*v_r];
    end



end

% Environment parameters (specific scenario)
function xi = env_parameters()

windParamsX=...
    -[ -0.0794e-6;
    0.5*(-0.1224e-6);
    -0.0147e-6;
    -1.0387e-3;
    -0.8834e-3;
    4* 1.8894];
windParamsY=...
    [ 0.1966e-6;
    0.5*(0.0780e-6);
    0.1080e-6;
    -0.4946e-3;
    -0.8716e-3;
    4*-1.8080];
waterParamsX=[...
    -0.0034e-6;
    0.5*(-0.0062e-6);
    0.0085e-6;
    0.0365e-3;
    -4*0.0214e-3;
    -0.0218];
waterParamsY=...
    -[-0.0150e-6;
    0.5*(-0.0037e-6);
    0.0048e-6;
    0.0040e-3;
    -0.0361e-3;
    0.0689];

xi=[windParamsX;windParamsY;waterParamsX;waterParamsY];
end

% Ferry parameters (described in paper)
function ferry_param = ferry_parameters()
% see paper for description
m = 35000;

c_x = 0.59 ;
A_FW = 7.3*8.11;
c_y = 0.84;
A_LW = 7.3*30;

X_uu = 409;
X_u = 799;
Y_vv = 2866;
Y_v = 5595;

c_p = 0.104;

ferry_param = [
    m;
    c_x;
    A_FW  ;
    c_y  ;
    A_LW  ;
    X_uu  ;
    X_u  ;
    Y_vv  ;
    Y_v;
    c_p ];

end
