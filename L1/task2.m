%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Adam Hallberg, Oscar Ljungdahl
% Date:    2025-09-19
% Status:  Complete
%
% Comments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup
close all;
clear;
rng('default');
format long
addpath('functions/');
clc;

%% Read OIS data
[OIS, option_data] = read_excel_file("ois_data.xlsx");

%% Setup variables
K = 2680;
S0 = 2600;
T = 1;
N = floor(T*365); % Discretization in time
M = 100; % 
r = riskfree(OIS, T);
[sigma, k_hat] = implied_volatility(option_data, K);
r = r(1:floor(T*365)); % only use the rates we need 
sigma = sigma(1:floor(T*365)); % same here
option = 'Call';

% Caculate S_low and S_high to satisfy P(S(T) not in [S_low, S_high]) =
% 0.999
[S_low, S_high] = price_bounds(S0, r(1), sigma(1), T, 1-0.999);
S_low = floor(S_low); S_high = ceil(S_high);

% Price boundaries, 
S0 = log(2600);
x_low = log(S_low);
x_high = log(S_high);

%% Calculate prices
% Our FD prices
clc
[F, x_grid, time, A, B, o] = anderson_ratcliffe(x_low, x_high, T, N, M, K,...
                                      r(1), sigma(1), option);


exp_price = exp(x_grid);

% Analytical solution
analytical_results = bsm_analytical(exp_price, K, T, r(1), sigma(1), option);

%% Comparision of FD and Analytical
figure;
plot(exp_price, F(:,1), 'b-', 'LineWidth', 2);
hold on;
plot(exp_price, analytical_results, 'r--', 'LineWidth', 2);
legend('Anderson-Ratcliffe FD-method', 'analytical solution', Location='best');
xlabel('Spot Price');
ylabel('Option Price');
title(sprintf('%s Option Price Comparison ', option));
grid on;


%% Surf plot
figure;
surf(time, exp_price, F);
hold on;
shading interp              
colormap jet              
colorbar                    
xlabel('Time');
ylabel('Spot Price');
zlabel('Value of Option')
title(sprintf('%s Option Price Over Time and Spot ', option));
view(45,30)                 
grid on;



%% Delta Calculation

% Approximation
delta_approx = option_delta(F, x_grid, true);
delta_bsm = delta_analytical(exp_price, K, T, r(end), sigma(end), option);

%% Delta Comparison
figure;
plot(exp_price, delta_approx(:, 1), 'b-', 'LineWidth', 2);
hold on;
plot(exp_price, delta_bsm, 'r--', 'LineWidth', 2);
legend('Anderson-Ratcliffe FD-method', 'analytical solution');
xlabel('Stock Price');
ylabel('Option Price');
title(sprintf(['Call Option Price Comparison (K=%d, T=%.1f, ...' ...
                    'σ=%.1f, r=%.2f)'], K, T, sigma, r));
grid on;

%% Surf plot of delta
figure;
surf(time, exp_price, delta_approx);
hold on;
shading interp              
colormap jet              
colorbar     
xlabel('Time');
ylabel('Spot Price');
zlabel('Delta of Option')
title(sprintf('%s Option Delta Over Time and Spot ', option));
grid on;



