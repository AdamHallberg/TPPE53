%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Adam Hallberg, Oscar Ljungdahl
% Date:    2025-09-19
% Status:  Incomplete
%
% Comments:
%   Delta calculations are "wrong". To get normal delta in [-1, 1] we need
%   to handle the discontinuity.
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
payout = 1000;
K = 2680;
S0 = 2600;
T = 1;
N = floor(T*365); % Discretization in time
M = 200; % 
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
[F_con, ~, ~, ~, ~] = cash_or_nothing(x_low, x_high, T, N, M, K,...
                                      r(1), sigma(1), option, payout);

[F_aon, x_grid, time, A, B] = asset_or_nothing(x_low, x_high, T, N, M, K,...
                                      r(1), sigma(1), option);


% Analytical solution
analytical_aon = 0;
analytical_con = 0;
%% Comparision of FD and Analytical
figure;
plot(exp(x_grid), F_aon(:,1), 'b-', 'LineWidth', 2);
hold on;
plot(exp(x_grid), analytical_aon, 'r--', 'LineWidth', 2);
legend('FD-method', 'analytical solution', Location='best');
xlabel('Spot Price');
ylabel('Option Price');
title(sprintf('%s Asset or Nothing Price Comparison ', option));
grid on;

figure;
plot(exp(x_grid), F_con(:,1), 'b-', 'LineWidth', 2);
hold on;
plot(exp(x_grid), analytical_con, 'r--', 'LineWidth', 2);
legend('FD-method', 'analytical solution', Location='best');
xlabel('Spot Price');
ylabel('Option Price');
title(sprintf('%s Cash or Nothing Price Comparison ', option));
grid on;



%% Surf plot
figure;
surf(time, exp(x_grid), F_aon);
hold on;
shading interp              
colormap jet              
colorbar                    
xlabel('Time');
ylabel('Spot Price');
zlabel('Value of Option')
title(sprintf('%s Asset or Nothing Price Over Time and Spot ', option));
view(45,30)                 
grid on;

figure;
surf(time, exp(x_grid), F_con);
hold on;
shading interp              
colormap jet              
colorbar                    
xlabel('Time');
ylabel('Spot Price');
zlabel('Value of Option')
title(sprintf('%s Cash or Nothing Price Over Time and Spot ', option));
view(45,30)                 
grid on;




%% Delta Calculation

% Approximation
delta_approx_aon = option_delta(F_aon, x_grid, true);
delta_approx_con = option_delta(F_con, x_grid, true);

%% Delta Comparison
figure;
plot(exp(x_grid), delta_approx_aon(:, 1), 'b-', 'LineWidth', 2);
hold on;
plot(exp(x_grid), delta_approx_con(:, 1), 'b-', 'LineWidth', 2);
legend('Asset or Nothing', 'Cash or Nothing');
xlabel('Stock Price');
ylabel('Delta');
title(sprintf('Comparision of AoN and CoN'));
grid on;

%% Surf plot of delta
figure;
surf(time, exp(x_grid), delta_approx_aon);
hold on;
shading interp              
colormap jet              
colorbar     
xlabel('Time');
ylabel('Spot Price');
zlabel('Delta')
title(sprintf('%s Asset or Nothing Delta Over Time and Spot ', option));
grid on;

figure;
surf(time, exp(x_grid), delta_approx_aon);
hold on;
shading interp              
colormap jet              
colorbar     
xlabel('Time');
ylabel('Spot Price');
zlabel('Delta')
title(sprintf('%s Cash or Nothing Delta Over Time and Spot ', option));
grid on;



