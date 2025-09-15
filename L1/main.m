%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Adam Hallberg, Oscar Ljungdahl
% Date:    2025-09-13
% Status:  Incomplete
%
% Comments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup
close all;
clear;
rng('default');
format long
clc;

%% Read OIS data
[OIS, option_data] = get_ois("ois_data.xlsx");

%% Setup variables
K = 100;
S0 = 100;
T = 1;
N = 100;
M = 100; % kanskeändra detta till beroende av options datan
r = riskfree(OIS, T);
r(1:floor(T*365)); % only use the rates we need 
sigma = 0.2;
q = 0;
option = 'call';
type = "eu";

% Caculate S_low and S_high to satisfy P(S(T) not in [S_low, S_high]) =
% 0.999
[S_low, S_high] = price_bounds(S0, r(end), sigma, T, 1-0.999);
S_low = floor(S_low); S_high = ceil(S_high);

%% Calculate prices
% Our FD prices
[F, price, time] = finite_differences(S_low, S_high, T, N, M, K,...
                                      r, sigma, option, q, type);

% Analytical solution
analytical_results = bsm_analytical(price, K, T, r, sigma, option);

%% Comparision of FD and Analytical
figure;
plot(price, F(:,1), 'b-', 'LineWidth', 2);
hold on;
plot(price, analytical_results, 'r--', 'LineWidth', 2);
legend('FD-method', 'analytical solution');
xlabel('Stock Price');
ylabel('Option Price');
title(sprintf(['Call Option Price Comparison (K=%d, T=%.1f, ...' ...
                    'σ=%.1f, r=%.2f)'], K, T, sigma, r));
grid on;


%% Surf plot
figure;
surf(time, price, F);
hold on;
xlabel('Time');
ylabel('Spot Price');
zlabel('Value of Option')
title(sprintf('%s Option Price Over Time and Spot ', option));
grid on;
