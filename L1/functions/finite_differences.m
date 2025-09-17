function [F, price, time] = finite_differences(S_low, S_high, T, N, M, K, r, sigma, option, q)
% FINITE_DIFFERENCES - Finite difference method for option pricing
    % Discretization
    dt = T / N;
    dS = (S_high - S_low) / M;
    time = linspace(0, T, N+1)';
    price = linspace(S_low, S_high, M+1)';

    % F is going to hold all of the option values
    F = zeros(M+1, N+1);
    
    % Terminal condition and Boundary Condition
    time_to_maturity = time(end) - time;
    if strcmp(option, 'Call') || strcmp(option, 'call')
        F(:, end) = max(price - K, 0);
        F(end, :) = S_high - K*exp(-r*time_to_maturity);
    else
        F(:, end) = max(K - price, 0);
        F(1, :) = K*exp(-r*time_to_maturity) - S_low;
    end

    % Get the matrices describing the dynamics
    [A, B, alpha1, gamma_end] = crank_nicholson_coeff(price(2:end-1), r, sigma, dt, dS, q);

    for n = N+1:-1:2
        % Build RHS: B * F_known(2:end-1)
        rhs = B * F(2:end-1, n);

        % A F_n = BF_{n+1} = rhs => F_n = A\rhs
        rhs(1) = rhs(1) + 0.5*dt*alpha1*(F(1, n) + F(1, n-1));
        rhs(end) = rhs(end) + 0.5*dt*gamma_end*(F(end, n) + F(end, n-1));
        
        % Solve for interior points at time n-1
        F_n_int = A \ rhs;

        % Store solution
        F(2:end-1, n-1) = F_n_int;
    end
end
