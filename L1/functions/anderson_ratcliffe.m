function [F, x_grid, time, A, B] = anderson_ratcliffe(x_low, x_high, T, N, M, K, r, sigma, option)
    % Discretization
    dt = T / N;
    dx = (x_high - x_low) / (M+1); 
    time = linspace(0, T, N+2)';
    x_grid = linspace(x_low, x_high, M+2)';     

    % F is going to hold all of the option values
    F = zeros(M+2, N+2);
    
    % Terminal Condition, simple for European option from 1c) 
    if strcmp(option, 'Call') || strcmp(option, 'call')
        F(:, end) = max(exp(x_grid) - K, 0);
    else
        F(:, end) = max(K - exp(x_grid), 0);
    end
    
    % Boundary Conditions
    time_to_maturity = time(end) - time;  
    if strcmp(option, 'Call') || strcmp(option, 'call')
        F(end, :) = exp(x_high) - K*exp(-r*time_to_maturity);
    else
        F(1, :) = K*exp(-r*time_to_maturity) - exp(x_low);
    end

    % Determine Coefficients and matrices
    v_j = sigma^2; 
    b_j = r - 0.5*v_j;
    [A, B, B1, Bend] = andersen_ratcliffe_coeff(v_j, b_j, r, dt, dx, M);

    for n = N+2:-1:2
        % Build RHS: B * F_known(2:end-1)
        rhs = B * F(2:end-1, n);
       
        % Handle BC 
        rhs(1) = rhs(1) + B1*(F(1, n) + F(1, n-1));
        rhs(end) = rhs(end) + Bend*(F(end, n) + F(end, n-1));
        
        % Solve for interior points at time n-1
        F_n_int = A \ rhs;

        % Store solution
        F(2:end-1, n-1) = F_n_int;
    end

end
