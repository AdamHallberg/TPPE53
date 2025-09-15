function [r_out] = riskfree(OIS, T, N)

    % Sort OIS contracts by maturity.
    [~, sorted_indices] = sort([OIS.Maturity]);
    OIS_sorted = OIS(sorted_indices);
    n_ois = length(OIS_sorted);
    
    % Get maturities (tau) and rates (y).
    tau = [OIS_sorted.Maturity];
    y = [OIS_sorted.Mid];
    
    % Initialize vectors for zero-coupon rates and discount factors.
    rates = zeros(1, n_ois);
    P = zeros(1, n_ois);
    
    % First OIS maturity is the base case for bootstrapping.
    rates(1) = y(1);
    P(1) = 1 / (1 + rates(1) * tau(1));

    % Loop through each subsequent OIS maturity to bootstrap the curve.
    for k = 2:n_ois
        
        % Calculate time difference between maturities.
        dt_k = tau(k) - tau(k-1);
        
        % Calculate the present value of fixed leg payments from previous maturities.
        sum_of_fixed_legs_pv = 0;
        for j = 1:k-1
            sum_of_fixed_legs_pv = sum_of_fixed_legs_pv + y(k) * (tau(j) - tau(j-1)) * P(j);
        end
        
        % Solve for the current discount factor using the zero-PV swap condition.
        P(k) = (1 - sum_of_fixed_legs_pv) / (1 + y(k) * dt_k);
        
        % Convert discount factor to a zero-coupon rate.
        rates(k) = -log(P(k)) / tau(k);
    end
    
    % Create an equally spaced time grid for the output.
    time_out = linspace(0, T, N);
    
    % Interpolate the calculated rates onto the output grid.
    r_out = interp1([0, tau], [0, rates], time_out, 'linear', 'extrap');
end
