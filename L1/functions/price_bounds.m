function [S_low, S_high] = price_bounds(S0, r, sigma, T, alpha)
    % from task 1b)
    % alpha is the probability level

    z_low = norminv(alpha/2);      % \Phi^-1(α/2)
    z_high = norminv(1 - alpha/2); % \Phi^-1(1 - α/2)
    
    ex  = (r - 0.5*sigma^2)*T;
    vol = sigma * sqrt(T);
    
    S_low = exp(vol * z_low + log(S0) + ex);
    S_high = exp(vol * z_high + log(S0) + ex);
end