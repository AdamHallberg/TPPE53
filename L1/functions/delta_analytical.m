function [bsm_delta] = delta_analytical(price, K, T, r, sigma, option_type)
% Black-Scholes analytical function
    d1 = (log(price/K) + (r + 0.5*sigma^2)*T) / (sigma*sqrt(T));
    
    if strcmp(option_type, 'call') || strcmp(option_type, 'Call')
        bsm_delta = normcdf(d1);
    else
        d1 = (log(price/K) + (r + 0.5*sigma^2)*T) / (sigma*sqrt(T));
        bsm_delta = normcdf(d1) - 1;
    end
end