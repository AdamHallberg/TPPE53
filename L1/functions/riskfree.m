function [r_out] = riskfree(OIS, T)

    % Sort OIS contracts by maturity
    OIS = sortrows(OIS, 1);
    n_ois = height(OIS);
    OIS = table2array(OIS);
    
    % Get maturities (tau in years) and OIS rates (y)
    tau = OIS(:,1);
    y   = OIS(:,2);
    
    % Initialize vectors for zero-coupon rates and discount factors
    rates = zeros(n_ois,1);
    P     = zeros(n_ois,1);
    
    % First maturity: direct bootstrap
    rates(1) = y(1);
    P(1)     = 1 / (1 + y(1)*tau(1));
    
    % Bootstrap subsequent maturities
    for k = 2:n_ois
        sum_of_fixed_legs_pv = 0;
        for j = 1:k-1
            dt_j = tau(j+1) - tau(j); % delta mellan betalningar
            sum_of_fixed_legs_pv = sum_of_fixed_legs_pv + y(k)*dt_j*P(j);
        end
        
        dt_k = tau(k) - tau(k-1);
        P(k) = (1 - sum_of_fixed_legs_pv) / (1 + y(k)*dt_k);
        
        rates(k) = -log(P(k)) / tau(k);
    end
    
    % Output grid (daily steps until last maturity)
    T = max(tau);
    time_out = (1:floor(T*360))'/360;
    
    % Lägg till punkt vid 0 manuellt
    tau_all = [0; tau];
    P_all   = [1; P];
    
    % Ta bort dubletter (behåll första förekomst)
    [tau_unique, idx] = unique(tau_all, 'stable');
    P_unique = P_all(idx);
    
    % Flat forward interpolation (log-linear på P)
    P_out = exp(interp1(tau_unique, log(P_unique), time_out, 'spline', 'extrap'));
    r_out = -log(P_out) ./ time_out;
    
    
    %plot(1:max(size(r_out)), r_out) % Ser nästan rimlig ut iallafall
    


end
