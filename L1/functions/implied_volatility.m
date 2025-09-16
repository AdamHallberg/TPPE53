function [vol_series, K_closest] = implied_volatility(option_data, K)

    strike = option_data(:,1);
    imp_vol = option_data(:,2);
    time = option_data(:,3);
    
    K_vec = linspace(min(strike), max(strike), max(size(strike)));
    T_vec = linspace(min(time), max(time), max(size(time)));
    
    [X, Y] = meshgrid(K_vec, T_vec);
    
    vol = griddata(strike, time, imp_vol, X, Y, 'linear');
    
    % Fill NaN holes with nearest interpolation
    nan_idx = isnan(vol);
    vol(nan_idx) = griddata(strike, time, imp_vol, X(nan_idx), Y(nan_idx), 'nearest');
    
    % for debugging
    if false
        % Plot surface
        figure;
        surf(X, Y, vol)
        
        % Make it pretty
        shading interp              % smooth shading, removes grid lines
        colormap turbo              % colorful colormap (try "parula" or "jet" too)
        colorbar                    % show vol scale
        xlabel('Strike')
        ylabel('Maturity (T)')
        zlabel('Implied Volatility')
        title('Implied Volatility Surface')
        view(45,30)                 % set viewing angl
    end

    % Get the implied volatility time series for our strike.
    [~, idx] = min(abs(K_vec - K));   
    K_closest = K_vec(idx);           
    vol_series = vol(:, idx)*0.01; % also convert to decimal         
end

