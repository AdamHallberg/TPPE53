function vol_surface = implied_volatility_surface(option_data, S0)
%IMPLIED_VOLATILITY_SURFACE Fit a smooth, differentiable surface to implied vols
%
% option_data: [strike, imp_vol(%), T(years)]
% S0         : spot price
%
% Returns:
%   vol_surface.fit   : tpaps spline of total variance w(k,T)
%   vol_surface.eval  : function @(K,T) implied vol σ(K,T)
%   vol_surface.S0    : spot
%   vol_surface.plot  : function handle to plot surfaceer

    % Extract columns
    strike = option_data(:,1);
    T      = option_data(:,3);
    sigma  = option_data(:,2) * 0.01;   % convert % to decimal

    % Work in log-moneyness (better behaved than raw moneyness)
    k = log(strike ./ S0);

    % Total implied variance
    w = (sigma.^2) .* T;

    % Build thin-plate smoothing spline: tpaps([x;y],z,p)
    % p = 1 means exact interpolation, <1 means smoothing
    p  = 0.999;
    xy = [k'; T'];       % 2 x N
    fn = tpaps(xy, w', p);

    % Build evaluation function (returns implied vol)
    eval_sigma = @(K,Tq) ...
        sqrt(max(fnval(fn, [log(K./S0)'; Tq']),0) ./ max(Tq',eps));

    % Package result
    vol_surface = struct();
    vol_surface.fit  = fn;          % spline for total variance
    vol_surface.eval = eval_sigma;  % σ(K,T)
    vol_surface.S0   = S0;

    % Plot helper
    vol_surface.plot = @() plot_surface(fn, S0, T);

end

function plot_surface(fn, S0, T_data)
    % Plot a nice surface of σ(K,T)
    k_vec = linspace(-0.5,0.5,50);   
    T_vec = linspace(0.01, max(T_data), 50);   % use data max instead of fn.knots
    
    [Kgrid,Tgrid] = meshgrid(exp(k_vec)*S0, T_vec);
    
    % Loop evaluate grid
    Wgrid = zeros(size(Kgrid));
    for i=1:numel(Kgrid)
        Wgrid(i) = fnval(fn, [log(Kgrid(i)/S0); Tgrid(i)]);
    end
    
    SIGgrid = sqrt(max(Wgrid,0)./max(Tgrid,eps));
    figure;
    surf(Kgrid, Tgrid, SIGgrid);
    shading interp; colormap turbo; colorbar;
    xlabel('Strike'); ylabel('Maturity T'); zlabel('Implied vol');
    title('Smoothed implied volatility surface');
    view(45,30);
end