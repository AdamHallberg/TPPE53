function delta = option_delta(F, price)
    % OPTION_DELTA - Option delta w central differens
    
    [Mp1, Np1] = size(F);
    delta = zeros(Mp1, Np1);
    dS = price(2) - price(1);
    
    % inner points
    for j = 2:Mp1-1
        delta(j, :) = (F(j+1, :) - F(j-1, :)) / (2 * dS);
    end
    
    % Forward diff at S_min
    delta(1, :) = (F(2, :) - F(1, :)) / dS;
    
    % Backward diff S_max  
    delta(end, :) = (F(end, :) - F(end-1, :)) / dS;
end