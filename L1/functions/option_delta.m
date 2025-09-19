function delta = option_delta(F, price, log)    
    if log
        % \frac{\partial F(t, s)}{\partial S} =/ x = ln(S) => dx =
        % \frac{dS}{S}  and S = exp(x)/=  exp(-x) \frac{\partial F(t,
        % x)}{dx}
        [Mp1, Np1] = size(F);
        delta = zeros(Mp1, Np1);
        dx = price(2)-price(1); %log space still uniform

        % inner points
        for j = 2:Mp1-1
            x = exp(-mean([price(j+1), price(j)]));
            delta(j, :) = (F(j+1, :) - F(j-1, :)) / (2 * dx);
            delta(j, :) = delta(j, :) * x;
        end
        
        % Forward diff at S_min
        x = exp(-mean([price(2), price(1)]));

        delta(1, :) = (F(2, :) - F(1, :)) / dx;
        delta(1, :) = delta(1, :) * x;
        
        % Backward diff S_max  
        x = exp(-mean([price(end), price(end-1)]));

        delta(end, :) = (F(end, :) - F(end-1, :)) / dx;
        delta(end, :) = delta(end, :) * x;
    else
        [Mp1, Np1] = size(F);
        delta = zeros(Mp1, Np1);
        dS = price(2) - price(1); % equidistant
        
        % inner points
        for j = 2:Mp1-1
            delta(j, :) = (F(j+1, :) - F(j-1, :)) / (2 * dS);
        end
        
        % Forward diff at S_min
        delta(1, :) = (F(2, :) - F(1, :)) / dS;
        
        % Backward diff S_max  
        delta(end, :) = (F(end, :) - F(end-1, :)) / dS;
    end
end