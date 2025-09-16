function [F, price, time] = anderson_ratcliffe(S_low, S_high, T, N, M, K, r, sigma, option, q, type)
% FINITE_DIFFERENCES - Finite difference method for option pricing
%
% In derivation we handle the system AF_n = BF_{n+1} but here we are
% looking at the n := n-1 basically. 

    % We start by defining our discretization
    dt = T / N;
    dS = (S_high - S_low) / M;
    time = linspace(0, T, N+1)';
    price = linspace(S_low, S_high, M+1)';

    % F is going to hold all of the option values
    F = zeros(M+1, N+1);
    
    % Terminal condition, simple for European option from 1c) (same for US)
    if strcmp(option, 'Call') || strcmp(option, 'call')
        F(:, end) = max(price - K, 0);
    else
        F(:, end) = max(K - price, 0);
    end
    
    % Internal grid points
    S_j = price(2:end-1);  % Points j=2 to M
  
    
    % Here, A and B are (M-1)x(M-1) and handle the inner points and alpha_1
    % and gamma_end are what is needed to handle the BC ad price(1) and
    % price(end)

    % Now we are to solve this system from time N+1 where we have knowledge
    % about the terminal conditions and move backward to "our time" or what
    % we could call t_0. We only step to n=2 because in step n=2 we obtain
    % the values for n=1 which is what we are intrested in finding out.

    for n = N+1:-1:2

        % Get the matrices describing the dynamics
        [A, B, alpha1, gamma_end] = andersen_ratcliffe_coeff(S_j, r, sigma, n, dt, dS);

        time_to_maturity = (n-2) * dt;  % Time to maturity at step n-1
        F_known = F(:, n); % F_known contains solution at time step n (known)
        
        % Boundary conditions at time step n (known), Ev diskonteringsfel
        if strcmp(option, 'Call') || strcmp(option, 'call')
            F_low_known = 0;
            F_high_known = S_high - K*exp(-r*(time_to_maturity + dt));
        else
            F_low_known = K*exp(-r*(time_to_maturity + dt)) - S_low;
            F_high_known = 0;
        end
        
        % Boundary conditions at time step n-1 (unknown) we talk about
        % these conditions in task 1c)
        if strcmp(option, 'Call') || strcmp(option, 'call') 
            F_low = 0;
            F_high = S_high - K*exp(-r*time_to_maturity);
        else
            F_low = K*exp(-r*time_to_maturity) - S_low;
            F_high = 0;
        end
        
        % Build RHS: B * F_known(2:end-1)
        rhs = B * F_known(2:end-1);
       
        % Handle BC 
        rhs(1) = rhs(1) + 0.5*dt*alpha1*(F_low_known + F_low);
        rhs(end) = rhs(end) + 0.5*dt*gamma_end*(F_high_known + F_high);
        
        % Solve for interior points at time n-1
        F_n_int = A \ rhs;

        if(strcmp(type, "us"))
            if(strcmp(option, "Call") || strcmp(option, "call)"))
                payoff_int = max(S_j - K, 0); 
            else
                payoff_int = max(K - S_j, 0);
            end

            % Fel med diskonterin här

            F_n_int = max(payoff_int, F_n_int);
        end

        % Store solution
        F(2:end-1, n-1) = F_n_int;
        F(1, n-1) = F_low;
        F(end, n-1) = F_high;
    end

end
