function [A, B, alpha_1, gamma_end] = crank_nicholson_coeff(S_j, r, sigma, dt, dS, q)
    % Create Matrices and Coefficients for the Crank-Nicolson scheme (theta = 0.5)
    %
    % Correct form: A*F^{n} = B*F^{n+1}

    M = length(S_j);
    
    % Coefficients for L without dt
    a = (sigma^2 * S_j.^2) / (2 * dS^2) - ((r-q) * S_j) / (2 * dS);
    b = - (sigma^2 * S_j.^2) / (dS^2) - r;
    c = (sigma^2 * S_j.^2) / (2 * dS^2) + ((r-q) * S_j) / (2 * dS);
    
    % Build L_mat (tridiagonal matrix for L)
    L_mat = diag(b) + diag(a(2:end), -1) + diag(c(1:end-1), 1);
    
    % Crank-Nicolson matrices
    A = eye(M) - 0.5 * dt * L_mat;
    B = eye(M) + 0.5 * dt * L_mat;

    % Coeficcients needed for handeling the boundary conditions. By ooking
    % at the matrices in a) it is obvious where these come from. 
    alpha_1 = (sigma^2 * S_j(1)^2) / (2 * dS^2) - ((r-q) * S_j(1)) / (2 * dS);
    gamma_end = (sigma^2 * S_j(end)^2) / (2 * dS^2) + ((r-q) * S_j(end)) / (2 * dS);
    
end
