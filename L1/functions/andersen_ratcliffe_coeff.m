function [A, B, alpha_1, gamma_end] = andersen_ratcliffe_coeff(S_j, b_hat, v_hat, dt, dS)
    % Create Matrices and Coefficients for the Crank-Nicolson scheme (theta = 0.5)
    %
    % Correct form: A*F^{n} = B*F^{n+1}

    M = length(S_j);

    alpha = dt/(dS^2);

    % Coefficients for L without dt
    u = 0.5 * alpha * (v_hat + dS*b_hat);
    c = - alpha * (v_hat);
    l = 0.5 * alpha * (v_hat - dS*b_hat);
    
    % Build L_mat (tridiagonal matrix for L)
    L_mat = diag(c) + diag(u(2:end), -1) + diag(l(1:end-1), 1);
    
    % Crank-Nicolson matrices
    A = eye(M) - 0.5 * dt * L_mat;
    B = eye(M) + 0.5 * dt * L_mat;

    % Coeficcients needed for handeling the boundary conditions. By looking
    % at the matrices in a) it is obvious where these come from. 
    alpha_1 = (sigma^2 * S_j(1)^2) / (2 * dS^2) - ((r-q) * S_j(1)) / (2 * dS);
    gamma_end = (sigma^2 * S_j(end)^2) / (2 * dS^2) + ((r-q) * S_j(end)) / (2 * dS);
    
end
