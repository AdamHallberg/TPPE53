function [A, B, alpha_1, gamma_end] = andersen_ratcliffe_coeff(S_j, r, sigma, i, dt, dS)
    % Create Matrices and Coefficients for the Crank-Nicolson scheme (theta = 0.5)
    %
    % Correct form: A*F^{n} = B*F^{n+1}
    v_hat = sigma(i)^2;       
    b_hat = r(i) - 0.5*v_hat;  % N long

    M = length(S_j);

    % Coefficients for L without dt
    a = (v_hat) / (2 * dS^2) - b_hat / (2 * dS);
    b = - (v_hat) / (dS^2) - r;
    c = (v_hat) / (2 * dS^2) + b_hat / (2 * dS);
    
    % Build L_mat (tridiagonal matrix for L)
    L_mat = diag(b) + diag(a(2:end), -1) + diag(c(1:end-1), 1);
    
    % Crank-Nicolson matrices
    A = eye(M) - 0.5 * dt * L_mat;
    B = eye(M) + 0.5 * dt * L_mat;

    % Coeficcients needed for handeling the boundary conditions. By ooking
    % at the matrices in a) it is obvious where these come from. 
    alpha_1 = a;
    gamma_end = c;
end
