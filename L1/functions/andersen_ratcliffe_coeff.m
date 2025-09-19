function [A, B, alpha_1, gamma_end] = andersen_ratcliffe_coeff(v_j, b_j, r_j, dt, dx, M)
    % vol: vector of local implicit volatility for internal points at time
    % j we do however not have a dependence on the price here.
    % b_j: defined as in article
    % r_j: interest rate for time j
    % dt: time step
    % dx: spacing in x-direction
    
    alpha = dt / (dx^2);
    
    c = ones(M, 1).*(-alpha * v_j);
    u = ones(M, 1).*(0.5 * alpha * (v_j + dx * b_j));
    l = ones(M, 1).*(0.5 * alpha * (v_j - dx * b_j));
    
    
    % Build the tridiagonal matrix M
    M_mat = diag(c) + diag(l(2:end), -1) + diag(u(1:end-1), 1);
    
    % Build matrices A and B for the system: A * H_j = B * H_{j+1} + boundary terms
    A = (1 + r_j * dt) * eye(M) - 0.5 * M_mat;
    B = 0.5 * M_mat + eye(M);
    
    % Boundary coefficients: l1 for left boundary, uN for right boundary
    alpha_1 = 0.5*l(1);
    gamma_end = 0.5*u(M);
end
