import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt


def thomas_algorithm(A, d):
    """
    Solve tridiagonal system A·x = d using Thomas algorithm.

    This algorithm is written by ChatGPT based on a old 
    thomas lagorithm implementation I had in matlab 
    
    Parameters:
    A: tridiagonal matrix (M-1 x M-1)
    d: right-hand side vector (length M-1)
    
    Returns:
    x: solution vector
    """
    n = len(d)
    
    # Extract diagonals
    a = np.diag(A, k=-1)  # sub-diagonal (length n-1)
    b = np.diag(A, k=0)   # main diagonal (length n)
    c = np.diag(A, k=1)   # super-diagonal (length n-1)
    
    # Make copies
    a = a.copy()
    b = b.copy()
    c = c.copy()
    d = d.copy()
    
    # Forward elimination
    for i in range(1, n):
        w = a[i-1] / b[i-1]
        b[i] -= w * c[i-1]
        d[i] -= w * d[i-1]
    
    # Back substitution
    x = np.zeros(n)
    x[-1] = d[-1] / b[-1]
    
    for i in range(n-2, -1, -1):
        x[i] = (d[i] - c[i] * x[i+1]) / b[i]
    
    return x


def finite_differences(S_low, S_high, T, N, M, K, r, sigma, option="call", theta = 0.5):
    """
    INPUT: 
        ...

    OUTPUT: 
        Vector of price of option for different S_0 for inputed K.

    NOTES:
        Defaults to Crank-Nicholson for call options. 
    """

    #   Discritization 
    dt = T / N
    dS = (S_high - S_low) / M
    time = np.linspace(0, T, N+1)
    price = np.linspace(S_low, S_high, M+1)

    # Setup terminal condition
    if option == "call":
        V = np.maximum(price - K, 0)
    else:
        V = np.maximum(K - price, 0)

    # Setup matrices, A*V^{n+1} = B*V^n
    # Problem: How do we get the current price? 
    # Solution:
    j_idx = np.arange(1, M)    # Gives [1, 2, ..., M-1]
    S_j = price[j_idx]          # Vector of all relevan prices between S_low and S_high

    # We construct help variables as we did in task 1a) 
    alpha = -0.5 * dt * theta * (sigma**2 * S_j**2 / dS**2 - r * S_j / dS)
    beta = 1 + dt * theta * (sigma**2 * S_j**2 / dS**2 + r)
    gamma = -0.5 * dt * theta * (sigma**2 * S_j**2 / dS**2 + r * S_j / dS)
    
    # For matrix B (right-hand side)
    delta = 0.5 * dt * (1-theta) * (sigma**2 * S_j**2 / dS**2 - r * S_j / dS)
    epsilon = 1 - dt * (1-theta) * (sigma**2 * S_j**2 / dS**2 + r)
    zeta = 0.5 * dt * (1-theta) * (sigma**2 * S_j**2 / dS**2 + r * S_j / dS)
    
    A = np.zeros((M-1, M-1))
    B = np.zeros((M-1, M-1))

    for k in range(M-1):
        if k > 0:
            A[k, k-1] = alpha[k]
            B[k, k-1] = delta[k]
        A[k, k]   = beta[k]
        B[k, k]   = epsilon[k]
        if k < M-2:
            A[k, k+1] = gamma[k]
            B[k, k+1] = zeta[k]

    # Now that we have constructed A and B
    
    # Solve sytem from BC to t_0


    for n in range(N-1, -1, -1):
        # Boundary conditions
        if option == "call":
            V[0] = 0
            V[-1] = price[-1] - K * np.exp(-r*(T-time[n]))
        else:
            V[0]  = K * np.exp(-r*(T-time[n]))
            V[-1] = 0

        #   Update solution
        V_internal = V[1:M].copy()  # internal points at current time
        d = B @ V_internal
        
        # Solve AV^{n+1} = d for internal points
        V_next_internal = thomas_algorithm(A, d)

        V[1:M] = V_next_internal



    return price, V


#   Constants
K = 100
S_low = 0
S_high = 200
T = 1.0  # 1 year
N = 100
M = 100
r = 0.01
sigma = 0.2



# Price the option
# price - stock price
# V     - option price as function of stock price
price, V = finite_differences(S_low, S_high, T, N, M, K, r, sigma, option="call", theta=0.5)

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(price, V)
plt.xlabel('Stock Price')
plt.ylabel('Option Price')
plt.title('Call Option Price vs. Stock Price')
plt.grid(True)
plt.show()


print(V)

print("Done.")