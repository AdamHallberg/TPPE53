import numpy as np
from scipy.stats import norm
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


def bsm_analytical(S, K, T, r, sigma, option_type="call"):
    """
    Black-Scholes-Merton analytical option pricing formula.
    
    Parameters:
    S: current stock price or array of prices
    K: strike price
    T: time to expiration (years)
    r: risk-free rate
    sigma: volatility
    option_type: "call" or "put"
    
    Returns:
    option_price: analytical BSM price(s)
    """
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    
    if option_type == "call":
        option_price = S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    else:
        option_price = K * np.exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)
    
    return option_price


def finite_differences(S_low, S_high, T, N, M, K, r, sigma, option="call", theta = 0.5):
    """
    INPUT: 
        ...

    OUTPUT: 
        Vector of price of option for different S_0 for inputed K.

    NOTES:
        theta = 0.5: ( Default )
        This is the Crank-Nicholson scheme. Crank-Nicholasson is unconditionally stable,
        so we do not need to consider the grid spaceing because of stability constraints.

        theta = 1: (Fullt implicit)
        This is the Fully-implicit scheme. This is also unconditionally stable, so we do
        not need to consider the grid spaceing because of stability constraints.

        theta = 0: ( Fully explicit)
        This is the Fully-explicit scheme. This is conditionally stable so we need to 
        choose Delta t leq fraq{(Delta S)^2}{(sigma S_{max})^2 + r Delta S}

        Other theta:
        I dont think there is a general solution for other thetas so this might or 
        might not work for other thetas.

    """

    #   Discritization 
    dt      = T / N
    dS      = (S_high - S_low) / M
    time    = np.linspace(0, T, N+1)
    price   = np.linspace(S_low, S_high, M+1)

    # Setup terminal condition, from 1c)
    if option == "call":
        V = np.maximum(price - K, 0)
    else:
        V = np.maximum(K - price, 0)
    

    # Setup matrices, A*V^{n+1} = B*V^n
    j_idx = np.arange(1, M)    # Gives [1, 2, ..., M-1]
    S_j = price[j_idx]         # Vector of all relevant prices between S_low and S_high

    # We construct help variables as we did in task 1a) 
    
    # For matrix A
    alpha = -theta * (r * S_j) / (2*dS) \
            + theta * 0.5 * (sigma**2 * S_j**2) / (dS**2)

    beta  = 1/dt \
            - theta * (sigma**2 * S_j**2) / (dS**2) \
            - theta * r   

    gamma = theta * (r * S_j) / (2*dS) \
            + theta * 0.5 * (sigma**2 * S_j**2) / (dS**2)

    # For matrix B  
    delta = (1-theta) * (r * S_j) / (2*dS) \
            - (1-theta) * 0.5 * (sigma**2 * S_j**2) / (dS**2)

    epsilon = 1/dt \
              + (1-theta) * (sigma**2 * S_j**2) / (dS**2) \
              + (1-theta) * r

    zeta = -(1-theta) * (r * S_j) / (2*dS) \
           - (1-theta) * 0.5 * (sigma**2 * S_j**2) / (dS**2)

        
    A = np.zeros((M-1, M-1))
    B = np.zeros((M-1, M-1))

    # This might be possible to do differently
    for k in range(M-1):
        if k > 0:
            A[k, k-1] = alpha[k]
            B[k, k-1] = delta[k]
        A[k, k]   = beta[k]
        B[k, k]   = epsilon[k]
        if k < M-2:
            A[k, k+1] = gamma[k]
            B[k, k+1] = zeta[k]

    #   Now that we have constructed A and B
    #   we can solve, V^n = inv(B)*A*V^{n+1}
    #   iteratively. With Boundary conditions
    #   at S_low and S_high (V[0] and V[-1])
    #   according to 1c).

    for n in range(N-1, -1, -1):
        # Boundary conditions
        current_time = time[n]
        time_to_maturity = T - current_time
    
        # BC are different depending on put or call.
        if option == "call":
            V[0] = 0  
            V[-1] = S_high - K * np.exp(-r * time_to_maturity) 
        else:
            V[0] = K * np.exp(-r * time_to_maturity)  
            V[-1] = 0  

        #   Update solution
        V_internal = V[1:M].copy()  # internal points at current time
        d = A @ V_internal
                
        d[0]  += alpha[0] * V[0]
        d[-1] += gamma[-1] * V[-1]

        # Solve V^{n} = inv(B)*A*V^{n+1} 
        #V_next_internal = thomas_algorithm(B, d)
        V_next_internal = np.linalg.solve(B, d)
        
        V[1:M] = V_next_internal



    return price, V


#   Constants
K = 100
T = 1.0     # Enhet?? År eler?
N = 100
M = 100
r = 0.01    # Enhet?
sigma = 0.2 # Enhet?
option="call"

# From 1.b)
S_low = 51.268
S_high = 191.191


# Price the option
# price - stock price
# V     - option price as function of stock price
price, V = finite_differences(S_low, S_high, T, N, M, K, r, sigma, option, theta=0.5)


if False:
    # Plot results
    plt.figure(figsize=(10, 6))
    plt.plot(price, V)
    plt.xlabel('Stock Price')
    plt.ylabel('Option Price')
    plt.title(f"{option} Option Price vs. Stock Price for K = {K} and T = {T}[y]")
    plt.grid(True)
    plt.show()

    # Wierd results for the endpoint? 


price   = np.linspace(S_low, S_high, M+1)
analytical_results = bsm_analytical(price, K, T, r, sigma, option_type="call")

if True:
    # Compare our prices and BSM prices
    plt.figure(figsize=(10, 6))
    plt.plot(price, V, 'b-', linewidth=2, label='Finite Difference')
    plt.plot(price, analytical_results, 'r--', linewidth=2, label='Analytical BSM')
    plt.legend(['FD-method','analytical solution'])
    plt.xlabel('Stock Price')
    plt.ylabel('Option Price')
    plt.title(f'{option} Option Price Comparison (K={K}, T={T} [y], σ={sigma}, r={r})', fontsize=14)    
    plt.grid(True)
    plt.show()

    # Wierd results for the endpoint? 
