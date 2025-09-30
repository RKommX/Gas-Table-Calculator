import math
import numpy as np

# Function to calculate the flow deflection angle given upstream Mach number and shock wave angle
def delta_cal(M1, beta, GAMMA):
    if beta == math.asin(1/M1) or beta == math.pi/2:    # Handle special cases to avoid divergence
        return 0
    X = (2 * (1 / math.tan(beta)) * ((M1 * math.sin(beta))**2 - 1)) / (2 + M1**2 *(GAMMA + math.cos(2 * beta)))
    return math.atan(X)

# Function to calculate the maximum deflection angle and corresponding shock wave angle for a given upstream Mach number
def delta_beta_max(M1, GAMMA):
    beta_min = math.asin(1/M1)  # Minimum wave angle
    beta_max = math.pi/2        # Maximum wave angle
    
    # Create the grid of betas and deltas, then get an initial guess of the beta that gives max delta
    betas = np.linspace(beta_min, beta_max, 1000)
    deltas = []
    for beta in betas:
        deltas.append(delta_cal(M1, beta, GAMMA))
    deltas = np.array(deltas)
    beta_guess = betas[np.argmax(deltas)]

    # Take a neighbourhood about beta_guess
    bound = 1e-3
    beta1 = beta_guess - bound
    beta4 = beta_guess + bound

    psi = (math.sqrt(5) - 1)/2  # Golden ratio
    tol = 1e-12                 # Convergence tolerance
    
    # Golden section search loop
    while(beta4 - beta1 < tol):
        beta2 = beta1 + psi * (beta4 - beta1)
        beta3 = beta4 - psi * (beta4 - beta1)
        
        if delta_cal(M1, beta2, GAMMA) > delta_cal(M1, beta3, GAMMA):
            beta4 = beta3
        else:
            beta1 = beta2

    beta_delta_max = (beta1 + beta4)/2  # Wave angle for maximum deflection
    delta_max = delta_cal(M1, beta_delta_max, GAMMA)    # Maximum deflection
    return beta_delta_max, delta_max

# Function to calculate shock parameters given upstream Mach number, turn angle and solution type
def obliqueshock_delta(M1, delta, weak_flag : bool, GAMMA):
    beta_max, _ = delta_beta_max(M1, GAMMA)            # Wave angle at maximum turn angle

    # Define bounds for wave angle depending on solution branch
    beta_a = math.asin(1/M1) if weak_flag else beta_max
    beta_b = beta_max if weak_flag else math.pi/2
    convergence = 1e-12
    
    beta = 0.0          # Initialising wave angle

    # Solve for beta via bisection method
    while True:
        beta_c = (beta_a + beta_b)/2
        re = delta - delta_cal(M1, beta_c, GAMMA)
        
        if abs(re) <= convergence:
            beta = beta_c
            break
        
        elif re < 0:
            if weak_flag:
                beta_b = beta_c
            else:
                beta_a = beta_c
        else:
            if weak_flag:
                beta_a = beta_c
            else:
                beta_b = beta_c

    # Calculate results
    Mn1 = M1*math.sin(beta)
    Mn2 = math.sqrt(((GAMMA - 1)*Mn1**2 + 2)/(2*GAMMA*Mn1**2 - (GAMMA - 1)))
    M2 = Mn2/math.sin(beta - delta)
    p2_p1 = (1 + GAMMA*Mn1**2)/(1 + GAMMA*Mn2**2)
    T2_T1 = (2 + (GAMMA - 1)*Mn1**2)/(2 + (GAMMA - 1)*Mn2**2)
    rho2_rho1 = p2_p1/T2_T1
    p02_p1 = ((1 + ((GAMMA - 1)/2)*Mn2**2)**(GAMMA/(GAMMA - 1)))*p2_p1
    p02_p01 = p02_p1/((1 + ((GAMMA - 1)/2)*Mn1**2)**(GAMMA/(GAMMA - 1)))
    
    return weak_flag, M2, delta, beta, p2_p1, T2_T1, rho2_rho1, p02_p01, Mn1, Mn2

# Function to calculate shock parameters given upstream Mach number and wave angle
def obliqueshock_beta(M1, beta, GAMMA):
    delta = delta_cal(M1, beta, GAMMA)  # Calculate turn angle

    # Obtain type of oblique shock by comparing with wave angle for maximum deflection
    beta_max, _ = delta_beta_max(M1, GAMMA)
    weak_flag = beta <= beta_max

    # Calculate results
    Mn1 = M1*math.sin(beta)
    Mn2 = math.sqrt(((GAMMA - 1)*Mn1**2 + 2)/(2*GAMMA*Mn1**2 - (GAMMA - 1)))
    M2 = Mn2/math.sin(beta - delta)
    p2_p1 = (1 + GAMMA*Mn1**2)/(1 + GAMMA*Mn2**2)
    T2_T1 = (2 + (GAMMA - 1)*Mn1**2)/(2 + (GAMMA - 1)*Mn2**2)
    rho2_rho1 = p2_p1/T2_T1
    p02_p1 = ((1 + ((GAMMA - 1)/2)*Mn2**2)**(GAMMA/(GAMMA - 1)))*p2_p1
    p02_p01 = p02_p1/((1 + ((GAMMA - 1)/2)*Mn1**2)**(GAMMA/(GAMMA - 1)))

    return weak_flag, M2, delta, beta, p2_p1, T2_T1, rho2_rho1, p02_p01, Mn1, Mn2

# Function to calculate shock parameters given upstream Mach number and its normal component
def obliqueshock_M1n(M1, M1n, GAMMA):
    beta = math.asin(M1n/M1)    # Calculate wave angle
    return obliqueshock_beta(M1, beta, GAMMA)   # Use previous function to calculate results