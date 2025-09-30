import math
from scipy.optimize import newton

# Function to calculate downstream and upstream Mach numbers since they are related by the same equation
def other_Mach(M1, GAMMA):
    M2 = math.sqrt(((GAMMA - 1)*M1**2 + 2)/(2*GAMMA*M1**2 - (GAMMA - 1)))   # Downstream Mach number
    return M2

# Functions for inverse calculations
def M1_p2_p1(p2_p1, GAMMA):
    M2 = lambda M1 : other_Mach(M1, GAMMA)  # Upstream Mach number
    res = lambda M : p2_p1 - (1 + GAMMA*M**2)/(1 + GAMMA*M2(M)**2)  # Residual function
    M1_guess = 5    # Initial guess
    M1 = newton(res, M1_guess)  # Solving for M1 using Newton-Raphson method
    return M1

def M1_T2_T1(T2_T1, GAMMA):
    M2 = lambda M1 : other_Mach(M1, GAMMA)  # Upstream Mach number
    res = lambda M : T2_T1 - (2 + (GAMMA - 1)*M**2)/(2 + (GAMMA - 1)*M2(M)**2)  # Residual function
    M1_guess = 5    # Initial guess
    M1 = newton(res, M1_guess)  # Solving for M1 using Newton-Raphson method
    return M1

def M1_rho2_rho1(rho2_rho1, GAMMA):
    M2 = lambda M1 : other_Mach(M1, GAMMA)  # Upstream Mach number
    p2_p1 = lambda M : (1 + GAMMA*M**2)/(1 + GAMMA*M2(M)**2)    # Pressure ratio
    T2_T1 = lambda M : (2 + (GAMMA - 1)*M**2)/(2 + (GAMMA - 1)*M2(M)**2)    # Temperature ratio
    res = lambda M : rho2_rho1 - p2_p1(M)/T2_T1(M)  # Residual function
    M1_guess = 2    # Initial guess
    M1 = newton(res, M1_guess)  # Solving for M1 using Newton-Raphson method
    return M1

def M1_p02_p1(p02_p1, GAMMA):
    f = lambda M : ((1 + ((GAMMA - 1)/2)*other_Mach(M, GAMMA)**2)**(GAMMA/(GAMMA - 1))) \
        * (1 + GAMMA*M**2)/(1 + GAMMA*other_Mach(M, GAMMA)**2)    # Total-to-static pressure ratio (downstream to upstream static)
    res = lambda x : p02_p1 - f(x)  # Residual function
    M1_guess = 2    # Initial guess
    M1 = newton(res, M1_guess)  # Solving for M1 using Newton-Raphson method
    return M1

def M1_p02_p01(p02_p01, GAMMA):
    p02_p1 = lambda M : ((1 + ((GAMMA - 1)/2)*other_Mach(M, GAMMA)**2)**(GAMMA/(GAMMA - 1))) \
        * (1 + GAMMA*M**2)/(1 + GAMMA*other_Mach(M, GAMMA)**2)  # Total-to-static pressure ratio (downstream to upstream static)
    res = lambda x : p02_p01 - p02_p1(x)/((1 + ((GAMMA - 1)/2)*x**2)**(GAMMA/(GAMMA - 1)))  # Residual function
    M1_guess = 2    # Initial guess
    M1 = newton(res, M1_guess)  # Solving for M1 using Newton-Raphson method
    return M1

# Defining the main normalshock function
def normalshock(M1, GAMMA):
    M2 = other_Mach(M1, GAMMA)
    p2_p1 = (1 + GAMMA*M1**2)/(1 + GAMMA*M2**2)                         # Pressure ratio
    T2_T1 = (2 + (GAMMA - 1)*M1**2)/(2 + (GAMMA - 1)*M2**2)             # Temperature ratio
    rho2_rho1 = p2_p1/T2_T1                                             # Density ratio
    p02_p1 = ((1 + ((GAMMA - 1)/2)*M2**2)**(GAMMA/(GAMMA - 1)))*p2_p1   # Total-to-static pressure ratio (downstream to upstream static)
    p02_p01 = p02_p1/((1 + ((GAMMA - 1)/2)*M1**2)**(GAMMA/(GAMMA - 1))) # Total pressure ratio (downstream to upstream total)
    
    return M2, p2_p1, T2_T1, rho2_rho1, p02_p1, p02_p01