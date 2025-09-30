import math
from scipy.optimize import newton

# Calculating Mach angle and Prandtl-Meyer expansion angle
def mach_angle(M):
    if M < 1:
        return None
    else:
        return math.asin(1/M)

def prandtl_meyer(M, GAMMA):
    if M < 1:
        return None
    else:
        fG = (GAMMA + 1)/(GAMMA - 1)
        pm = math.sqrt(fG)*math.atan(math.sqrt((M**2 - 1)/fG)) - math.atan(math.sqrt(M**2 - 1))
        return pm

# Functions for isentropic flow calculations
def isentropic(M1, M2, GAMMA):
    T2_T1 = (2 + (GAMMA - 1)*M1**2)/(2 + (GAMMA - 1)*M2**2) # Temperature ratio
    p2_p1 = (T2_T1)**(GAMMA/(GAMMA - 1))                    # Pressure ratio
    rho2_rho1 = (T2_T1)**(1/(GAMMA - 1))                    # Density ratio

    return T2_T1, p2_p1, rho2_rho1

def A_Ax(M, GAMMA):
    A_Ax = (1/M)*((2/(GAMMA + 1))*(1 + ((GAMMA - 1)/2)*M**2))**((GAMMA + 1)/(2*(GAMMA - 1)))    # Area ratio
    return A_Ax

# Functions for inverse calculations
def isentropic_solve_p0_p(p0_p, GAMMA):
    M = math.sqrt((p0_p**((GAMMA - 1)/GAMMA) - 1)*2/(GAMMA - 1))
    return M

def isentropic_solve_T0_T(T0_T, GAMMA):
    M = math.sqrt((T0_T - 1)*2/(GAMMA - 1))
    return M

def isentropic_solve_rho0_rho(rho0_rho, GAMMA):
    M = math.sqrt((rho0_rho**(GAMMA - 1) - 1)*2/(GAMMA - 1))
    return M

def isentropic_solve_A_Ax(A_Ax_in, GAMMA, sub_flag : bool):
    f = lambda M : A_Ax_in - A_Ax(M, GAMMA) # Defining residual function
    M_guess = 0.01 if sub_flag else 20.0    # Set initial guess based on subsonic/supersonic flow
    M = newton(f, M_guess)                  # Using Newton-Raphson method to solve for M
    return M

def isentropic_solve_mach_angle(m_angle):
    M = 1/math.sin(m_angle)   # m_angle is in degrees
    return M

def isentropic_solve_prandtl_meyer(pm_angle, GAMMA):
    f = lambda M : pm_angle - prandtl_meyer(M, GAMMA) # Defining residual function
    M_guess = 2.0                                     # Initial guess
    M = newton(f, M_guess)                            # Using Newton-Raphson method to solve for M
    return M