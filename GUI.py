import streamlit as st
import math
from isentropic_flow import *
from normal_shock import *
from oblique_shock import *

# Define handler functions that write calculation results into dictionaries
def isentropic_flow_handler(M, GAMMA, result : dict):
    result["M"] = M
    result["M angle (deg)"] = math.degrees(mach_angle(M)) if M >= 1 else None
    result["PM angle (deg)"] = math.degrees(prandtl_meyer(M, GAMMA)) if M >= 1 else None
    result["T0/T"], result["p0/p"], result["rho0/rho"] = isentropic(M1=M, M2=0, GAMMA=GAMMA)
    result["T*/T"], result["p*/p"], result["rho*/rho"] = isentropic(M1=M, M2=1, GAMMA=GAMMA)
    result["A/A*"] = A_Ax(M, GAMMA)

def normal_shock_handler(M1, GAMMA, result : dict):
    result["M1"] = M1
    result["M2"], result["p2/p1"], result["T2/T1"], result["rho2/rho1"], result["p02/p1"], result["p02/p01"] = normalshock(M1, GAMMA)

def oblique_shock_handler(M1, delta_in, beta_in, Mn1_in, weak_flag_in, GAMMA, result : dict):
    if delta_in is not None:
        weak_flag, M2, delta, beta, p2_p1, T2_T1, rho2_rho1, p02_p01, Mn1, Mn2 = obliqueshock_delta(M1, delta_in, weak_flag=weak_flag_in, GAMMA=GAMMA)
    elif beta_in is not None:
        weak_flag, M2, delta, beta, p2_p1, T2_T1, rho2_rho1, p02_p01, Mn1, Mn2 = obliqueshock_beta(M1, beta_in, GAMMA)
    else:
        weak_flag, M2, delta, beta, p2_p1, T2_T1, rho2_rho1, p02_p01, Mn1, Mn2 = obliqueshock_M1n(M1, Mn1_in, GAMMA)
    result["M1"] = M1
    result["M2"] = M2
    result["beta (deg)"] = math.degrees(beta)
    result["delta (deg)"] = math.degrees(delta)
    result["p2/p1"] = p2_p1
    result["T2/T1"] = T2_T1
    result["rho2/rho1"] = rho2_rho1
    result["p02/p01"] = p02_p01
    result["Mn1"] = Mn1
    result["Mn2"] = Mn2
    result["type"] = "Weak" if weak_flag else "Strong"

    
# GUI for isentropic flow calculations
def isentropic_gui():
    st.header("Isentropic Flow Calculator")

    gamma_i = st.number_input("GAMMA = ", value=1.4, key="isentropic_gamma", format="%.2f") # Specific heat ratio input

    options = ["M", "p0/p", "T0/T", "rho0/rho", "A/A* (subsonic)", "A/A* (supersonic)", "Mach Angle (deg)", "Prandtl-Meyer Expansion Angle (deg)"]
    x = st.selectbox("Input:", options) # Selection box for type of input

    keys = ["M",
            "M angle (deg)", "PM angle (deg)",
            "p0/p", "T0/T", "rho0/rho",
            "p*/p", "T*/T", "rho*/rho",
            "A/A*"] # Result parameter keys
    result = dict.fromkeys(keys)    # Initialize empty result directory
    M = None                        # Initialize Mach number

    # Control block to handle all kinds of input
    if x == "M":
        M = st.number_input("M = ", value=2.0, format="%.4f")
        if M <= 0:
            st.error("Mach number must be greater than 0")
            return
    elif x == "p0/p":
        p0_p = st.number_input("p0/p = ", value=2.0, format="%.4f")
        if p0_p <= 1:
            st.error("p0/p must be greater than 1")
            return
        else:
            M = isentropic_solve_p0_p(p0_p, gamma_i)
    elif x == "T0/T":
        T0_T = st.number_input("T0/T = ", value=2.0, format="%.4f")
        if T0_T <= 1:
            st.error("T0/T must be greater than 1")
            return
        else:
            M = isentropic_solve_T0_T(T0_T, gamma_i)
    elif x == "rho0/rho":
        rho0_rho = st.number_input("rho0/rho = ", value=2.0, format="%.4f")
        if rho0_rho <= 1:
            st.error("rho0/rho must be greater than 1")
            return
        else:
            M = isentropic_solve_rho0_rho(rho0_rho, gamma_i)
    elif x == "A/A* (subsonic)":
        A_Ax_in = st.number_input("A/A* = ", value=2.0, format="%.4f")
        if A_Ax_in <= 1:
            st.error("A/A* must be greater than 1")
            return
        else:
            M = isentropic_solve_A_Ax(A_Ax_in, gamma_i, sub_flag=True)
    elif x == "A/A* (supersonic)":
        A_Ax_in = st.number_input("A/A* = ", value=2.0, format="%.4f")
        if A_Ax_in <= 1:
            st.error("A/A* must be greater than 1")
            return
        else:
            M = isentropic_solve_A_Ax(A_Ax_in, gamma_i, sub_flag=False)
    elif x == "Mach Angle (deg)":
        m_angle = st.number_input("Mach Angle (deg) = ", value=30.0, format="%.4f")
        if m_angle <= 0 or m_angle >= 90:
            st.error("Mach angle must be between 0 and 90 degrees")
            return
        else:
            m_angle = math.radians(m_angle)
            M = isentropic_solve_mach_angle(m_angle)
    elif x == "Prandtl-Meyer Expansion Angle (deg)":
        pm_angle = st.number_input("Prandtl-Meyer Expansion Angle (deg) = ", value=30.0, format="%.4f")
        pm_angle_max = 90*(math.sqrt((gamma_i + 1)/(gamma_i - 1)) - 1)
        if pm_angle <= 0 or pm_angle >= pm_angle_max:
            st.error(f"Prandtl-Meyer expansion angle must be between 0 and {pm_angle_max:.4f} degrees")
            return
        else:
            pm_angle = math.radians(pm_angle)
            M = isentropic_solve_prandtl_meyer(pm_angle, gamma_i)

    # Output block
    if M is not None:
        isentropic_flow_handler(M, gamma_i, result)
        st.subheader("Results")
        # Display Mach number and angles as metrics
        col1, col2, col3 = st.columns(3)
        col1.metric(label="$M$", value=f"{result['M']:.4f}")
        col2.metric(label="Mach Angle $\\mu$ (deg)", value=f"{result['M angle (deg)']:.4f}" if result['M angle (deg)'] is not None else "-")
        col3.metric(label="Prandtl-Meyer $\\nu$ (deg)", value=f"{result['PM angle (deg)']:.4f}" if result['PM angle (deg)'] is not None else "-")

        st.markdown("---")
        # Display isentropic ratios with LaTeX
        st.latex(r"\text{Isentropic Ratios}")
        col4, col5, col6 = st.columns(3)
        col4.latex(r"\frac{T_0}{T} = " + f"{result['T0/T']:.4f}")
        col5.latex(r"\frac{p_0}{p} = " + f"{result['p0/p']:.4f}")
        col6.latex(r"\frac{{\rho_0}}{{\rho}} = " + f"{result['rho0/rho']:.4f}")

        col7, col8, col9 = st.columns(3)
        col7.latex(r"\frac{T^*}{T} = " + f"{result['T*/T']:.4f}")
        col8.latex(r"\frac{p^*}{p} = " + f"{result['p*/p']:.4f}")
        col9.latex(r"\frac{{\rho^*}}{{\rho}} = " + f"{result['rho*/rho']:.4f}")
        st.latex(r"\frac{A}{A^*} = " + f"{result['A/A*']:.4f}")


# GUI for normal shock calculations
def normal_shock_gui():
    st.header("Normal Shock Calculator")

    gamma_n = st.number_input("GAMMA = ", value=1.4, key="normal_shock_gamma", format="%.2f")   # Specific heat ratio input

    options = ["M1", "M2", "p2/p1", "T2/T1", "rho2/rho1", "p02/p1", "p02/p01"]
    x = st.selectbox("Input:", options) # Selection box for type of input

    keys = ["M1", "M2",
            "p2/p1", "T2/T1", "rho2/rho1",
            "p02/p1", "p02/p01"]    # Result parameter keys
    result = dict.fromkeys(keys)    # Create empty result dictionary
    M1 = None                       # Initialize upstream Mach number

    # Control block to handle all types of input
    if x == "M1":
        M1 = st.number_input("M1 = ", value=2.0, format="%.4f")
        if M1 <= 1:
            st.error("Upstream Mach number must be greater than 1")
            return
    elif x == "M2":
        M2_min = math.sqrt((gamma_n - 1)/(2*gamma_n))
        M2 = st.number_input("M2 = ", value=0.5, format="%.4f")
        if M2 <= M2_min or M2 >= 1:
            st.error(f"Downstream Mach number must be between {M2_min} and 1")
            return
        else:
            M1 = other_Mach(M2, gamma_n)
    elif x == "p2/p1":
        p2_p1 = st.number_input("p2/p1 = ", value=2.0, format="%.4f")
        if p2_p1 <= 1:
            st.error("p2/p1 must be greater than 1")
            return
        else:
            M1 = M1_p2_p1(p2_p1, gamma_n)
    elif x == "T2/T1":
        T2_T1 = st.number_input("T2/T1 = ", value=2.0, format="%.4f")
        if T2_T1 <= 1:
            st.error("T2/T1 must be greater than 1")
            return
        else:
            M1 = M1_T2_T1(T2_T1, gamma_n)
    elif x == "rho2/rho1":
        rho2_rho1_max = round(((gamma_n + 1)/(gamma_n - 1)), 5)
        rho2_rho1 = st.number_input("rho2/rho1 = ", value = 1.5, format="%.4f")
        if rho2_rho1 <= 1 or rho2_rho1 >= rho2_rho1_max:
            st.error(f"rho2/rho1 must be between 1 and {rho2_rho1_max}")
            return
        else:
            M1 = M1_rho2_rho1(rho2_rho1, gamma_n)
    elif x == "p02/p1":
        p02_p1_min = round(((gamma_n + 1)/2)**(gamma_n/(gamma_n - 1)), 5)
        p02_p1 = st.number_input("p02/p1 = ", value = 5.0, format="%.4f")
        if p02_p1 <= p02_p1_min:
            st.error(f"p02/p1 must be greater than {p02_p1_min}")
            return
        else:
            M1 = M1_p02_p1(p02_p1, gamma_n)
    elif x == "p02/p01":
        p02_p01 = st.number_input("p02/p01 = ", value = 0.5, format="%.4f")
        if p02_p01 <= 0 or p02_p01 >= 1:
            st.error("p02/p01 must be between 0 and 1")
            return
        else:
            M1 = M1_p02_p01(p02_p01, gamma_n)

    # Calculate and display result
    if M1 is not None:
        normal_shock_handler(M1, gamma_n, result)
        st.subheader("Results")
        # Mach numbers
        col1, col2 = st.columns(2)
        col1.metric(label="$M_1$ (Upstream)", value=f"{result['M1']:.4f}")
        col2.metric(label="$M_2$ (Downstream)", value=f"{result['M2']:.4f}")

        st.markdown("---")
        # Static and stagnation ratios
        st.latex(r"\text{Normal Shock Ratios}")
        col3, col4, col5 = st.columns(3)
        col3.latex(r"\frac{p_2}{p_1} = " + f"{result['p2/p1']:.4f}")
        col4.latex(r"\frac{T_2}{T_1} = " + f"{result['T2/T1']:.4f}")
        col5.latex(r"\frac{\rho_2}{\rho_1} = " + f"{result['rho2/rho1']:.4f}")

        col6, col7 = st.columns(2)
        col6.latex(r"\frac{p_{02}}{p_1} = " + f"{result['p02/p1']:.4f}")
        col7.latex(r"\frac{p_{02}}{p_{01}} = " + f"{result['p02/p01']:.4f}")

    
# GUI for oblique shock calculations
def oblique_shock_gui():
    st.header("Oblique Shock Calculator")
    gamma_o = st.number_input("GAMMA = ", value=1.4, key="oblique_shock_gamma", format="%.2f")  # Specific ratio input

    # Upstream Mach number input and handling
    M1 = st.number_input("M1 = ", value=5.0, format="%.4f")
    if M1 <= 1:
        st.error("Upstream Mach number must be greater than 1")
        return
    
    options = ["Turn angle (weak) (deg)", "Turn angle (strong) (deg)", "Wave angle (deg)", "M1n"]
    x = st.selectbox("Input:", options) # Selection box for type of input

    keys = ["type",
            "M2", "delta (deg)", "beta (deg)",
            "p2/p1", "T2/T1", "rho2/rho1", "p02/p01",
            "Mn1", "Mn2"]   # Result parameter keys
    result = dict.fromkeys(keys)    # Create empty result directory

    # Control block to handle all types of input
    if x == "Turn angle (weak) (deg)":
        delta = st.number_input("Delta (deg) = ", value=10.0, format="%.4f")
        _, delta_max = delta_beta_max(M1, gamma_o)
        if delta <= 0 or delta >= 90:
            st.error("Turn angle out of bounds")
            return
        elif delta >= math.degrees(delta_max):
            st.warning(f"Shock Detached")
            return
        else:
            delta = math.radians(delta)
            weak_flag = True
            oblique_shock_handler(M1, delta, None, None, weak_flag, gamma_o, result)
    elif x == "Turn angle (strong) (deg)":
        delta = st.number_input("Delta (deg) = ", value=10.0, format="%.4f")
        _, delta_max = delta_beta_max(M1, gamma_o)
        if delta <= 0 or delta >= 90:
            st.error("Turn angle out of bounds")
            return
        elif delta >= math.degrees(delta_max):
            st.warning(f"Shock Detached")
            return
        else:
            delta = math.radians(delta)
            weak_flag = False
            oblique_shock_handler(M1, delta, None, None, weak_flag, gamma_o, result)
    elif x == "Wave angle (deg)":
        beta = st.number_input("Beta (deg) = ", value=30.0, format="%.4f")
        beta_min = math.asin(1/M1)
        if beta <= math.degrees(beta_min) or beta >= 90:
            st.error(f"Wave angle must be between {math.degrees(beta_min):.4f} and 90 degrees")
            return
        else:
            beta = math.radians(beta)
            oblique_shock_handler(M1, None, beta, None, None, gamma_o, result)
    elif x == "M1n":
        Mn1 = st.number_input("M1n = ", value=2.0, format="%.4f")
        if Mn1 <= 1 or Mn1 >= M1:
            st.error(f"M1n must be between 1 and {M1}")
            return
        else:
            oblique_shock_handler(M1, None, None, Mn1, None, gamma_o, result)
    
    # Display results
    if result["M2"] is not None:
        st.subheader("Results")
        # Upstream Mach numbers
        col1, col2 = st.columns(2)
        col1.metric(label="$M_1$ (Upstream)", value=f"{result['M1']:.4f}")
        col2.metric(label="$M_{1n}$ (Upstream Normal)", value=f"{result['Mn1']:.4f}")
        # Downstream Mach numbers
        col3, col4 = st.columns(2)
        col3.metric(label="$M_2$ (Downstream)", value=f"{result['M2']:.4f}")
        col4.metric(label="$M_{2n}$ (Downstream Normal)", value=f"{result['Mn2']:.4f}")
        st.markdown("---")
        # Shock type and deflection angle
        col5, col6, col7 = st.columns(3)
        col5.metric(label="Shock Type", value=f"{result['type']}")
        col6.metric(label="Deflection Angle $\\delta$ (deg)", value=f"{result['delta (deg)']:.4f}")
        col7.metric(label="Wave Angle $\\beta$ (deg)", value=f"{result['beta (deg)']:.4f}")
        st.markdown("---")
        # Static and stagnation ratios
        st.latex(r"\text{Oblique Shock Ratios}")
        col8, col9, col10, col11 = st.columns(4)
        col8.latex(r"\frac{p_2}{p_1} = " + f"{result['p2/p1']:.4f}")
        col9.latex(r"\frac{T_2}{T_1} = " + f"{result['T2/T1']:.4f}")
        col10.latex(r"\frac{\rho_2}{\rho_1} = " + f"{result['rho2/rho1']:.4f}")
        col11.latex(r"\frac{p_{02}}{p_{01}} = " + f"{result['p02/p01']:.4f}")


st.title("Gas Table Calculator")
st.markdown("")

isentropic_gui()
st.markdown("---")

normal_shock_gui()
st.markdown("---")

oblique_shock_gui()
st.markdown("---")