from helpers import abc_from_vtheta, dq_from_abc
import math
def initialize_history(hist, p):
    # 1. Calculate what the physics demands at t=0
    delta0 = p["theta_inv_0"] - p["theta_grid_init"]
    P0_phys = (p["v_inv_0"]**2 * p["G_pu"] - 
               p["v_inv_0"] * p["Vmag_pu"] * (p["G_pu"] * math.cos(delta0) + p["B_pu"] * math.sin(delta0)))
    Q0_phys = (-(p["v_inv_0"]**2 * p["B_pu"]) - 
               p["v_inv_0"] * p["Vmag_pu"] * (p["G_pu"] * math.sin(delta0) - p["B_pu"] * math.cos(delta0)))
    
    # 2. Basic States
    hist["V_inv"][0] = p["v_inv_0"]
    hist["theta_inv"][0] = p["theta_inv_0"]
    hist["theta_pll"][0] = p["theta_inv_0"]
    hist["vq_fil"][0] = 0.0
    hist["I_pll"][0] = 0.0
    hist["delta_pll"][0] = 0.0
    hist["Id"][0] = (2 / 3 * P0_phys) / p["v_inv_0"]
    hist["Iq"][0] = (2 / 3 * Q0_phys) / p["v_inv_0"]
    hist["P_fil"][0] = P0_phys
    hist["Q_fil"][0] = Q0_phys
    hist["P_expr"][0] = P0_phys
    hist["Q_expr"][0] = Q0_phys
    hist["I_pqd"][0] = hist["Id"][0]
    hist["I_pqq"][0] = hist["Iq"][0]
    hist["Id_ref"][0] = hist["Id"][0]
    hist["Iq_ref"][0] = hist["Iq"][0]
    hist["Id_fil"][0] = hist["Id"][0]
    hist["Iq_fil"][0] = hist["Iq"][0]
    hist["I_ccd"][0] = p["Rf_pu"] * hist["Id"][0]
    hist["I_ccq"][0] = p["Rf_pu"] * hist["Iq"][0]
    hist["err_d"][0] = hist["I_ccd"][0]
    hist["err_q"][0] = hist["I_ccq"][0]
    va0, vb0, vc0 = abc_from_vtheta(hist["V_inv"][0], hist["theta_inv"][0], p["time_steps"][0], p["omega"])
    vd0, vq0 = dq_from_abc(va0, vb0, vc0, hist["theta_pll"][0], p["time_steps"][0], p["omega"])
    hist["vd"][0] = vd0
    hist["vq"][0] = vq0
    hist["omega_pll_feedback"][0] = p["omega"] * (1 + hist["delta_pll"][0])
    hist["e_ac_d"][0] = (hist["err_d"][0] - hist["omega_pll_feedback"][0] * hist["Iq"][0] * p["Lf_pu"] + vd0)
    hist["e_ac_q"][0] = (hist["err_q"][0] + hist["omega_pll_feedback"][0] * hist["Id"][0] * p["Lf_pu"] + vq0)

    va0, vb0, vc0 = abc_from_vtheta(
        hist["V_inv"][0],
        hist["theta_inv"][0],
        p["time_steps"][0],
        p["omega"],
    )
    vd0, vq0 = dq_from_abc(
        va0, vb0, vc0,
        hist["theta_pll"][0],
        p["time_steps"][0],
        p["omega"],
    )

    hist["vd"][0] = vd0
    hist["vq"][0] = vq0
    hist["P_expr"][0] = (3 / 2) * (vd0 * hist["Id"][0] + vq0 * hist["Iq"][0])
    hist["Q_expr"][0] = (-3 / 2) * (vq0 * hist["Id"][0] - vd0 * hist["Iq"][0])
    hist["omega_pll_feedback"][0] = p["omega"] * (1 + hist["delta_pll"][0])

    hist["e_ac_d"][0] = (
        hist["err_d"][0]
        - hist["omega_pll_feedback"][0] * hist["Iq"][0] * p["Lf_pu"]
        + vd0
    )
    hist["e_ac_q"][0] = (
        hist["err_q"][0]
        + hist["omega_pll_feedback"][0] * hist["Id"][0] * p["Lf_pu"]
        + vq0
    )