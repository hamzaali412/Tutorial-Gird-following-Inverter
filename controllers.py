from pyomo.environ import *


def pll_residuals(m, p):
    c1 = p["h"] / (2 * p["Tf_pll"] + p["h"])
    c2 = (2 * p["Tf_pll"] - p["h"]) / (2 * p["Tf_pll"] + p["h"])

    delta_prev = m.prev["delta_pll"] * p["omega"]
    delta_curr = m.delta_pll * p["omega"]

    return {
        "delta_pll_eq": m.delta_pll - (p["kp_pll"] * m.vq_fil + m.I_pll),
        "theta_pll_eq": m.theta_pll - (m.prev["theta_pll"] + (p["h"] / 2) * (delta_prev + delta_curr)),
        "pll_vq_fil_eq": m.vq_fil - (c1 * (m.vq + m.prev["vq"]) + c2 * m.prev["vq_fil"]),
        "I_pll_eq": m.I_pll - (
            m.prev["I_pll"] + (p["ki_pll"] * p["h"] / 2) * (m.vq_fil + m.prev["vq_fil"])
        ),
    }

def smax(a, b, eps=1e-7):
    return (a + b + sqrt((a - b)**2 + eps)) / 2

def smin(a, b, eps=1e-7):
    return (a + b - sqrt((a - b)**2 + eps)) / 2

def power_controller_residuals(m, p):
    c1pq = p["h"] / (2 * p["Tf_pq"] + p["h"])
    c2pq = (2 * p["Tf_pq"] - p["h"]) / (2 * p["Tf_pq"] + p["h"])
    mode = p.get("mode", "PQ").upper()
    p_sup = 0.0
    q_sup = 0.0
    if mode == "FS":
        f_dev = m.delta_pll 
        f_db = p.get("f_db", 0.0005)
        k_f = p.get("K_droop_f", 20.0)
        p_sup = (-k_f * smax(0, f_dev - f_db)) + (-k_f * smin(0, f_dev + f_db))

    if mode == "VOLT-VAR":
        v_target = p.get("V_target", 1.0686727498480264) 
        dV = m.V_inv - v_target
        v_db = p.get("Vdb", 0.01)
        k_v = p.get("K_droop_v", 20.0)
        q_max = p.get("Qmax_sup", 0.4)
        q_val = (-k_v * smax(0, dV - v_db)) + (-k_v * smin(0, dV + v_db))
        q_sup = smax(-q_max, smin(q_max, q_val))
    current_pref = p["Pref_const"] + p_sup
    current_qref = p["Qref_const"] + q_sup

    return {
        "power_flow_P_eq": m.P_expr - m.P_flow,
        "power_flow_Q_eq": m.Q_expr - m.Q_flow,
        "P_fil_eq": m.P_fil - (c1pq * (m.P_expr + m.prev["P_expr"]) + c2pq * m.prev["P_fil"]),
        "Q_fil_eq": m.Q_fil - (c1pq * (m.Q_expr + m.prev["Q_expr"]) + c2pq * m.prev["Q_fil"]),
        "I_pqd_eq": m.I_pqd - (
            m.prev["I_pqd"]
            + (p["ki_pq"] * p["h"] / 2) * ((current_pref - m.P_fil) + (current_pref - m.prev["P_fil"]))
        ),
        "I_pqq_eq": m.I_pqq - (
            m.prev["I_pqq"]
            + (p["ki_pq"] * p["h"] / 2) * ((current_qref - m.Q_fil) + (current_qref - m.prev["Q_fil"]))
        ),
        "Id_ref_eq": m.Id_ref - (p["Kp_pq"] * (current_pref - m.P_fil) + m.I_pqd),
        "Iq_ref_eq": m.Iq_ref - (p["Kp_pq"] * (current_qref - m.Q_fil) + m.I_pqq),
    }

# def power_controller_residuals(m, p):
#     c1pq = p["h"] / (2 * p["Tf_pq"] + p["h"])
#     c2pq = (2 * p["Tf_pq"] - p["h"]) / (2 * p["Tf_pq"] + p["h"])

#     return {
#         "power_flow_P_eq": m.P_expr - m.P_flow,
#         "power_flow_Q_eq": m.Q_expr - m.Q_flow,
#         "P_fil_eq": m.P_fil - (c1pq * (m.P_expr + m.prev["P_expr"]) + c2pq * m.prev["P_fil"]),
#         "Q_fil_eq": m.Q_fil - (c1pq * (m.Q_expr + m.prev["Q_expr"]) + c2pq * m.prev["Q_fil"]),
#         "I_pqd_eq": m.I_pqd - (
#             m.prev["I_pqd"]
#             + (p["ki_pq"] * p["h"] / 2) * ((p["Pref_const"] - m.P_fil) + (p["Pref_const"] - m.prev["P_fil"]))
#         ),
#         "I_pqq_eq": m.I_pqq - (
#             m.prev["I_pqq"]
#             + (p["ki_pq"] * p["h"] / 2) * ((p["Qref_const"] - m.Q_fil) + (p["Qref_const"] - m.prev["Q_fil"]))
#         ),
#         "Id_ref_eq": m.Id_ref - (p["Kp_pq"] * (p["Pref_const"] - m.P_fil) + m.I_pqd),
#         "Iq_ref_eq": m.Iq_ref - (p["Kp_pq"] * (p["Qref_const"] - m.Q_fil) + m.I_pqq),
#     }


def current_controller_residuals(m, p):
    c1ig = p["h"] / (2 * p["Tf_ig"] + p["h"])
    c2ig = (2 * p["Tf_ig"] - p["h"]) / (2 * p["Tf_ig"] + p["h"])

    rhs_prev_d = (
        -p["Rf_pu"] * m.prev["Id"]
        + m.prev["omega_pll_feedback"] * p["Lf_pu"] * m.prev["Iq"]
        + m.prev["e_ac_d"] - m.prev["vd"]
    )
    rhs_curr_d = (
        -p["Rf_pu"] * m.Id
        + m.omega_pll_feedback * p["Lf_pu"] * m.Iq
        + m.e_ac_d - m.vd
    )

    rhs_prev_q = (
        -p["Rf_pu"] * m.prev["Iq"]
        - m.prev["omega_pll_feedback"] * p["Lf_pu"] * m.prev["Id"]
        + m.prev["e_ac_q"] - m.prev["vq"]
    )
    rhs_curr_q = (
        -p["Rf_pu"] * m.Iq
        - m.omega_pll_feedback * p["Lf_pu"] * m.Id
        + m.e_ac_q - m.vq
    )

    return {
        "Id_fil_eq": m.Id_fil - (c1ig * (m.Id + m.prev["Id"]) + c2ig * m.prev["Id_fil"]),
        "Iq_fil_eq": m.Iq_fil - (c1ig * (m.Iq + m.prev["Iq"]) + c2ig * m.prev["Iq_fil"]),
        "I_ccd_eq": m.I_ccd - (
            m.prev["I_ccd"]
            + (p["ki_ig"] * p["h"] / 2) * ((m.Id_ref - m.Id_fil) + (m.prev["Id_ref"] - m.prev["Id_fil"]))
        ),
        "I_ccq_eq": m.I_ccq - (
            m.prev["I_ccq"]
            + (p["ki_ig"] * p["h"] / 2) * ((m.Iq_ref - m.Iq_fil) + (m.prev["Iq_ref"] - m.prev["Iq_fil"]))
        ),
        "err_d_eq": m.err_d - (p["kp_ig"] * (m.Id_ref - m.Id_fil) + m.I_ccd),
        "err_q_eq": m.err_q - (p["kp_ig"] * (m.Iq_ref - m.Iq_fil) + m.I_ccq),
        "e_ac_d_eq": m.e_ac_d - (m.err_d - m.omega_pll_feedback * m.Iq * p["Lf_pu"] + m.vd),
        "e_ac_q_eq": m.e_ac_q - (m.err_q + m.omega_pll_feedback * m.Id * p["Lf_pu"] + m.vq),
        "Id_eq": m.Id - (m.prev["Id"] + (p["h"] / 2) * (rhs_prev_d + rhs_curr_d) / p["Lf_pu"]),
        "Iq_eq": m.Iq - (m.prev["Iq"] + (p["h"] / 2) * (rhs_prev_q + rhs_curr_q) / p["Lf_pu"]),
    }