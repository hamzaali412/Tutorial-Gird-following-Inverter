import math
from pyomo.environ import ConcreteModel, Param, Var, Expression, cos, sin

from helpers import STATE_KEYS
from controllers import pll_residuals, power_controller_residuals, current_controller_residuals


def build_reusable_step_model(p):
    m = ConcreteModel()

    m.t_now = Param(initialize=0.0, mutable=True)
    m.theta_grid_now = Param(initialize=0.0, mutable=True)
    m.Vmag_now = Param(initialize=1.025, mutable=True)
    m.prev = Param(STATE_KEYS, initialize=0.0, mutable=True)

    m.V_inv = Var(initialize=p["v_inv_0"], bounds=(0.5, 1.5))
    m.theta_inv = Var(initialize=p["theta_inv_0"], bounds=(-1.5, 1.5))

    m.theta_pll = Var(initialize=p["theta_inv_0"])
    m.vq_fil = Var(initialize=0.0)
    m.I_pll = Var(initialize=0.0)
    m.delta_pll = Var(initialize=0.0)

    m.P_fil = Var(initialize=p["Pref_const"])
    m.Q_fil = Var(initialize=p["Qref_const"])
    m.I_pqd = Var(initialize=0.0)
    m.I_pqq = Var(initialize=0.0)
    m.Id_ref = Var(initialize=0.0)
    m.Iq_ref = Var(initialize=0.0)

    m.Id = Var(initialize=0.0)
    m.Iq = Var(initialize=0.0)
    m.Id_fil = Var(initialize=0.0)
    m.Iq_fil = Var(initialize=0.0)
    m.I_ccd = Var(initialize=0.0)
    m.I_ccq = Var(initialize=0.0)
    m.err_d = Var(initialize=0.0)
    m.err_q = Var(initialize=0.0)
    m.e_ac_d = Var(initialize=0.0)
    m.e_ac_q = Var(initialize=0.0)

    m.va = Expression(expr=m.V_inv * cos(p["omega"] * m.t_now + m.theta_inv))
    m.vb = Expression(expr=m.V_inv * cos(p["omega"] * m.t_now + m.theta_inv - 2 * math.pi / 3))
    m.vc = Expression(expr=m.V_inv * cos(p["omega"] * m.t_now + m.theta_inv + 2 * math.pi / 3))

    m.omega_pll_angle = Expression(expr=(p["omega"] * m.t_now) + m.theta_pll)

    m.vd = Expression(expr=(2 / 3) * (
        m.va * cos(m.omega_pll_angle)
        + m.vb * cos(m.omega_pll_angle - 2 * math.pi / 3)
        + m.vc * cos(m.omega_pll_angle + 2 * math.pi / 3)
    ))

    m.vq = Expression(expr=-(2 / 3) * (
        m.va * sin(m.omega_pll_angle)
        + m.vb * sin(m.omega_pll_angle - 2 * math.pi / 3)
        + m.vc * sin(m.omega_pll_angle + 2 * math.pi / 3)
    ))

    m.P_expr = Expression(expr=(3 / 2) * (m.vd * m.Id + m.vq * m.Iq))
    m.Q_expr = Expression(expr=(-3 / 2) * (m.vq * m.Id - m.vd * m.Iq))
    m.omega_pll_feedback = Expression(expr=p["omega"] * (1 + m.delta_pll))

    delta = m.theta_inv - m.theta_grid_now
    m.P_flow = Expression(
        expr=m.V_inv**2 * p["G_pu"] - m.V_inv * m.Vmag_now * (p["G_pu"] * cos(delta) + p["B_pu"] * sin(delta))
    )
    m.Q_flow = Expression(
        expr=-(m.V_inv**2 * p["B_pu"]) - m.V_inv * m.Vmag_now * (p["G_pu"] * sin(delta) - p["B_pu"] * cos(delta))
    )

    residual_map = {}
    residual_map.update(pll_residuals(m, p))
    residual_map.update(power_controller_residuals(m, p))
    residual_map.update(current_controller_residuals(m, p))

    residual_order = [
        "delta_pll_eq", "theta_pll_eq", "power_flow_P_eq", "power_flow_Q_eq",
        "pll_vq_fil_eq", "I_pll_eq", "P_fil_eq", "Q_fil_eq",
        "I_pqd_eq", "I_pqq_eq", "Id_ref_eq", "Iq_ref_eq",
        "Id_fil_eq", "Iq_fil_eq", "I_ccd_eq", "I_ccq_eq",
        "err_d_eq", "err_q_eq", "e_ac_d_eq", "e_ac_q_eq",
        "Id_eq", "Iq_eq"
    ]

    vars_order = [
        m.V_inv, m.theta_inv, m.theta_pll, m.delta_pll,
        m.vq_fil, m.I_pll, m.P_fil, m.Q_fil, m.I_pqd, m.I_pqq,
        m.Id_ref, m.Iq_ref, m.Id, m.Iq, m.Id_fil, m.Iq_fil,
        m.I_ccd, m.I_ccq, m.err_d, m.err_q, m.e_ac_d, m.e_ac_q
    ]

    residual_exprs = [residual_map[k] for k in residual_order]
    return m, vars_order, residual_exprs