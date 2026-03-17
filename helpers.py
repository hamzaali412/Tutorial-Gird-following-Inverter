import math
import numpy as np
from pyomo.environ import *
from pyomo.environ import value
from pyomo.core.expr.calculus.derivatives import differentiate, Modes
import scipy.sparse as sp
import scipy.sparse.linalg as spla


STATE_KEYS = [
    "V_inv", "theta_inv", "theta_pll", "vq_fil", "I_pll", "delta_pll",
    "P_fil", "Q_fil", "I_pqd", "I_pqq", "Id_ref", "Iq_ref",
    "Id", "Iq", "Id_fil", "Iq_fil", "I_ccd", "I_ccq",
    "err_d", "err_q", "e_ac_d", "e_ac_q",
    "P_expr", "Q_expr", "vd", "vq", "omega_pll_feedback"
]

SOLVE_KEYS = [
    "V_inv", "theta_inv", "theta_pll", "delta_pll",
    "vq_fil", "I_pll", "P_fil", "Q_fil", "I_pqd", "I_pqq",
    "Id_ref", "Iq_ref", "Id", "Iq", "Id_fil", "Iq_fil",
    "I_ccd", "I_ccq", "err_d", "err_q", "e_ac_d", "e_ac_q"
]


# def safe_solve(A, b):
#     try:
#         return np.linalg.solve(A, b)
#     except np.linalg.LinAlgError:
#         return np.linalg.lstsq(A, b, rcond=None)[0]

def safe_solve(A, b, use_sparse=False):
    if use_sparse:
        A_sparse = sp.csc_matrix(A)
        try:
            return spla.spsolve(A_sparse, b)
        except Exception:
            return np.linalg.lstsq(A, b, rcond=None)[0]
    else:
        try:
            return np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            return np.linalg.lstsq(A, b, rcond=None)[0]


def eval_expr(expr):
    val = value(expr, exception=False)
    if val is None:
        raise ValueError(f"Pyomo expression evaluated to None: {expr}")
    return float(val)


def abc_from_vtheta(V_inv, theta_inv, t, omega):
    va = V_inv * math.cos(omega * t + theta_inv)
    vb = V_inv * math.cos(omega * t + theta_inv - 2 * math.pi / 3)
    vc = V_inv * math.cos(omega * t + theta_inv + 2 * math.pi / 3)
    return va, vb, vc


def dq_from_abc(va, vb, vc, theta_pll, t, omega):
    omega_pll = omega * t + theta_pll
    vd = (2 / 3) * (
        va * math.cos(omega_pll)
        + vb * math.cos(omega_pll - 2 * math.pi / 3)
        + vc * math.cos(omega_pll + 2 * math.pi / 3)
    )
    vq = -(2 / 3) * (
        va * math.sin(omega_pll)
        + vb * math.sin(omega_pll - 2 * math.pi / 3)
        + vc * math.sin(omega_pll + 2 * math.pi / 3)
    )
    return vd, vq


def create_history(nfe):
    return {k: np.zeros(nfe) for k in STATE_KEYS}


def get_prev_state(hist, n):
    return {k: hist[k][n - 1] for k in STATE_KEYS}


def set_prev_params(model, prev, t_now, theta_grid_now, Vmag_now):
    model.t_now.set_value(float(t_now))
    model.theta_grid_now.set_value(float(theta_grid_now))
    model.Vmag_now.set_value(float(Vmag_now))
    for k, v in prev.items():
        if k in model.prev:
            model.prev[k].set_value(float(v))


def make_initial_guess(prev):
    return np.array([prev[k] for k in SOLVE_KEYS], dtype=float)


def set_var_vector(model, vars_order, x):
    for i, v in enumerate(vars_order):
        v.set_value(float(x[i]))


def eval_residual_and_jacobian(model, residuals, vars_order):
    F = np.array([eval_expr(r) for r in residuals], dtype=float)
    J = np.zeros((len(residuals), len(vars_order)))
    for i, r in enumerate(residuals):
        J[i, :] = np.array(
            differentiate(r, wrt_list=vars_order, mode=Modes.reverse_numeric),
            dtype=float
        )
    return F, J