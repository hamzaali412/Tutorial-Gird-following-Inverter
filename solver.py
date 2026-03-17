import numpy as np
from pyomo.environ import value

from helpers import (
    safe_solve,
    eval_expr,
    eval_residual_and_jacobian,
    set_prev_params,
    set_var_vector,
    make_initial_guess,
)


def extract_solution(model):
    return {
        "V_inv": value(model.V_inv),
        "theta_inv": value(model.theta_inv),
        "theta_pll": value(model.theta_pll),
        "vq_fil": value(model.vq_fil),
        "I_pll": value(model.I_pll),
        "delta_pll": value(model.delta_pll),
        "P_fil": value(model.P_fil),
        "Q_fil": value(model.Q_fil),
        "I_pqd": value(model.I_pqd),
        "I_pqq": value(model.I_pqq),
        "Id_ref": value(model.Id_ref),
        "Iq_ref": value(model.Iq_ref),
        "Id": value(model.Id),
        "Iq": value(model.Iq),
        "Id_fil": value(model.Id_fil),
        "Iq_fil": value(model.Iq_fil),
        "I_ccd": value(model.I_ccd),
        "I_ccq": value(model.I_ccq),
        "err_d": value(model.err_d),
        "err_q": value(model.err_q),
        "e_ac_d": value(model.e_ac_d),
        "e_ac_q": value(model.e_ac_q),
        "vd": value(model.vd),
        "vq": value(model.vq),
        "P_expr": value(model.P_expr),
        "Q_expr": value(model.Q_expr),
        "omega_pll_feedback": value(model.omega_pll_feedback),
    }


def line_search_update(model, vars_order, residual_exprs, x, dx, line_search_max_iter):
    alpha = 1.0
    old_F = np.array([eval_expr(r) for r in residual_exprs], dtype=float)
    old_norm = np.linalg.norm(old_F, ord=2)
    accepted = False

    for _ in range(line_search_max_iter):
        x_trial = x + alpha * dx
        set_var_vector(model, vars_order, x_trial)

        F_trial = np.array([eval_expr(r) for r in residual_exprs], dtype=float)
        trial_norm = np.linalg.norm(F_trial, ord=2)

        if trial_norm < old_norm:
            return x_trial, True

        alpha *= 0.5

    x_trial = x + 0.1 * dx
    set_var_vector(model, vars_order, x_trial)
    return x_trial, False


def full_newton_step(model, vars_order, residual_exprs, p, n, prev):
    t_now = float(p["time_steps"][n])
    theta_grid_now = float(p["theta_grid_profile"][n])
    Vmag_now = float(p["Vmag_pu_profile"][n])

    set_prev_params(model, prev, t_now, theta_grid_now, Vmag_now)

    x = make_initial_guess(prev)
    set_var_vector(model, vars_order, x)
    for it in range(p["max_iter"]):
        F, J = eval_residual_and_jacobian(model, residual_exprs, vars_order)
        full_norm = np.linalg.norm(F, ord=2)

        if p["verbose"]:
            print(f"[Full NR] Step {n:5d}, iter {it:2d}, ||F||={full_norm:.3e}")

        if full_norm < p["tol"]:
            break

        # dx = safe_solve(J, -F)
        dx = safe_solve(J, -F, use_sparse=p["use_sparse"])
        x, _ = line_search_update(
            model=model,
            vars_order=vars_order,
            residual_exprs=residual_exprs,
            x=x,
            dx=dx,
            line_search_max_iter=p["line_search_max_iter"],
        )

    return extract_solution(model)


def schur_newton_step(model, vars_order, residual_exprs, p, n, prev):
    t_now = float(p["time_steps"][n])
    theta_grid_now = float(p["theta_grid_profile"][n])
    Vmag_now = float(p["Vmag_pu_profile"][n])

    set_prev_params(model, prev, t_now, theta_grid_now, Vmag_now)

    x = make_initial_guess(prev)
    set_var_vector(model, vars_order, x)

    nN = 4
    idx_N = list(range(nN))
    idx_L = list(range(nN, len(vars_order)))

    for it in range(p["max_iter"]):
        F, J = eval_residual_and_jacobian(model, residual_exprs, vars_order)

        F_N = F[:nN]
        F_L = F[nN:]

        D = J[np.ix_(idx_N, idx_N)]
        B = J[np.ix_(idx_N, idx_L)]
        C = J[np.ix_(idx_L, idx_N)]
        A = J[np.ix_(idx_L, idx_L)]

        # A_inv_F_L = safe_solve(A, F_L)
        A_inv_F_L = safe_solve(A, F_L, use_sparse=p["use_sparse"])
        # A_inv_C = np.column_stack([safe_solve(A, C[:, k]) for k in range(C.shape[1])])
        # A_inv_C = safe_solve(A, C)
        A_inv_C = safe_solve(A, C, use_sparse=p["use_sparse"])

        S = D - B @ A_inv_C
        F_red = F_N - B @ A_inv_F_L

        full_norm = np.linalg.norm(F, ord=2)
        red_norm = np.linalg.norm(F_red, ord=2)

        if p["verbose"]:
            print(f"[Schur]   Step {n:5d}, iter {it:2d}, ||F_red||={red_norm:.3e}, ||F||={full_norm:.3e}")

        if full_norm < p["tol"]:
            break

        # dx_N = safe_solve(S, -F_red)
        dx_N = safe_solve(S, -F_red, use_sparse=p["use_sparse"])
        # dx_L = -safe_solve(A, F_L + C @ dx_N)
        dx_L = -safe_solve(A, F_L + C @ dx_N, use_sparse=p["use_sparse"])

        dx = np.zeros_like(x)
        dx[idx_N] = dx_N
        dx[idx_L] = dx_L

        x, _ = line_search_update(
            model=model,
            vars_order=vars_order,
            residual_exprs=residual_exprs,
            x=x,
            dx=dx,
            line_search_max_iter=p["line_search_max_iter"],
        )

    return extract_solution(model)


def solve_one_step(model, vars_order, residual_exprs, p, n, prev, method="schur"):
    if method.lower() == "schur":
        return schur_newton_step(model, vars_order, residual_exprs, p, n, prev)
    elif method.lower() in {"full", "full_nr", "newton", "nr"}:
        return full_newton_step(model, vars_order, residual_exprs, p, n, prev)
    else:
        raise ValueError(f"Unknown solve method: {method}")

