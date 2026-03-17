from helpers import create_history, get_prev_state
from initialization import initialize_history
from model_builder import build_reusable_step_model
from solver import solve_one_step


def run_simulation(p, method="schur"):
    hist = create_history(p["nfe"])
    initialize_history(hist, p)

    model, vars_order, residual_exprs = build_reusable_step_model(p)

    for n in range(1, p["nfe"]):
        prev = get_prev_state(hist, n)

        sol = solve_one_step(
            model=model,
            vars_order=vars_order,
            residual_exprs=residual_exprs,
            p=p,
            n=n,
            prev=prev,
            method=method,
        )

        for k in hist:
            hist[k][n] = sol[k]

        if n % 500 == 0:
            print(f"[{method}] Completed step {n}/{p['nfe']}")

    return hist