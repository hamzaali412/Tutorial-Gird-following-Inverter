import json
import math
import numpy as np


def load_parameters(json_file="params.json", settings=None):
    with open(json_file, "r") as f:
        cfg = json.load(f)

    p = {}

    sys_cfg = cfg["system"]
    ctrl_cfg = cfg["controller"]
    dst_cfg = cfg["disturbance"]
    ref_cfg = cfg["references"]
    ic_cfg = cfg["initial_conditions"]
    slv_cfg = cfg["solver"]
    
    p.update(sys_cfg)
    p.update(ctrl_cfg)
    p.update(dst_cfg)
    p.update(ref_cfg)
    p.update(ic_cfg)
    p.update(slv_cfg)

    p["h"] = settings["h"]
    p["T_end"] = settings["T_end"]
    p["max_iter"] = settings["max_iter"]
    p["tol"] = settings["tol"]
    p["use_sparse"] = settings.get("use_sparse", False)


    p["omega"] = 2 * math.pi * p["f"]
    p["fbase"] = p["f"]
    p["nfe"] = int(p["T_end"] / p["h"])
    p["time_steps"] = np.linspace(0.0, p["T_end"], p["nfe"], endpoint=False)

    p["V_base"] = (math.sqrt(2) / math.sqrt(3)) * p["VN"]
    p["Zbase"] = p["V_base"] / (math.sqrt(2) * p["Sn"] / (math.sqrt(3) * p["VN"]))
    p["omega_b"] = 2 * math.pi * p["f"]

    p["Rf_pu"] = p["Rf"] / p["Zbase"]
    p["Lf_pu"] = p["Lf"] / (p["Zbase"] / p["omega_b"])

    den = p["R_line_pu"] ** 2 + p["X_line_pu"] ** 2
    p["G_pu"] = p["R_line_pu"] / den
    p["B_pu"] = -p["X_line_pu"] / den

    p["kp_pll"] = 1 / (2 * math.pi * p["fbase"] * p["Tf_pll"])
    p["ki_pll"] = p["kp_pll"] / (4 * p["Tf_pll"])

    p["ki_pq"] = p["Kp_pq"] / (10 * p["Tf_pq"])

    p["kp_ig"] = p["Lf_pu"] / (2 * p["Tf_ig"])
    p["ki_ig"] = p["Rf_pu"] / (2 * p["Tf_ig"])

    p["theta_grid_profile"] = np.where(
        p["time_steps"] >= p["phase_jump_time"],
        p["phase_jump_angle"],
        p["theta_grid_init"],
    )

    # Enable to check the volt-var
    p["Vmag_pu_profile"] = np.where(
        p["time_steps"] >= p["phase_jump_time"],
        p["Vmag_dist"],
        p["Vmag_pu"],
    )

    return p