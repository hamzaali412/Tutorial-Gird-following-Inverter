from config import load_parameters
from simulation import run_simulation
from plotting import plot_results
import time

settings = {"T_end": 0.4, "h": 0.0001, "max_iter": 12, "tol": 1e-9, "use_sparse": True,}

p = load_parameters("params.json", settings)
# Choose method: "schur" or "full"
method = "full"
# method = "schur"
hist = run_simulation(p, method=method)
plot_results(hist, p)

############# Comparsion ####################
# Full Newton
# t0 = time.perf_counter()
# hist_full = run_simulation(p, method="full")
# t1 = time.perf_counter()
# full_time = t1 - t0

# # Schur Newton
# t2 = time.perf_counter()
# hist_schur = run_simulation(p, method="schur")
# t3 = time.perf_counter()
# schur_time = t3 - t2

# print(f"Full Newton time : {full_time:.6f} s")
# print(f"Schur time       : {schur_time:.6f} s")
# print(f"Speedup (Full/Schur) = {full_time / schur_time:.3f}")

