from typing import Dict, Tuple
import numpy as np
import math
import csv
from collections import Counter
import gurobipy as gp
import os
from gurobipy import GRB
#Halamadrid_10

try:
    import pulp as pl
except Exception as e:
    raise RuntimeError("Please install PuLP first: pip install pulp") from e

from TresBolillosSpace import TresBolillosSpace  # <-- adjust filename if needed

# ---------- helper: proportional downscale of TARGETS to any node count ----------
# ---- helpers for smaller runs with the SAME rules ----
def downscale_targets_exact(original: Dict[str, int], new_total: int) -> Dict[str, int]:
    """Apportion to sum exactly new_total (Hamilton, largest remainder)."""
    items = list(original.items())
    tot_old = sum(v for _, v in items)
    shares = [(k, v * new_total / tot_old) for k, v in items]
    base = {k: int(math.floor(s)) for k, s in shares}
    remainder = new_total - sum(base.values())
    fracs = sorted(((k, s - base[k]) for k, s in shares), key=lambda x: x[1], reverse=True)
    for i in range(remainder):
        base[fracs[i][0]] += 1
    return base

def _pre_counts(space: TresBolillosSpace) -> Dict[str, int]:
    """Count pre-planted species in the current space state."""
    cnt = Counter()
    for u in space.N:
        i = int(space.species_init[u])
        if i >= 0 and int(space.y_init[u]) == 1:
            cnt[space.species[i]] += 1
    return cnt

def sample_preplant_until_feasible(space: TresBolillosSpace,
                                   targets_small: Dict[str, int],
                                   tol: float,
                                   seed0: int,
                                   tries: int = 200) -> bool:
    """
    Use the SAME preplant sampler, but reject samples that already violate
    small targets (e.g., preplant counts > target when TOL=0).
    """
    for t in range(tries):
        space.sample_initial(seed=seed0 + t)

        # hard guard: no species can exceed its (tight) band from preplant alone
        pc = _pre_counts(space)
        if any(pc[s] > targets_small[s] for s in targets_small):
            continue

        rem, checks = space.remaining_bands(targets_small, tol=tol)
        if checks.get("min_feasible_given_free?") and checks.get("max_feasible_given_free?"):
            return True
    return False

def run_smaller_with_same_rules(rows: int, cols: int,
                                tol: float,
                                comp_seed: int,
                                space_seed: int,
                                preplant_seed: int,
                                comp_sparsity: float,
                                comp_diag_bonus: float,
                                comp_offdiag_scale: float,
                                alpha_comp: float,
                                beta_surv: float):
    """Build smaller grid and solve with identical constraints & locks."""
    space = TresBolillosSpace.from_rect(rows=rows, cols=cols, seed=space_seed)

    # Keep species list identical to big instance
    space.species = list(TARGETS.keys())

    # Downscale quotas to new node count (TOL preserved)
    total_small = rows * cols
    TARGETS_SMALL = downscale_targets_exact(TARGETS, total_small)
    assert sum(TARGETS_SMALL.values()) == total_small
    P_NOT_SURV_SMALL = dict(P_NOT_SURV)

    # Preplant with identical mechanism; reject samples that break tight bands
    if not sample_preplant_until_feasible(space, TARGETS_SMALL, tol=tol,
                                          seed0=preplant_seed, tries=500):
        raise RuntimeError(
            "No feasible preplant found for the smaller grid with given tol. "
            "Try a few more nodes (e.g., 7×7 or 8×8) or increase tries."
        )

    # Sanity print to PROVE you're on the small instance
    rem, checks = space.remaining_bands(TARGETS_SMALL, tol=tol)
    print(f"SMALL | nodes={space.n} free_slots={checks['free_slots']} "
          f"sum_L_rem={checks['sum_L_rem']} sum_U_rem={checks['sum_U_rem']} "
          f"min_feas?={checks['min_feasible_given_free?']} "
          f"max_feas?={checks['max_feasible_given_free?']}")

    # --------- PATCH A: pasar targets/es y tol del small ----------
    model, x_val, y_val, summary =  build_and_solve_gurobi(
            space=space,
            targets=TARGETS_SMALL,         # <-- small targets
            p_not_surv=P_NOT_SURV_SMALL,   # <-- same survival map
            tol=tol,                       # <-- tol from small test
            alpha_comp=alpha_comp,
            beta_surv=beta_surv,
            comp_seed=comp_seed,
            comp_sparsity=comp_sparsity,
            comp_diag_bonus=comp_diag_bonus,
            comp_offdiag_scale=comp_offdiag_scale
        )
    # --------- END PATCH A ----------

    print("\n=== SMALL RUN (same rules, preplant locked) ===")
    print(f"Grid: {rows}×{cols} | Status: {summary['status']} | Obj: {round(summary['objective'], 6)}")
    print(f"Nodes: {summary['n_nodes']} | Edges: {summary['n_edges']} | Avg degree: {summary['avg_degree']}")
    print(f"Vars: x={summary['n_x_vars']}  y={summary['n_y_vars']} | Constraints: {summary['n_constraints']}")
    print("\nSpecies counts (solution) vs bands:")
    for name in space.species:
        c = summary["species_counts"][name]
        L = summary["bands"][name]["L"]
        U = summary["bands"][name]["U"]
        print(f"  {name:25s}: {c:3d}  (band: [{L:2d}, {U:2d}])")

    return space, model, x_val, y_val, summary


# -------------------------- Inputs --------------------------

TARGETS = {
    "Agave lechuguilla":        42,
    "Agave salmiana":          196,
    "Agave scabra":             42,
    "Agave striata":            42,
    "Opuntia cantabrigiensis":  49,
    "Opuntia engelmannii":      38,
    "Opuntia robusta":          73,
    "Opuntia streptacantha":    64,
    "Prosopis laevigata":       86,
    "Yucca filifera":           26,
}
assert sum(TARGETS.values()) == 658

P_NOT_SURV = {
    "Agave lechuguilla":        0.03675,
    "Agave salmiana":           0.04590,
    "Agave scabra":             0.12735,
    "Agave striata":            0.04310,
    "Opuntia cantabrigiensis":  0.08565,
    "Opuntia engelmannii":      0.09030,
    "Opuntia robusta":          0.04010,
    "Opuntia streptacantha":    0.07375,
    "Prosopis laevigata":       0.09105,
    "Yucca filifera":           0.13685,
}

# Tolerance ±5% on quotas
TOL = 0.05

# Objective trade-off weights
ALPHA_COMP = 1.0   # weight on neighbor competition (penalize)
BETA_SURV  = 1.0   # weight on expected survival (reward)

# Random competition matrix controls
SEED = 1234
COMP_SPARSITY = 0.10       # 0 → dense, 1 → keep only diagonal (same-species)
COMP_DIAG_BONUS = 1.0      # extra penalty for same-species neighbors
COMP_OFFDIAG_SCALE = 1.0   # scale off-diagonal penalties


def make_competition_matrix(species: list,
                            seed: int = 0,
                            sparsity: float = 0.0,
                            diag_bonus: float = 1.0,
                            offdiag_scale: float = 1.0) -> np.ndarray:
    """
    Symmetric nonnegative competition matrix C (bigger => worse).
    Diagonal can be larger to discourage clumping of same species.
    Some off-diagonals can be pruned by 'sparsity' to reduce model size.
    """
    rng = np.random.default_rng(seed)
    m = len(species)
    A = rng.random((m, m))
    C = (A + A.T) / 2.0
    np.fill_diagonal(C, 0.0)

    if sparsity > 0.0:
        mask = rng.random((m, m)) < sparsity
        mask = np.triu(mask, 1)
        mask = mask + mask.T
        C = np.where(mask, 0.0, C)

    C *= offdiag_scale
    np.fill_diagonal(C, C.diagonal() + diag_bonus)
    return C

def build_and_solve_gurobi(space,
                           targets,
                           p_not_surv,
                           tol=0.05,
                           alpha_comp=1.0,
                           beta_surv=1.0,
                           comp_seed=0,
                           comp_sparsity=0.0,
                           comp_diag_bonus=1.0,
                           comp_offdiag_scale=1.0):
    """
    Pure Gurobi version of build_and_solve().
    Handles both large and small grids, and returns best incumbent solution.
    """

    # --- Connect to Gurobi with WLS credentials ---
    env = gp.Env(params={
        "WLSAccessID": "18f52dbc-253b-4d74-adc3-6e2d1c4b89d2",
        "WLSSecret":   "3cab943f-8d9c-4e49-92a4-18e8b6d90024",
        "LicenseID":   2727523,
        "LogToConsole": 1
    })

    N, edges, species = space.to_pulp()
    n_nodes = len(N)
    n_species = len(species)
    total_targets = sum(targets.values())
    if total_targets != n_nodes:
        # Adjust the first species to fill the remaining nodes
        first_species = list(targets.keys())[0]
        targets[first_species] += n_nodes - total_targets
        print(f"Adjusted targets: sum now = {sum(targets.values())}, targets = {targets}")

    assert sum(targets.values()) == n_nodes, "Targets must sum to number of nodes"

    # --- Parameters ---
    bands = space.quota_bands_total(targets, tol=tol)
    surv = {name: 1.0 - float(p_not_surv[name]) for name in species}

    # --- Fixed competition matrix ---
    C = np.array([
        [1.0, 0.8, 0.8, 0.8, 0.3, 0.3, 0.3, 0.3, 0.5, 0.2],
        [0.8, 1.0, 0.8, 0.8, 0.3, 0.3, 0.3, 0.3, 0.5, 0.2],
        [0.8, 0.8, 1.0, 0.8, 0.3, 0.3, 0.3, 0.3, 0.5, 0.2],
        [0.8, 0.8, 0.8, 1.0, 0.3, 0.3, 0.3, 0.3, 0.5, 0.2],
        [0.3, 0.3, 0.3, 0.3, 1.0, 0.8, 0.8, 0.8, 0.5, 0.2],
        [0.3, 0.3, 0.3, 0.3, 0.8, 1.0, 0.8, 0.8, 0.5, 0.2],
        [0.3, 0.3, 0.3, 0.3, 0.8, 0.8, 1.0, 0.8, 0.5, 0.2],
        [0.3, 0.3, 0.3, 0.3, 0.8, 0.8, 0.8, 1.0, 0.5, 0.2],
        [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.3],
        [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 1.0],
    ])
    assert C.shape == (len(species), len(species)), "Competition matrix size mismatch!"

    name_to_idx = {name: i for i, name in enumerate(species)}

    # --- Model setup ---
    model = gp.Model("Reforestation_Layout", env=env)

    # Time & performance parameters
    model.Params.TimeLimit = 600 if n_nodes <= 64 else 1800
    model.Params.Threads = min(8, os.cpu_count())
    model.Params.MIPGap = 0.02
    model.Params.MIPFocus = 1
    model.Params.Heuristics = 0.2
    model.Params.Presolve = 2
    model.Params.Cuts = 2
    model.Params.Symmetry = 2

    # --- Decision variables ---
    x = model.addVars(n_nodes, n_species, vtype=GRB.BINARY, name="x")

    E = [(u, v) for (u, v) in edges]
    pair_index = [((u, v), i, j)
                  for (u, v) in E
                  for i in range(n_species)
                  for j in range(n_species)
                  if C[i, j] > 0.0]
    y = model.addVars(pair_index, vtype=GRB.BINARY, name="y")

    # --- Constraints ---
    # (1) One species per node
    for u in N:
        model.addConstr(gp.quicksum(x[u, i] for i in range(n_species)) == 1)

    # (2) Lock preplanted
    for u in N:
        pre_i = int(space.species_init[u])
        if pre_i >= 0 and int(space.y_init[u]) == 1:
            for i in range(n_species):
                model.addConstr(x[u, i] == (1 if i == pre_i else 0))

    # (3) Quotas within bands
    for name, bu in bands.items():
        i = name_to_idx[name]
        total_i = gp.quicksum(x[u, i] for u in N)
        model.addConstr(total_i >= bu["L"])
        model.addConstr(total_i <= bu["U"])

    # (4) Linearization
    for ((u, v), i, j) in pair_index:
        model.addConstr(y[((u, v), i, j)] <= x[u, i])
        model.addConstr(y[((u, v), i, j)] <= x[v, j])
        model.addConstr(y[((u, v), i, j)] >= x[u, i] + x[v, j] - 1)

    # --- Objective ---
    surv = {name: 1.0 - float(p_not_surv[name]) for name in species}
    surv_factor = {(i, j): 1.0 / (1.0 - surv[species[i]] * surv[species[j]]) for i in range(n_species) for j in range(n_species)}
    adj_term = gp.quicksum(0.5 * C[i, j] * surv_factor[(i, j)] * y[((u, v), i, j)]
                           for ((u, v), i, j) in pair_index)
    model.setObjective(alpha_comp * adj_term, GRB.MINIMIZE)

    # --- Solve ---
    print(f"Solving model with {n_nodes} nodes and {len(x)} x-vars, {len(y)} y-vars...")
    model.optimize()

    # --- Extract results safely ---
    if model.SolCount == 0:
        print("⚠️ No feasible solution found.")
        return model, {}, {}, {"status": model.Status, "objective": None}

    x_val = {(u, i): int(x[u, i].X > 0.5) if x[u, i].X is not None else 0 for u in N for i in range(n_species)}

    counts = {name: 0 for name in species}
    for u in N:
        for i, name in enumerate(species):
            if x_val[(u, i)] == 1:
                counts[name] += 1
                break

    summary = {
        "status": model.Status,
        "objective": model.ObjVal if model.SolCount > 0 else None,
        "species_counts": counts,
        "bands": bands,
        "n_nodes": n_nodes,
        "n_edges": len(E),
        "avg_degree": round(2 * len(E) / n_nodes, 3),
        "n_x_vars": len(x),
        "n_y_vars": len(y),
        "n_constraints": model.NumConstrs,
    }

    print(f"\n✅ Model done: Obj={summary['objective']:.3f}, Gap={model.MIPGap*100:.2f}% | Status={model.Status}")
    return model, x_val, {}, summary

if __name__ == "__main__":
    USE_SMALL = True  # ← set True to run the smaller space, False for your 14×47

    if USE_SMALL:
        # e.g., 10×10 = 100 nodes with the SAME rules (bands via tol, preplant locked, same objective)
        space, model, x_val, y_val, summary = run_smaller_with_same_rules(
            rows=10, cols=10,
            tol=TOL,                      # usa el tol que definiste arriba
            comp_seed=SEED,
            space_seed=42,
            preplant_seed=43,
            comp_sparsity=COMP_SPARSITY,
            comp_diag_bonus=COMP_DIAG_BONUS,
            comp_offdiag_scale=COMP_OFFDIAG_SCALE,
            alpha_comp=ALPHA_COMP,
            beta_surv=BETA_SURV
        )

    else:
        # ----- original big run (unchanged) -----
        space = TresBolillosSpace.from_rect(rows=14, cols=47, seed=42)
        space.sample_initial(seed=43)

        rem, checks = space.remaining_bands(TARGETS, tol=TOL)
        print("free_slots:", checks["free_slots"])
        print("sum_L_rem:", checks["sum_L_rem"], " sum_U_rem:", checks["sum_U_rem"])
        print("min_feasible_given_free?:", checks["min_feasible_given_free?"])
        print("max_feasible_given_free?:", checks["max_feasible_given_free?"])
        print("species_over_cap:", checks["species_over_cap"])

        model, x_val, y_val, summary = build_and_solve_gurobi(
            space=space,
            targets=TARGETS,
            p_not_surv=P_NOT_SURV,
            tol=TOL,
            alpha_comp=ALPHA_COMP,
            beta_surv=BETA_SURV,
            comp_seed=SEED,
            comp_sparsity=COMP_SPARSITY,
            comp_diag_bonus=COMP_DIAG_BONUS,
            comp_offdiag_scale=COMP_OFFDIAG_SCALE
        )

    # ---------- common reporting (works for BOTH small and big) ----------
    print("Status:", summary["status"])
    print("Objective:", round(summary["objective"], 6))
    print("Nodes:", summary["n_nodes"], "| Edges:", summary["n_edges"], "| Avg degree:", summary["avg_degree"])
    print("Vars: x =", summary["n_x_vars"], " y =", summary["n_y_vars"], "| Constraints:", summary["n_constraints"])
    print("\nSpecies counts (solution) vs bands:")
    for name in space.species:
        c = summary["species_counts"][name]
        L = summary["bands"][name]["L"]
        U = summary["bands"][name]["U"]
        print(f"  {name:25s}: {c:3d}  (band: [{L:3d}, {U:3d}])")

    # ----- the rest of your saving/plotting block stays exactly as-is -----
    # (chosen_idx build, CSVs, NPZ, optional plot)
    # ...

    # ----- Build chosen species index per node (0..m-1) -----
    chosen_idx = np.full(space.n, -1, dtype=int)
    for u in space.N:
        for i in range(len(space.species)):
            if x_val[(u, i)] == 1:
                chosen_idx[u] = i
                break

    # ---------------- SAVE SOLUTION TO DISK -----------------
    # Positions for plotting / GIS-ish export
    if space.mode == "rect":
        pos = space._positions_rect(spacing=1.0)   # (x,y) for visualization
        rows, cols = space.rows, space.cols
        def rc(u): return u // cols, u % cols
        has_rc = True
    else:
        pos = space._positions_disk(spacing=1.0)
        has_rc = False

    # Node-level CSV
    csv_nodes = "assignment_nodes.csv"
    with open(csv_nodes, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        header = ["node", "species_idx", "species_name", "preplanted", "pre_species_name", "x", "y"]
        if has_rc:
            header.extend(["row", "col"])
        w.writerow(header)

        for u in space.N:
            s_idx = int(chosen_idx[u])
            s_name = space.species[s_idx]
            pre_i = int(space.species_init[u])
            preplanted = 1 if pre_i >= 0 else 0
            pre_name = space.species[pre_i] if pre_i >= 0 else ""
            row = [u, s_idx, s_name, preplanted, pre_name, float(pos[u,0]), float(pos[u,1])]
            if has_rc:
                r, c = rc(u)
                row.extend([r, c])
            w.writerow(row)
    print(f"\nSaved per-node assignment to: {csv_nodes}")

    # Species totals CSV
    counts = Counter(space.species[i] for i in chosen_idx)
    csv_counts = "assignment_species_counts.csv"
    with open(csv_counts, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["species_name", "count", "band_L", "band_U"])
        for name in space.species:
            L = summary["bands"][name]["L"]
            U = summary["bands"][name]["U"]
            w.writerow([name, counts[name], L, U])
    print(f"Saved species totals to: {csv_counts}")

    # Numpy bundle for quick reload
    np.savez("assignment_solution.npz",
             chosen_idx=chosen_idx,
             chosen_name=np.array([space.species[i] for i in chosen_idx], dtype=object),
             pos=pos,
             species=np.array(space.species, dtype=object))
    print("Saved numpy bundle to: assignment_solution.npz")

    # Optional quick plot using your class' plot() (requires matplotlib)
    try:
        import matplotlib.pyplot as plt
        space.y_init = np.ones(space.n, dtype=int)  # mark all as occupied for plotting
        space.species_init = chosen_idx.copy()      # show chosen assignment
        space.plot(figsize=(12, 5), point_size_occ=24, point_size_empty=5)
    except Exception:
        pass


#--- Solvers ---
# # Gurobi (academic free)
# pl.LpSolverDefault = pl.GUROBI_CMD(msg=1, timeLimit=1800, threads=8)

# # or CPLEX
# pl.LpSolverDefault = pl.CPLEX_PY(msg=1, timelimit=1800)

# # or SCIP
# pl.LpSolverDefault = pl.SCIP_CMD(msg=True, limits_time=1800, threads=8)
