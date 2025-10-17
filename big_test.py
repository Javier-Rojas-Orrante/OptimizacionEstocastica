from typing import Dict, Tuple
import numpy as np
import math
import csv
from collections import Counter

try:
    import pulp as pl
except Exception as e:
    raise RuntimeError("Please install PuLP first: pip install pulp") from e

from TresBolillosSpace import TresBolillosSpace  # <-- adjust filename if needed


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


def build_and_solve(space: TresBolillosSpace,
                    targets: Dict[str, int],
                    p_not_surv: Dict[str, float],
                    tol: float = 0.05,
                    alpha_comp: float = 1.0,
                    beta_surv: float = 1.0,
                    comp_seed: int = 0,
                    comp_sparsity: float = 0.0,
                    comp_diag_bonus: float = 1.0,
                    comp_offdiag_scale: float = 1.0,
                    solver: pl.LpSolver_CMD | None = None
                    ) -> Tuple[pl.LpProblem, dict, dict, dict]:
    """
    Returns:
      model: PuLP model
      x_val: {(u,i): 0/1} assignments
      y_val: (empty dict by default; can extract if needed)
      summary: status, objective, counts, bands, sizes
    """
    # Ensure we have a pre-planted state
    if space.y_init is None or space.species_init is None:
        space.sample_initial(seed=comp_seed)

    N, edges, species = space.to_pulp()
    n_nodes = len(N)
    n_species = len(species)

    assert sum(targets.values()) == n_nodes, "Targets must sum to number of nodes."
    assert set(targets.keys()) == set(species), "Targets keys must match species list."
    assert set(p_not_surv.keys()) == set(species), "Survival table keys must match species."

    # Quota bands
    bands = space.quota_bands_total(targets, tol=tol)

    # Survival reward per species
    surv = {name: 1.0 - float(p_not_surv[name]) for name in species}

    # Competition matrix
    C = make_competition_matrix(
        species,
        seed=comp_seed,
        sparsity=comp_sparsity,
        diag_bonus=comp_diag_bonus,
        offdiag_scale=comp_offdiag_scale
    )

    model = pl.LpProblem("Reforestation_Layout", pl.LpMinimize)

    # Decision variables
    x = pl.LpVariable.dicts("x",
                            ((u, i) for u in N for i in range(n_species)),
                            lowBound=0, upBound=1, cat=pl.LpBinary)

    # Only create y for pairs with C[i,j] > 0
    E = [(u, v) for (u, v) in edges]  # undirected with u < v assumed by class
    pair_index = [((u, v), i, j)
                  for (u, v) in E
                  for i in range(n_species)
                  for j in range(n_species)
                  if C[i, j] > 0.0]

    y = pl.LpVariable.dicts("y", pair_index, lowBound=0, upBound=1, cat=pl.LpBinary)

    # Constraints

    # (1) One species per node
    for u in N:
        model += pl.lpSum(x[(u, i)] for i in range(n_species)) == 1, f"one_species_node_{u}"

    # (2) Lock pre-planted nodes exactly
    space.add_preplanted_constraints_pulp(model, x, N, space.species_init, n_species)

    # (3) Species totals within bands
    name_to_idx = {name: i for i, name in enumerate(species)}
    for name, bu in bands.items():
        i = name_to_idx[name]
        total_i = pl.lpSum(x[(u, i)] for u in N)
        model += total_i >= bu["L"], f"quota_L_{name}"
        model += total_i <= bu["U"], f"quota_U_{name}"

    # (4) Linearization for neighbor competition y_uvij = x_ui * x_vj
    for ((u, v), i, j) in pair_index:
        model += y[((u, v), i, j)] <= x[(u, i)]
        model += y[((u, v), i, j)] <= x[(v, j)]
        model += y[((u, v), i, j)] >= x[(u, i)] + x[(v, j)] - 1

    # Objective = alpha * competition - beta * survival
    comp_term = pl.lpSum(C[i, j] * y[((u, v), i, j)] for ((u, v), i, j) in pair_index)
    surv_term = pl.lpSum(surv[species[i]] * x[(u, i)] for u in N for i in range(n_species))
    model += alpha_comp * comp_term - beta_surv * surv_term

    # Solve
    if solver is None:
        solver = pl.PULP_CBC_CMD(
    msg=1,                 # print to console
    options=[
        'log', '5',        # CBC verbosity (0..5; 4–5 is very chatty)
        'threads', '8',    # if you want
        # 'ratio','0.01',  # mip gap target (optional)
    ],
    # keepFiles=True,      # keep temp files (see #3)
)
    status = model.solve(solver)

    # Results
    x_val = {(u, i): int(pl.value(x[(u, i)]) > 0.5) for u in N for i in range(n_species)}
    y_val = {}  # huge; skip unless needed

    # Species counts
    counts = {name: 0 for name in species}
    for u in N:
        for i, name in enumerate(species):
            if x_val[(u, i)] == 1:
                counts[name] += 1
                break

    summary = {
        "status": pl.LpStatus[status],
        "objective": float(pl.value(model.objective)),
        "species_counts": counts,
        "bands": bands,
        "n_nodes": n_nodes,
        "n_edges": len(E),
        "avg_degree": round(2 * len(E) / n_nodes, 3),
        "n_x_vars": len(x),
        "n_y_vars": len(y),
        "n_constraints": len(model.constraints),
    }
    return model, x_val, y_val, summary


if __name__ == "__main__":
    # Build your space from the class that lives in another file
    # Example with rectangular grid 14x47 = 658 nodes
    space = TresBolillosSpace.from_rect(rows=14, cols=47, seed=42)

    # Ensure pre-planting is sampled (or if you already sampled elsewhere, skip)
    space.sample_initial(seed=43)

    rem, checks = space.remaining_bands(TARGETS, tol=TOL)
    print("free_slots:", checks["free_slots"])
    print("sum_L_rem:", checks["sum_L_rem"], " sum_U_rem:", checks["sum_U_rem"])
    print("min_feasible_given_free?:", checks["min_feasible_given_free?"])
    print("max_feasible_given_free?:", checks["max_feasible_given_free?"])
    print("species_over_cap:", checks["species_over_cap"])


    # Solve
    model, x_val, y_val, summary = build_and_solve(
        space=space,
        targets=TARGETS,
        p_not_surv=P_NOT_SURV,
        tol=TOL,
        alpha_comp=ALPHA_COMP,
        beta_surv=BETA_SURV,
        comp_seed=SEED,
        comp_sparsity=COMP_SPARSITY,
        comp_diag_bonus=COMP_DIAG_BONUS,
        comp_offdiag_scale=COMP_OFFDIAG_SCALE,
        solver=None  # default CBC
    )

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
