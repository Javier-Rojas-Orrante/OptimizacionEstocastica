# -*- coding: utf-8 -*-
# Reduced-space MILP with ALL 10 species, targets downscaled proportionally,
# SAME restrictions (incl. pre-plant locks), FIXED competition matrix.

from typing import Dict, List, Tuple
import math
import numpy as np
import pulp as pl
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import get_cmap

# ======================= SPACE / GRAPH =======================
try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None

class TresBolillosSpace:
    DEFAULT_SPECIES = [
        "Agave lechuguilla",       # AL
        "Agave salmiana",          # AS
        "Agave scabra",            # ASc
        "Agave striata",           # ASt
        "Opuntia cantabrigiensis", # OC
        "Opuntia engelmannii",     # OE
        "Opuntia robusta",         # OR
        "Opuntia streptacantha",   # OS
        "Prosopis laevigata",      # PL
        "Yucca filifera",          # YF
    ]

    DEFAULT_SPECIES_PROBS = {
        "Agave lechuguilla":        0.06326106326,
        "Agave salmiana":           0.295133438,
        "Agave scabra":             0.06203063346,
        "Agave striata":            0.05956977386,
        "Opuntia cantabrigiensis":  0.07772922059,
        "Opuntia engelmannii":      0.06190334762,
        "Opuntia robusta":          0.1165089737,
        "Opuntia streptacantha":    0.09711909712,
        "Prosopis laevigata":       0.1280495566,
        "Yucca filifera":           0.03869489584,
    }

    def __init__(self, mode="rect", rows=14, cols=47, n_sites=658,
                 p_occupied=0.1963257453, species=None, species_probs=None, seed=None):
        self.mode = mode
        self.rows = rows
        self.cols = cols
        self.n_sites = n_sites
        self.p_occupied = float(p_occupied)
        self.species = list(species) if species is not None else list(self.DEFAULT_SPECIES)
        self.species_probs = dict(species_probs) if species_probs is not None else dict(self.DEFAULT_SPECIES_PROBS)
        self.seed = seed

        self.N, self.edges = self._build_graph()
        self.n = len(self.N)

        self.y_init = None        # np.ndarray (0/1)
        self.species_init = None  # np.ndarray (-1 or species index)
        self.counts = None        # dict name->count

    @classmethod
    def from_rect(cls, rows=14, cols=47, **kwargs):
        return cls(mode="rect", rows=rows, cols=cols, n_sites=rows*cols, **kwargs)

    @classmethod
    def from_disk(cls, n_sites=658, **kwargs):
        return cls(mode="disk", n_sites=n_sites, **kwargs)

    def _build_graph(self):
        if self.mode == "rect":
            return self._build_hex_edges_rect(self.rows, self.cols)
        elif self.mode == "disk":
            return self._build_hex_edges_n(self.n_sites)
        raise ValueError("mode must be 'rect' or 'disk'.")

    @staticmethod
    def _build_hex_edges_rect(rows, cols):
        n_sites = rows * cols
        N = list(range(n_sites))
        E = set()
        for r in range(rows):
            for c in range(cols):
                u = r*cols + c
                cand = [(r, c-1), (r, c+1), (r-1, c), (r+1, c)]
                cand += [(r-1, c-1), (r+1, c-1)] if (r % 2 == 0) else [(r-1, c+1), (r+1, c+1)]
                for rr, cc in cand:
                    if 0 <= rr < rows and 0 <= cc < cols:
                        v = rr*cols + cc
                        if u < v:
                            E.add((u, v))
        return N, sorted(E)

    @staticmethod
    def _build_hex_edges_n(n_sites):
        def hex_count(R): return 1 + 3*R*(R+1)
        R = 0
        while hex_count(R) < n_sites:
            R += 1
        coords = []
        for q in range(-R, R+1):
            rmin = max(-R, -q-R)
            rmax = min(R, -q+R)
            for r in range(rmin, rmax+1):
                coords.append((q, r))
        coords = coords[:n_sites]
        coord_to_id = {coords[i]: i for i in range(n_sites)}
        N = list(range(n_sites))
        dirs = [(1,0),(1,-1),(0,-1),(-1,0),(-1,1),(0,1)]
        E = set()
        for u in N:
            q,r = coords[u]
            for dq,dr in dirs:
                v = coord_to_id.get((q+dq, r+dr))
                if v is not None and u < v:
                    E.add((u, v))
        return N, sorted(E)

    def _normalized_probs(self):
        p = np.array([self.species_probs[name] for name in self.species], dtype=float)
        s = p.sum()
        if s <= 0:
            raise ValueError("Sum of species probabilities must be > 0.")
        return p / s

    def sample_initial(self, seed=None):
        rng = np.random.default_rng(self.seed if seed is None else seed)
        y = (rng.random(self.n) < self.p_occupied).astype(int)
        species_init = np.full(self.n, -1, dtype=int)
        k = int(y.sum())
        if k > 0:
            draws = rng.choice(len(self.species), size=k, p=self._normalized_probs())
            species_init[y == 1] = draws
        self.y_init = y
        self.species_init = species_init
        self.counts = {self.species[i]: int((species_init == i).sum()) for i in range(len(self.species))}

    def to_pulp(self):
        return self.N, self.edges, self.species

    @staticmethod
    def add_preplanted_constraints_pulp(model, x, N, species_init, n_species):
        for u in N:
            i_star = int(species_init[u])
            if i_star >= 0:
                model += x[(u, i_star)] == 1
                for j in range(n_species):
                    if j != i_star:
                        model += x[(u, j)] == 0

    # ===== quotas / bands =====
    def pre_counts(self):
        if self.species_init is None:
            return {name: 0 for name in self.species}
        counts = {name: 0 for name in self.species}
        for idx in self.species_init:
            if idx >= 0:
                counts[self.species[idx]] += 1
        return counts

    def quota_bands_total(self, targets, tol=0.05):
        bands = {}
        for name in self.species:
            Ti = int(targets[name])
            Li = math.floor((1.0 - tol) * Ti)
            Ui = math.ceil((1.0 + tol) * Ti)
            bands[name] = {"T": Ti, "L": Li, "U": Ui}
        return bands

    def remaining_bands(self, targets, tol=0.05):
        bands = self.quota_bands_total(targets, tol)
        pre = self.pre_counts()
        free_slots = int(self.n - (self.y_init.sum() if self.y_init is not None else 0))
        rem = {}
        infeasible = []
        for name, bu in bands.items():
            P = pre.get(name, 0)
            Lp = max(0, bu["L"] - P)
            Up = max(0, bu["U"] - P)
            if P > bu["U"]:
                infeasible.append((name, f"pre={P} > U={bu['U']}"))
            rem[name] = {"pre": P, "L_rem": Lp, "U_rem": Up}
        min_needed = sum(v["L_rem"] for v in rem.values())
        max_allow = sum(v["U_rem"] for v in rem.values())
        checks = {
            "free_slots": free_slots,
            "sum_L_rem": min_needed,
            "sum_U_rem": max_allow,
            "min_feasible_given_free?": (min_needed <= free_slots),
            "max_feasible_given_free?": (free_slots <= max_allow),
            "species_over_cap": infeasible,
        }
        return rem, checks

# ======================= DATA (FULL) =======================
SPECIES_TARGETS_FULL = {
    "Agave lechuguilla":        42,
    "Agave salmiana":           196,
    "Agave scabra":             42,
    "Agave striata":            42,
    "Opuntia cantabrigiensis":  49,
    "Opuntia engelmannii":      38,
    "Opuntia robusta":          73,
    "Opuntia streptacantha":    64,
    "Prosopis laevigata":       86,
    "Yucca filifera":           26,
}
PROB_NOT_SURVIVE = {
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

# Fixed 10x10 competition matrix (order = DEFAULT_SPECIES)
# Abbrev order: [AL, AS, ASc, ASt, OC, OE, OR, OS, PL, YF]
C_FIXED = np.array([
    [1.0, 0.8, 0.8, 0.8, 0.3, 0.3, 0.3, 0.3, 0.5, 0.2],  # AL
    [0.8, 1.0, 0.8, 0.8, 0.3, 0.3, 0.3, 0.3, 0.5, 0.2],  # AS
    [0.8, 0.8, 1.0, 0.8, 0.3, 0.3, 0.3, 0.3, 0.5, 0.2],  # ASc
    [0.8, 0.8, 0.8, 1.0, 0.3, 0.3, 0.3, 0.3, 0.5, 0.2],  # ASt
    [0.3, 0.3, 0.3, 0.3, 1.0, 0.8, 0.8, 0.8, 0.5, 0.2],  # OC
    [0.3, 0.3, 0.3, 0.3, 0.8, 1.0, 0.8, 0.8, 0.5, 0.2],  # OE
    [0.3, 0.3, 0.3, 0.3, 0.8, 0.8, 1.0, 0.8, 0.5, 0.2],  # OR
    [0.3, 0.3, 0.3, 0.3, 0.8, 0.8, 0.8, 1.0, 0.5, 0.2],  # OS
    [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.3],  # PL
    [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 1.0],  # YF
], dtype=float)

# ======================= HELPERS =======================
def apportion_targets_proportionally(original: Dict[str, int], new_total: int, species_order: List[str]) -> Dict[str, int]:
    """
    Hamilton (largest remainder) apportionment to keep all species proportions and sum exactly new_total.
    """
    items = [(k, original[k]) for k in species_order]
    tot_old = sum(v for _, v in items)
    shares = [(k, v * new_total / tot_old) for k, v in items]
    base = {k: int(math.floor(s)) for k, s in shares}
    remainder = new_total - sum(base.values())
    # sort by fractional part descending
    fracs = sorted(((k, s - base[k]) for k, s in shares), key=lambda x: x[1], reverse=True)
    for i in range(remainder):
        base[fracs[i][0]] += 1
    return base

def sample_preplant_until_feasible(space: TresBolillosSpace,
                                   targets_small: Dict[str, int],
                                   tol: float,
                                   seed0: int,
                                   tries: int = 200) -> bool:
    """
    Uses the SAME preplant mechanism and keeps hard locks; only varies seed to find a feasible draw for bands (tol).
    """
    for t in range(tries):
        space.sample_initial(seed=seed0 + t)
        _, checks = space.remaining_bands(targets_small, tol=tol)
        if checks.get("min_feasible_given_free?") and checks.get("max_feasible_given_free?") and not checks.get("species_over_cap"):
            return True
    return False


#===================== me da flojera mergear
def fitnessCompetencia_edges(space, plantas_matrix, C):
    """Sum C[su, sv] over every undirected edge in space.edges."""
    # IMPORTANT: your helpers already return matrices flipped with np.flipud,
    # so DO NOT flip again here, or you'll misalign indices.
    mat = plantas_matrix
    competencia = 0.0
    rows, cols = mat.shape

    for (u, v) in space.edges:
        ru, cu = divmod(u, cols)
        rv, cv = divmod(v, cols)
        su = mat[ru, cu]
        sv = mat[rv, cv]
        if su >= 0 and sv >= 0:
            competencia += C[su, sv]   # <- was 'total' by mistake
    return competencia



def fitnessDiversidad(plantas, num_especies):
    """
    Calcula la diversidad de especies en la matriz de plantas.
    La suma de las desviaciones de una distribución uniforme de especies.
    Valor ideal es 0 (todas las especies están igualmente representadas).
    """
    freq_especies = np.zeros(num_especies)
    total_plantas = 0
    for i in range(plantas.shape[0]):
        for j in range(plantas.shape[1]):
            especie_nodo = plantas[i, j]
            if especie_nodo >= 0:  # Solo contar nodos con plantas
                freq_especies[especie_nodo] += 1
                total_plantas += 1

    diversidad = 0
    for i in range(num_especies):
        proporcion_actual = freq_especies[i] / total_plantas if total_plantas > 0 else 0
        proporcion_ideal = 1 / num_especies if num_especies > 0 else 0
        diversidad += abs(proporcion_actual - proporcion_ideal)

    return diversidad



# ======================= MILP =======================
def build_and_solve(
    space: TresBolillosSpace,
    targets: Dict[str, int],
    p_not_survive: Dict[str, float],
    C: np.ndarray,
    tol: float = 0.00,      # bands at 0 as requested
    w_comp: float = 1.0,    # penalize competition
    w_surv: float = 1.0,    # reward survival (implemented as -survival in objective)
    time_limit_sec: int | None = 120,
) -> Tuple[Dict[int, int], Dict[str, int], float, dict]:
    # ensure preplant present
    if space.y_init is None or space.species_init is None:
        space.sample_initial()

    N, E, species = space.to_pulp()
    n_nodes, n_edges = len(N), len(E)
    S = len(species)

    # Map inputs to index order
    surv_reward = np.array([1.0 - float(p_not_survive[s]) for s in species], dtype=float)    

    # Bands & feasibility check
    bands = space.quota_bands_total(targets, tol=tol)
    _, checks = space.remaining_bands(targets, tol=tol)
    if checks["species_over_cap"]:
        over = ", ".join(f"{nm} ({msg})" for nm, msg in checks["species_over_cap"])
        raise ValueError("Infeasible (preplant exceeds cap): " + over)
    if not checks["min_feasible_given_free?"] or not checks["max_feasible_given_free?"]:
        raise ValueError(
            f"Infeasible bands vs free slots: free={checks['free_slots']}, "
            f"min_needed={checks['sum_L_rem']}, max_allow={checks['sum_U_rem']}"
        )

    # Model
    model = pl.LpProblem("Reforestation_ReducedSpace", pl.LpMinimize)

    # Decision variables
    x = pl.LpVariable.dicts("x", ((u, i) for u in N for i in range(S)), 0, 1, cat="Binary")
    y = pl.LpVariable.dicts("y", ((u, v, i, j) for (u, v) in E for i in range(S) for j in range(S)), 0, 1, cat="Binary")

    # (1) one species per node
    for u in N:
        model += pl.lpSum(x[(u, i)] for i in range(S)) == 1, f"one_species_node_{u}"

    # (2) preplanted locks
    TresBolillosSpace.add_preplanted_constraints_pulp(model, x, N, space.species_init, S)

    # (3) quotas (bands at tol)
    for i, sname in enumerate(species):
        L_i, U_i = bands[sname]["L"], bands[sname]["U"]
        tot_i = pl.lpSum(x[(u, i)] for u in N)
        model += tot_i >= L_i, f"quota_L_{sname}"
        model += tot_i <= U_i, f"quota_U_{sname}"

    # (4) linearization y[u,v,i,j] = x[u,i] AND x[v,j]
    for (u, v) in E:
        for i in range(S):
            for j in range(S):
                model += y[(u, v, i, j)] <= x[(u, i)]
                model += y[(u, v, i, j)] <= x[(v, j)]
                model += y[(u, v, i, j)] >= x[(u, i)] + x[(v, j)] - 1

   # (5) objective = w_comp * sum(C[i,j]*y)  -  w_surv * sum(surv_reward[i]*x)
    comp_term = pl.lpSum(C[i, j] * y[(u, v, i, j)] for (u, v) in E for i in range(S) for j in range(S))
    surv_term = pl.lpSum(surv_reward[i] * x[(u, i)] for u in N for i in range(S))
    model += w_comp * comp_term

    # Solve
    solver = pl.PULP_CBC_CMD(msg=1, timeLimit=time_limit_sec) if time_limit_sec else pl.PULP_CBC_CMD(msg=1)
    status = model.solve(solver)
    status_str = pl.LpStatus[status]
    if status_str not in ("Optimal", "Integer Feasible", "Feasible"):
        raise RuntimeError(f"Solver status: {status_str}. No integral solution to extract.")


    # Extract assignment
    assign_idx: Dict[int, int] = {}
    for u in N:
        # pick the argmax in case of numerical ties
        vals = [pl.value(x[(u, i)]) for i in range(S)]
        assign_idx[u] = int(np.argmax(vals))

    counts_by_name: Dict[str, int] = {species[i]: 0 for i in range(S)}
    for u, i in assign_idx.items():
        counts_by_name[species[i]] += 1

    summary = {
        "status": pl.LpStatus[status],
        "objective": float(pl.value(model.objective)),
        "n_nodes": n_nodes,
        "n_edges": n_edges,
        "avg_degree": round(2 * n_edges / max(1, n_nodes), 3),
        "n_x_vars": len(x),
        "n_y_vars": len(y),
        "n_constraints": len(model.constraints),
        "bands": bands,
        "species": species,
    }
    return assign_idx, counts_by_name, float(pl.value(model.objective)), summary

# ---------- Robust plot helper (no class internals needed) ----------
def _hex_positions_rect(rows: int, cols: int, spacing: float = 1.0):
    import math
    dx = spacing * 0.5
    dy = spacing * math.sqrt(3) / 2.0
    pts = []
    for r in range(rows):
        for c in range(cols):
            x = c * spacing + (dx if (r % 2) else 0.0)
            y = r * dy
            pts.append((x, y))
    return np.array(pts, dtype=float)

def _grid_positions(n: int, spacing: float = 1.0):
    # fallback when we don't know rows/cols: near-square grid
    import math
    rows = int(math.floor(math.sqrt(n)))
    cols = int(math.ceil(n / rows))
    pts = []
    for u in range(n):
        r, c = divmod(u, cols)
        pts.append((c * spacing, r * spacing))
    return np.array(pts, dtype=float)


def plot_solution(space,
                  assignment: dict | None = None,
                  title: str = "Layout",
                  spacing: float = 1.0,
                  show_edges: bool = True,
                  point_size_occ: int = 26,
                  point_size_empty: int = 6):
    """
    Hex (offset) plot. If `assignment` is provided (dict node->species_idx), plot that;
    else plot the preplant in space.y_init/space.species_init.
    """
    import math
    import numpy as np
    import matplotlib.pyplot as plt

    n = getattr(space, "n", len(getattr(space, "N", [])))
    rows = getattr(space, "rows", None)
    cols = getattr(space, "cols", None)
    mode = getattr(space, "mode", "rect")
    species_names = list(getattr(space, "species", []))
    m = len(species_names)

    # hex positions if rectangular
    if mode == "rect" and rows and cols and rows * cols == n:
        dx = 0.5 * spacing
        dy = math.sqrt(3) / 2.0 * spacing
        pos = np.zeros((n, 2), dtype=float)
        for r in range(rows):
            for c in range(cols):
                u = r * cols + c
                pos[u, 0] = c * spacing + (dx if (r % 2) else 0.0)
                pos[u, 1] = r * dy
    else:
        # simple grid fallback
        side = int(math.ceil(math.sqrt(n)))
        pos = np.column_stack(np.unravel_index(np.arange(n), (side, side)))[0:n]
        pos = pos[:, ::-1].astype(float) * spacing

    # Decide labels to plot
    if assignment is not None:
        lab = np.array([assignment[u] for u in range(n)], dtype=int)
        occ_mask = np.ones(n, dtype=bool)
        empty_mask = ~occ_mask
    else:
        if space.y_init is None or space.species_init is None:
            print("Nothing to plot: provide `assignment` or sample preplant first.")
            return
        lab = np.array(space.species_init, dtype=int)
        occ_mask = (np.array(space.y_init, dtype=int) == 1)
        empty_mask = ~occ_mask

    fig, ax = plt.subplots(figsize=(10, 5))

    if show_edges and hasattr(space, "edges"):
        for (u, v) in space.edges:
            x1, y1 = pos[u]; x2, y2 = pos[v]
            ax.plot([x1, x2], [y1, y2], lw=0.5, alpha=0.25, color="#8aa6c1")

    if empty_mask.any():
        ax.scatter(pos[empty_mask, 0], pos[empty_mask, 1],
                   s=point_size_empty, alpha=0.35, color="#c9d6e4", edgecolors="none")

    # Discrete colormap: one bin per species
    cmap = get_cmap("tab10", max(10, m))
    boundaries = np.arange(-0.5, m + 0.5, 1.0)
    norm = BoundaryNorm(boundaries, cmap.N)

    sc = ax.scatter(pos[occ_mask, 0], pos[occ_mask, 1],
                    s=point_size_occ, c=lab[occ_mask],
                    cmap=cmap, norm=norm, alpha=0.95, edgecolors="none")

    cbar = fig.colorbar(sc, ax=ax, ticks=np.arange(m), fraction=0.035, pad=0.02)
    if species_names and len(species_names) == m:
        cbar.ax.set_yticklabels(species_names)
    cbar.set_label("Especie", rotation=90)
    cbar.ax.tick_params(labelsize=9)

    ax.set_aspect("equal")
    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()
    plt.show()

def assignment_dict_to_matrix(assignment: dict, rows: int, cols: int) -> np.ndarray:
    """
    Turn MILP assignment dict (u->species_idx) into a rows x cols matrix,
    then flip vertically to match how your GA scorer reads it.
    """
    mat = np.full((rows, cols), -1, dtype=int)
    for u, sp in assignment.items():
        r, c = divmod(u, cols)
        mat[r, c] = int(sp)
    return np.flipud(mat)


def preplant_to_matrix(space) -> np.ndarray:
    """
    Build a rows x cols matrix from the preplant (fixed nodes).
    Empty nodes become -1. Flip vertically to match GA scorer.
    """
    rows, cols = space.rows, space.cols
    mat = np.full((rows, cols), -1, dtype=int)
    for u in range(space.n):
        r, c = divmod(u, cols)
        if space.y_init[u] == 1:
            mat[r, c] = int(space.species_init[u])
    return np.flipud(mat)


def score_layout_like_GA(matrix: np.ndarray, C: np.ndarray, n_species: int, space) -> tuple[float, float]:
    comp_sum = fitnessCompetencia_edges(space, matrix, C)
    n_occupied = (matrix >= 0).sum()
    f1_avg = comp_sum / max(n_occupied, 1)
    f2_div = fitnessDiversidad(matrix, n_species)
    return f1_avg, f2_div



# ======================= MAIN (REDUCED SPACE) =======================
if __name__ == "__main__":
    # --- Choose a smaller grid here ---
    rows, cols = 5, 6           # e.g., 6x6 = 36 nodes (change as you like)
    TOL = 0.05                  # bands at 0 (exact targets)
    SEED_SPACE = 42             # geometry / RNG base
    SEED_PREPLANT = 43          # preplant RNG (we will iterate seeds to hit feasibility)
    P_OCCUPIED = 0.1963257453   # keep same occupied probability to preserve rule

    # Build reduced space with ALL species (same order as matrix)
    space = TresBolillosSpace.from_rect(
        rows=rows, cols=cols,
        species=TresBolillosSpace.DEFAULT_SPECIES,
        species_probs=TresBolillosSpace.DEFAULT_SPECIES_PROBS,
        p_occupied=P_OCCUPIED,
        seed=SEED_SPACE
    )

    # Downscale targets proportionally to rows*cols (ALL species kept)
    total_sites = rows * cols
    TARGETS_SMALL = apportion_targets_proportionally(SPECIES_TARGETS_FULL, total_sites, space.species)
    assert sum(TARGETS_SMALL.values()) == total_sites, "Downscaled targets must sum to grid size."

    # Try preplant seeds until bands (TOL=0) are feasible with locks
    ok = sample_preplant_until_feasible(space, TARGETS_SMALL, tol=TOL, seed0=SEED_PREPLANT, tries=400)
    if not ok:
        raise RuntimeError("Could not find a feasible preplant for the reduced bands with TOL=0. Try a slightly larger grid or increase tries.")

    plot_solution(space, assignment=None, title="Preplant (locked nodes)")

    # Solve (same rules; fixed competition matrix)
    assignment, counts, obj, summary = build_and_solve(
        space=space,
        targets=TARGETS_SMALL,
        p_not_survive=PROB_NOT_SURVIVE,
        C=C_FIXED,
        tol=0.05,
        w_comp=1.0,
        w_surv=1.0,
        time_limit_sec=180
    )

    # Report
    # --- GA-style scoring of the MILP solution ---
    milp_matrix = assignment_dict_to_matrix(assignment, rows, cols)
    f1_milp, f2_milp = score_layout_like_GA(milp_matrix, C_FIXED, len(space.species), space)

    print("\n[Scoring with GA metrics]")
    print(f"[MILP] competencia_avg = {f1_milp:.6f}   diversidad_dev = {f2_milp:.6f}")

    # --- Baseline: preplant ---

    preplant_matrix = preplant_to_matrix(space)
    f1_pre, f2_pre = score_layout_like_GA(preplant_matrix, C_FIXED, len(space.species), space)

    print(f"[BASE] competencia_avg = {f1_pre:.6f}   diversidad_dev = {f2_pre:.6f}")


    print("\n=== REDUCED SPACE RUN (all species, proportional targets, TOL=0.05) ===")
    print(f"Grid: {rows}x{cols}  Nodes: {summary['n_nodes']}  Edges: {summary['n_edges']}  AvgDeg: {summary['avg_degree']}")
    print("Status:", summary["status"], " | Objective:", round(summary["objective"], 6))
    print("Vars: x =", summary["n_x_vars"], " y =", summary["n_y_vars"], " | Constraints:", summary["n_constraints"])
    print("\nSpecies counts (solution) vs exact targets:")
    for name in space.species:
        T = TARGETS_SMALL[name]
        c = counts.get(name, 0)
        L = summary["bands"][name]["L"]; U = summary["bands"][name]["U"]
        print(f"  {name:25s}: {c:3d}  (target {T}, band [{L},{U}])")

        # Optional quick plot of final assignment (requires matplotlib)
    try:
        plot_solution(space,
                    assignment=assignment,   # dict node->species_idx
                    title="MILP assignment (reduced grid)",
                    spacing=1.0,
                    show_edges=False,
                    point_size_occ=28)
    except Exception as e:
        print("Plot skipped:", e)


