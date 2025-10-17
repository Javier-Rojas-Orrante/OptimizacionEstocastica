# -*- coding: utf-8 -*-
from typing import Dict, List, Tuple, Optional
import math
import numpy as np

# ---------------- TresBolillosSpace (as provided) ----------------
try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None

class TresBolillosSpace:
    DEFAULT_SPECIES = [
        "Agave lechuguilla",
        "Agave salmiana",
        "Agave scabra",
        "Agave striata",
        "Opuntia cantabrigiensis",
        "Opuntia engelmannii",
        "Opuntia robusta",
        "Opuntia streptacantha",
        "Prosopis laevigata",
        "Yucca filifera",
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
        self.species_init = None  # np.ndarray (-1 o índice)
        self.counts = None        # dict nombre->conteo

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
        raise ValueError("mode debe ser 'rect' o 'disk'.")

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
            raise ValueError("La suma de probabilidades de especies debe ser > 0.")
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
        try:
            import pulp as pl  # noqa: F401
        except Exception:
            raise RuntimeError("Necesitas 'pulp' instalado.")
        for u in N:
            i_star = int(species_init[u])
            if i_star >= 0:
                model += x[(u, i_star)] == 1
                for j in range(n_species):
                    if j != i_star:
                        model += x[(u, j)] == 0

    def _positions_rect(self, spacing=1.0):
        dx = spacing * 0.5
        dy = spacing * math.sqrt(3) / 2
        pts = []
        for r in range(self.rows):
            for c in range(self.cols):
                x = c * spacing + (dx if r % 2 else 0.0)
                y = r * dy
                pts.append((x, y))
        return np.array(pts)

    def _positions_disk(self, spacing=1.0):
        def hex_count(R): return 1 + 3*R*(R+1)
        R = 0
        while hex_count(R) < self.n:
            R += 1
        coords = []
        for q in range(-R, R+1):
            rmin = max(-R, -q-R)
            rmax = min(R, -q+R)
            for r in range(rmin, rmax+1):
                coords.append((q, r))
        coords = coords[:self.n]
        pts = []
        rt3_2 = math.sqrt(3)/2
        for (q, r) in coords:
            x = q + 0.5*r
            y = rt3_2 * r
            pts.append((spacing * x, spacing * y))
        return np.array(pts)

    def plot(self, spacing=1.0, show_edges=True, figsize=(12, 5),
             point_size_empty=10, point_size_occ=28):
        if plt is None:
            raise RuntimeError("Matplotlib no está disponible.")
        from matplotlib.colors import BoundaryNorm
        pos = self._positions_rect(spacing) if self.mode == "rect" else self._positions_disk(spacing)
        fig, ax = plt.subplots(figsize=figsize)
        if show_edges:
            for (u, v) in self.edges:
                x1, y1 = pos[u]; x2, y2 = pos[v]
                ax.plot([x1, x2], [y1, y2], lw=0.4, alpha=0.25, color="#6a9fb5")
        if self.y_init is None or self.species_init is None:
            ax.scatter(pos[:, 0], pos[:, 1], s=18, alpha=0.9, color="#4f83cc")
        else:
            empty = (self.y_init == 0)
            if np.any(empty):
                ax.scatter(pos[empty, 0], pos[empty, 1], s=point_size_empty, alpha=0.5, color="#b0c4de")
            occ = (self.y_init == 1)
            if np.any(occ):
                cmap = plt.get_cmap("tab10")
                boundaries = np.arange(-0.5, len(self.species) + 0.5, 1.0)
                norm = BoundaryNorm(boundaries, cmap.N)
                sc = ax.scatter(pos[occ, 0], pos[occ, 1], s=point_size_occ,
                                c=self.species_init[occ], cmap=cmap, norm=norm, alpha=0.95, edgecolors="none")
                cbar = fig.colorbar(sc, ax=ax, ticks=np.arange(len(self.species)), fraction=0.03, pad=0.02)
                cbar.ax.set_yticklabels(self.species)
                cbar.set_label("Especie", rotation=90)
                cbar.ax.tick_params(labelsize=9)
        ax.set_aspect("equal")
        title = f"Tres bolillos — {self.mode} — nodos={self.n}"
        if self.y_init is not None:
            title += f" | ocupados={int(self.y_init.sum())}"
        ax.set_title(title)
        ax.axis("off")
        fig.tight_layout()
        plt.show()

    # ===== cuotas / bandas =====
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
                infeasible.append((name, "pre={} > U={}".format(P, bu["U"])))
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
# ----------------------------------------------------------------

import pulp as pl

# === Original full targets and mortality ===
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

# ---------- helpers ----------
def make_scaled_targets(original_targets: Dict[str, int], total_sites: int, subset: List[str] | None = None) -> Dict[str, int]:
    items = [(k, original_targets[k]) for k in (subset if subset else original_targets.keys())]
    total_orig = sum(v for _, v in items)
    if total_orig <= 0:
        raise ValueError("Original targets must sum > 0")
    floats = [(k, v * total_sites / total_orig) for k, v in items]
    rounded = {k: int(round(x)) for k, x in floats}
    drift = total_sites - sum(rounded.values())
    # adjust by residual magnitude
    residuals = sorted([(k, x - round(x)) for k, x in floats], key=lambda t: -abs(t[1]))
    idx = 0
    while drift != 0 and residuals:
        k, _ = residuals[idx % len(residuals)]
        if drift > 0:
            rounded[k] += 1
            drift -= 1
        else:
            if rounded[k] > 0:
                rounded[k] -= 1
                drift += 1
        idx += 1
    return rounded

def subset_dict(d: Dict[str, float], keys: List[str]) -> Dict[str, float]:
    return {k: d[k] for k in keys}

def build_random_competition(species: list, seed: int = 123) -> np.ndarray:
    rng = np.random.default_rng(seed)
    S = len(species)
    M = rng.uniform(0.1, 1.0, size=(S, S))
    C = 0.5 * (M + M.T)
    for i in range(S):
        C[i, i] = max(C[i, i], 1.05)  # stronger same-species penalty
    return C

# ---------- MILP ----------
def build_and_solve(
    space: TresBolillosSpace,
    targets: Dict[str, int],
    p_not_survive: Dict[str, float],
    tol: float = 0.05,
    w_surv: float = 1.0,
    w_comp: float = 1.0,
    comp_seed: int = 123,
    time_limit_sec: int | None = 60,
) -> Tuple[Dict[int, int], Dict[str, int], float]:
    if space.y_init is None or space.species_init is None:
        space.sample_initial()  # uses space.seed if provided

    N, E, species = space.to_pulp()
    n, m = len(N), len(E)
    S = len(species)

    # map inputs to index order
    T = np.array([int(targets[s]) for s in species], dtype=int)
    p_die = np.array([float(p_not_survive[s]) for s in species], dtype=float)

    # bands and feasibility quick-checks
    bands = space.quota_bands_total(targets, tol=tol)
    rem, checks = space.remaining_bands(targets, tol=tol)
    if checks["species_over_cap"]:
        over = ", ".join(f"{nm} ({msg})" for nm, msg in checks["species_over_cap"])
        raise ValueError("Infeasible: pre-planted counts exceed cap: " + over)
    if not checks["min_feasible_given_free?"] or not checks["max_feasible_given_free?"]:
        raise ValueError(
            f"Infeasible: free={checks['free_slots']}, "
            f"min_needed={checks['sum_L_rem']}, max_allow={checks['sum_U_rem']}"
        )

    C = build_random_competition(species, seed=comp_seed)

    model = pl.LpProblem("ReforestationAssignment_Small", pl.LpMinimize)

    # x[u,i]
    x = pl.LpVariable.dicts("x", ((u, i) for u in N for i in range(S)), lowBound=0, upBound=1, cat="Binary")
    # y[u,v,i,j]
    y = pl.LpVariable.dicts(
        "y",
        ((u, v, i, j) for (u, v) in E for i in range(S) for j in range(S)),
        lowBound=0, upBound=1, cat="Binary"
    )

    # one species per node
    for u in N:
        model += pl.lpSum(x[(u, i)] for i in range(S)) == 1

    # honor pre-planted
    TresBolillosSpace.add_preplanted_constraints_pulp(model, x, N, space.species_init, S)

    # quotas
    for i, sname in enumerate(species):
        L_i = bands[sname]["L"]
        U_i = bands[sname]["U"]
        model += pl.lpSum(x[(u, i)] for u in N) >= L_i
        model += pl.lpSum(x[(u, i)] for u in N) <= U_i

    # AND linearization
    for (u, v) in E:
        for i in range(S):
            for j in range(S):
                model += y[(u, v, i, j)] <= x[(u, i)]
                model += y[(u, v, i, j)] <= x[(v, j)]
                model += y[(u, v, i, j)] >= x[(u, i)] + x[(v, j)] - 1

    # objective
    scale_surv = 1.0 / max(1, n)
    scale_comp = 1.0 / max(1, m)
    obj_surv = scale_surv * pl.lpSum(p_die[i] * x[(u, i)] for u in N for i in range(S))
    obj_comp = scale_comp * pl.lpSum(C[i, j] * y[(u, v, i, j)] for (u, v) in E for i in range(S) for j in range(S))
    model += w_surv * obj_surv + w_comp * obj_comp

    solver = pl.PULP_CBC_CMD(msg=1, timeLimit=time_limit_sec) if time_limit_sec else pl.PULP_CBC_CMD(msg=1)
    model.solve(solver)

    assign_idx: Dict[int, int] = {}
    for u in N:
        chosen = None
        for i in range(S):
            if pl.value(x[(u, i)]) > 0.5:
                chosen = i
                break
        if chosen is None:
            chosen = int(np.argmax([pl.value(x[(u, i)]) for i in range(S)]))
        assign_idx[u] = chosen

    counts_by_name: Dict[str, int] = {species[i]: 0 for i in range(S)}
    for u, i in assign_idx.items():
        counts_by_name[species[i]] += 1

    return assign_idx, counts_by_name, pl.value(model.objective)

# --------------------------- TINY TEST MAIN ---------------------------
if __name__ == "__main__":
    # Tiny grid: 6x8 = 48 nodes
    rows, cols = 6, 8
    total_sites = rows * cols

    # Use 4 species (subset of your 10)
    species_subset = [
        "Agave salmiana",
        "Opuntia robusta",
        "Prosopis laevigata",
        "Yucca filifera",
    ]

    # Build tiny space with NO pre-planted nodes (feasibility safe)
    space = TresBolillosSpace.from_rect(
        rows=rows, cols=cols,
        species=species_subset,
        species_probs=subset_dict(TresBolillosSpace.DEFAULT_SPECIES_PROBS, species_subset),
        p_occupied=0.0,
        seed=1
    )
    space.sample_initial(seed=1)

    # Scale targets exactly to 48 nodes for the chosen subset
    SPECIES_TARGETS_TINY = make_scaled_targets(SPECIES_TARGETS_FULL, total_sites, subset=species_subset)
    PROB_NOT_SURVIVE_TINY = subset_dict(PROB_NOT_SURVIVE, species_subset)

    # Solve
    assignment, counts, obj = build_and_solve(
        space=space,
        targets=SPECIES_TARGETS_TINY,
        p_not_survive=PROB_NOT_SURVIVE_TINY,
        tol=0.05,
        w_surv=1.0,
        w_comp=1.0,
        comp_seed=2025,
        time_limit_sec=60
    )

    print("\n[TINY] Objective value:", obj)
    print("[TINY] Assigned counts (±5% bands):")
    for k in species_subset:
        T = SPECIES_TARGETS_TINY[k]
        L, U = math.floor(0.95*T), math.ceil(1.05*T)
        print(f" - {k:24s}: {counts[k]:3d}  (target {T}, band [{L},{U}])")

    # Optional: quick plot
    try:
        if plt is not None:
            backup_y, backup_species = space.y_init.copy(), space.species_init.copy()
            space.y_init = np.ones(space.n, dtype=int)
            space.species_init = np.array([assignment[u] for u in range(space.n)], dtype=int)
            space.plot(spacing=1.0, show_edges=False, figsize=(8, 4), point_size_empty=0, point_size_occ=28)
            space.y_init, space.species_init = backup_y, backup_species
        else:
            print("Plot skipped: Matplotlib no disponible.")
    except Exception as e:
        print("Plot skipped:", e)
