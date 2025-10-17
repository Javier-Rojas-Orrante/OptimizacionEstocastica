from typing import Dict, List, Tuple, Optional
import math
import numpy as np

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


