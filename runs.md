# Run Configuration Summary 1 

No feasible solution found

Lower bound:                    -418.039

Enumerated nodes:               0

Total iterations:               137712

Time (CPU seconds):             7183.92

Time (Wallclock seconds):       7451.44

## Problem data
- **Space/class:** `TresBolillosSpace.from_rect(rows=14, cols=47, seed=42)`
- **Pre-plant sampling:** `space.sample_initial(seed=43)` (locks some nodes)
- **Species (|S| = 10) & targets (sum = 658):**
  - Agave lechuguilla: 42  
  - Agave salmiana: 196  
  - Agave scabra: 42  
  - Agave striata: 42  
  - Opuntia cantabrigiensis: 49  
  - Opuntia engelmannii: 38  
  - Opuntia robusta: 73  
  - Opuntia streptacantha: 64  
  - Prosopis laevigata: 86  
  - Yucca filifera: 26
- **Not-survival probabilities** $P_{\text{not}}(i)$. Survival: $s_i = 1 - P_{\text{not}}(i)$.

## Decision variables
- $x_{u,i}\in\{0,1\}$: node $u$ assigned to species $i$
- $y_{uvij}\in\{0,1\}$: linearization of $x_{u,i}\cdot x_{v,j}$ for each edge $(u,v)$ and species pair $(i,j)$ with $C_{ij}>0$

## Objective (minimize)
$$
\min\;
\alpha \sum_{(u,v)\in E}\sum_{i\in S}\sum_{j\in S} C_{ij}\, y_{uvij}
\;-\;
\beta \sum_{u\in N}\sum_{i\in S} s_i\, x_{u,i}
$$

- Weights: `ALPHA_COMP = 1.0`, `BETA_SURV = 1.0`  
- $s_i=1-P_{\text{not}}(i)$ (expected survival reward)  
- $C$ = competition matrix (penalizes neighboring assignments)

## Constraints
1. **One species per node**  
   $\sum_{i\in S} x_{u,i} = 1 \quad \forall u\in N$

2. **Pre-planted locks**  
   If $\mathrm{species\_init}[u]=i^\*\ge 0$: $x_{u,i^\*}=1$ and $x_{u,j}=0\;\forall j\neq i^\*$

3. **Species totals within ±5% of original targets**  
   $L_i=\lfloor(1-0.05)T_i\rfloor,\; U_i=\lceil(1+0.05)T_i\rceil$ and  
   $L_i \le \sum_{u\in N} x_{u,i} \le U_i \quad \forall i\in S$

4. **McCormick linearization for neighbor pairs** (only where $C_{ij}>0$)  
   $y_{uvij} \le x_{u,i},\;\; y_{uvij} \le x_{v,j},\;\; y_{uvij} \ge x_{u,i}+x_{v,j}-1$

## Competition matrix
`C = make_competition_matrix(seed=1234, sparsity=0.10, diag_bonus=1.0, offdiag_scale=1.0)`

## Solver (PuLP + CBC)
```python
solver = pl.PULP_CBC_CMD(
    msg=1,
    options=['log','5','threads','8'],
)
status = model.solve(solver)



# El bueno 

Result - Stopped on time limit

Objective value:                2.27055000
Lower bound:                    -18.092
Gap:                            1.13
Enumerated nodes:               1919
Total iterations:               158123
Time (CPU seconds):             171.45
Time (Wallclock seconds):       173.19

Option for printingOptions changed from normal to all
Total time (CPU seconds):       171.50   (Wallclock seconds):       173.26


[Scoring with GA metrics]
[MILP] competencia_avg = 1.006667   diversidad_dev = 0.333333
[BASE] competencia_avg = 0.100000   diversidad_dev = 1.200000

=== REDUCED SPACE RUN (all species, proportional targets, TOL=0.05) ===
Grid: 5x6  Nodes: 30  Edges: 69  AvgDeg: 4.6
Status: Optimal  | Objective: 2.27055
Vars: x = 300  y = 6900  | Constraints: 20810

Species counts (solution) vs exact targets:
  Agave lechuguilla        :   3  (target 2, band [1,3])
  Agave salmiana           :   8  (target 9, band [8,10])
  Agave scabra             :   2  (target 2, band [1,3])
  Agave striata            :   2  (target 2, band [1,3])
  Opuntia cantabrigiensis  :   3  (target 2, band [1,3])
  Opuntia engelmannii      :   2  (target 2, band [1,3])
  Opuntia robusta          :   3  (target 3, band [2,4])
  Opuntia streptacantha    :   2  (target 3, band [2,4])
  Prosopis laevigata       :   3  (target 4, band [3,5])
  Yucca filifera           :   2  (target 1, band [0,2])