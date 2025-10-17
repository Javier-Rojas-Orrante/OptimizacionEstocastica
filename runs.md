# Run Configuration Summary 1 

No feasible solution found
Lower bound:                    -418.039
Enumerated nodes:               0
Total iterations:               137712
Time (CPU seconds):             7183.92
Time (Wallclock seconds):       7451.44

## Problem data
- **Space/class:** `TresBolillosSpace.from_rect(rows=14, cols=47, seed=42)`
- **Pre-plant sampling:** `space.sample_initial(seed=43)` (locks some nodes to specific species)
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
- **Not-survival probabilities** \(P_{\text{not}}(i)\). Survival: \(s_i = 1 - P_{\text{not}}(i)\).

## Decision variables
- \(x_{u,i}\in\{0,1\}\): node \(u\) assigned to species \(i\)
- \(y_{uvij}\in\{0,1\}\): linearization of \(x_{u,i}\cdot x_{v,j}\) for each edge \((u,v)\) and species pair \((i,j)\) with \(C_{ij}>0\)

## Objective (minimize)
$$
\min\;\; 
\alpha\sum_{(u,v)\in E}\sum_{i\in S}\sum_{j\in S} C_{ij}\,y_{uvij}
\;-\;
\beta\sum_{u\in N}\sum_{i\in S} s_i\,x_{u,i}
$$

- Weights: `ALPHA_COMP = 1.0`, `BETA_SURV = 1.0`  
- \(s_i=1-P_{\text{not}}(i)\) (expected survival reward)  
- \(C\) = competition matrix (penalizes neighboring assignments)

## Constraints
1. **One species per node**  
   \(\sum_{i\in S} x_{u,i}=1\quad\forall u\in N\)

2. **Pre-planted locks**  
   If \(\text{species\_init}[u]=i^\*\ge 0\): \(x_{u,i^\*}=1\) and \(x_{u,j}=0\;\forall j\neq i^\*\)

3. **Species totals within ±5% of original targets**  
   \(L_i=\lfloor(1-0.05)T_i\rfloor,\; U_i=\lceil(1+0.05)T_i\rceil\) and  
   \(L_i \le \sum_{u\in N} x_{u,i} \le U_i\quad\forall i\in S\)

4. **McCormick linearization for neighbor pairs** (only where \(C_{ij}>0\))  
   \(y_{uvij}\le x_{u,i},\;\; y_{uvij}\le x_{v,j},\;\; y_{uvij}\ge x_{u,i}+x_{v,j}-1\)

## Competition matrix \(C=\texttt{make\_competition\_matrix}(\cdot)\)
- Params: `seed=1234`, `sparsity=0.10`, `diag_bonus=1.0`, `offdiag_scale=1.0`  
- Symmetric, nonnegative; diagonal boosted to discourage same-species clumping

## Solver (PuLP + CBC)
```python
solver = pl.PULP_CBC_CMD(
    msg=1,
    options=['log','5','threads','8'],
)
status = model.solve(solver)
