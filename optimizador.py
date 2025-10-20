def fitnessCompetencia(plantas, matrizCompetencia):
    """
    Calcula la competencia total entre las distribución de plantas en una
    matriz, se asume una distribución de tres bolillos. La competencia se calcula segun
    la matriz de competencia.
    """
    competencia = 0

    rows = len(plantas)
    cols = len(plantas[0])

    for i in range(rows):
        for j in range(cols):
            especie_nodo = plantas[i, j]
            
            # Skip empty nodes (-1)
            if especie_nodo == -1:
                continue
                
            if rows % 2 == 0: # Distribución de tres bolillos - filas pares
                if i == 0: # Primera fila
                    if j == 0: # Primera columna
                        if plantas[i, j+1] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i, j+1]]
                        if plantas[i+1, j] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i+1, j]]
                        if plantas[i+1, j+1] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i+1, j+1]]
                    elif j == cols-1: # Ultima columna (CORREGIDO: era rows-1)
                        if plantas[i, j-1] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i, j-1]]
                        if plantas[i+1, j] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i+1, j]]
                    else:
                        if plantas[i, j+1] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i, j+1]]
                        if plantas[i, j-1] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i, j-1]]
                        if plantas[i+1, j] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i+1, j]]
                        if plantas[i+1, j+1] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i+1, j+1]]   
                elif i == rows-1: # Ultima fila
                    if j == 0: # Primera columna
                        if plantas[i-1, j] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i-1, j]]
                        if plantas[i, j+1] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i, j+1]]
                    elif j == cols-1: # Ultima columna
                        if plantas[i, j-1] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i, j-1]]
                        if plantas[i-1, j-1] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i-1, j-1]]
                        if plantas[i-1, j] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i-1, j]]
                    else:
                        if plantas[i-1, j-1] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i-1, j-1]]
                        if plantas[i-1, j] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i-1, j]]
                        if plantas[i, j-1] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i, j-1]]
                        if plantas[i, j+1] >= 0:
                            competencia += matrizCompetencia[especie_nodo][plantas[i, j+1]]
                else: # Filas intermedias
                    if i % 2 == 0: # Fila par
                        if j == 0: # Primera columna
                            if plantas[i-1, j] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i-1, j]]
                            if plantas[i-1, j+1] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i-1, j+1]]
                            if plantas[i+1, j] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i+1, j]]
                            if plantas[i+1, j+1] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i+1, j+1]]
                            if plantas[i, j+1] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i, j+1]]
                        elif j == cols-1: # Ultima columna (CORREGIDO: era rows-1)
                            if plantas[i-1, j] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i-1, j]]
                            if plantas[i+1, j] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i+1, j]]
                            if plantas[i, j-1] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i, j-1]]
                        else:
                            if plantas[i-1, j] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i-1, j]]
                            if plantas[i-1, j+1] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i-1, j+1]]
                            if plantas[i+1, j] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i+1, j]]
                            if plantas[i+1, j+1] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i+1, j+1]]
                            if plantas[i, j+1] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i, j+1]]
                            if plantas[i, j-1] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i, j-1]] 
                    else: # Fila impar
                        if j == 0: # Primera columna
                            if plantas[i-1, j] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i-1, j]]
                            if plantas[i+1, j] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i+1, j]]
                            if plantas[i, j+1] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i, j+1]]
                        elif j == cols-1: # Ultima columna (CORREGIDO: era rows-1)
                            if plantas[i-1, j] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i-1, j]]
                            if plantas[i-1, j-1] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i-1, j-1]]
                            if plantas[i+1, j] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i+1, j]]
                            if plantas[i+1, j-1] >= 0:  # CORREGIDO: era duplicado
                                competencia += matrizCompetencia[especie_nodo][plantas[i+1, j-1]]
                            if plantas[i, j-1] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i, j-1]]
                        else:
                            if plantas[i-1, j] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i-1, j]]
                            if plantas[i-1, j-1] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i-1, j-1]]
                            if plantas[i+1, j] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i+1, j]]
                            if plantas[i+1, j-1] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i+1, j-1]]
                            if plantas[i, j+1] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i, j+1]]
                            if plantas[i, j-1] >= 0:
                                competencia += matrizCompetencia[especie_nodo][plantas[i, j-1]]                         
                                                  
    return competencia

import numpy as np

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

import numpy as np
from pymoo.core.problem import Problem

class PlantDistributionProblem(Problem):
    """
    Problema de optimización multi-objetivo para distribución de plantas.
    
    Objetivos:
    1. Minimizar competencia entre plantas vecinas
    2. Minimizar desviación de distribución uniforme (maximizar diversidad)
    
    Representación:
    - Cada solución es un vector 1D de longitud n_nodes
    - Cada gen representa la especie en ese nodo (0-9 para especies, -1 para vacío)
    """
    
    def __init__(self, space, competencia_matrix, targets, tol=0.05, enforce_constraints=False):
        """
        Parameters:
        -----------
        space : TresBolillosSpace
            Espacio de tres bolillos con configuración inicial
        competencia_matrix : np.ndarray
            Matriz de competencia entre especies (10x10)
        targets : dict
            Diccionario con objetivos por especie
        tol : float
            Tolerancia para bandas de objetivos (ej: 0.15 = ±15%)
        enforce_constraints : bool
            Si True, las bandas son restricciones duras. Si False, son suaves (penalizadas en objetivos)
        """
        self.space = space
        self.competencia_matrix = competencia_matrix
        self.targets = targets
        self.tol = tol
        self.rows = space.rows
        self.cols = space.cols
        self.n_nodes = space.n
        self.n_species = len(space.species)
        self.enforce_constraints = enforce_constraints
        
        # Calcular bandas permitidas por especie
        self.bands_rem, self.checks = space.remaining_bands(targets, tol=tol)
        
        # Nodos que ya tienen plantas (pre-sembrados) - estos no cambiarán
        self.fixed_nodes = np.where(space.y_init == 1)[0]
        self.fixed_species = space.species_init[self.fixed_nodes].copy()
        
        # Nodos vacíos que podemos modificar
        self.empty_nodes = np.where(space.y_init == 0)[0]
        self.n_variables = len(self.empty_nodes)
        
        print(f"Problema configurado:")
        print(f"  - Nodos totales: {self.n_nodes}")
        print(f"  - Nodos fijos (pre-sembrados): {len(self.fixed_nodes)}")
        print(f"  - Nodos a optimizar (vacíos): {self.n_variables}")
        print(f"  - Especies: {self.n_species}")
        print(f"  - Tolerancia de bandas: ±{tol*100:.0f}%")
        print(f"  - Restricciones: {'DURAS (obligatorias)' if enforce_constraints else 'SUAVES (penalizadas)'}")
        
        # Definir el problema
        n_constr = self.n_species if enforce_constraints else 0
        super().__init__(
            n_var=self.n_variables,  # Solo optimizamos los nodos vacíos
            n_obj=2,  # Dos objetivos: competencia y diversidad
            n_constr=n_constr,  # Restricciones opcionales
            xl=0,  # Límite inferior: especie 0
            xu=self.n_species - 1,  # Límite superior: especie 9
            type_var=int
        )
    
    def _evaluate(self, X, out, *args, **kwargs):
        """
        Evalúa las soluciones.
        
        X: array de shape (pop_size, n_variables)
           Cada fila es una solución con especies solo para nodos vacíos
        """
        pop_size = X.shape[0]
        f1 = np.zeros(pop_size)  # Competencia
        f2 = np.zeros(pop_size)  # Diversidad
        
        for i in range(pop_size):
            # Reconstruir la solución completa
            full_solution = np.full(self.n_nodes, -1, dtype=int)
            
            # Mantener nodos fijos
            full_solution[self.fixed_nodes] = self.fixed_species
            
            # Asignar especies a nodos vacíos desde X[i]
            full_solution[self.empty_nodes] = X[i]
            
            # Convertir a matriz para evaluación
            matrix = full_solution.reshape(self.rows, self.cols)
            matrix = np.flipud(matrix)  # Ajustar orientación
            
            # Calcular objetivos
            competencia_raw = fitnessCompetencia(matrix, self.competencia_matrix)
            diversidad_raw = fitnessDiversidad(matrix, self.n_species)
            
            # Normalizar competencia dividiéndola por el número de nodos ocupados
            # Esto hace que ambos objetivos estén en escalas comparables
            n_occupied = (full_solution >= 0).sum()
            f1[i] = competencia_raw / max(n_occupied, 1)  # Competencia promedio por planta
            f2[i] = diversidad_raw  # Diversidad ya está en escala pequeña
            
            # Si las restricciones están desactivadas, penalizar violaciones en los objetivos
            if not self.enforce_constraints:
                species_counts = np.bincount(full_solution[full_solution >= 0], 
                                            minlength=self.n_species)
                
                penalty = 0
                for s_idx, species_name in enumerate(self.space.species):
                    count = species_counts[s_idx]
                    band = self.bands_rem[species_name]
                    L_total = band['pre'] + band['L_rem']
                    U_total = band['pre'] + band['U_rem']
                    
                    # Penalizar desviaciones de la banda
                    if count < L_total:
                        penalty += (L_total - count) ** 2
                    elif count > U_total:
                        penalty += (count - U_total) ** 2
                
                # Agregar penalización normalizada a ambos objetivos
                f1[i] += penalty * 0.5  # Penalización moderada en competencia
                f2[i] += penalty * 0.02  # Penalización suave en diversidad
        
        out["F"] = np.column_stack([f1, f2])
        
        # Si las restricciones están activadas, calcularlas
        if self.enforce_constraints:
            g = np.zeros((pop_size, self.n_species))
            
            for i in range(pop_size):
                full_solution = np.full(self.n_nodes, -1, dtype=int)
                full_solution[self.fixed_nodes] = self.fixed_species
                full_solution[self.empty_nodes] = X[i]
                
                species_counts = np.bincount(full_solution[full_solution >= 0], 
                                            minlength=self.n_species)
                
                for s_idx, species_name in enumerate(self.space.species):
                    count = species_counts[s_idx]
                    band = self.bands_rem[species_name]
                    L_total = band['pre'] + band['L_rem']
                    U_total = band['pre'] + band['U_rem']
                    
                    if count < L_total:
                        g[i, s_idx] = L_total - count
                    elif count > U_total:
                        g[i, s_idx] = count - U_total
                    else:
                        g[i, s_idx] = 0
            
            out["G"] = g

from pymoo.core.crossover import Crossover

class MidRowCrossover(Crossover):
    """
    Operador de cruza personalizado: cruza de un punto a la mitad de las filas.
    """
    def __init__(self):
        super().__init__(2, 2)  # 2 padres, 2 hijos
    
    def _do(self, problem, X, **kwargs):
        """
        X: array de shape que puede venir como (n_parents, n_matings, n_vars) O (n_matings, n_parents, n_vars)
        Debe retornar: array de shape (n_offsprings, n_matings, n_variables)
        """
        # DETECTAR LA ORIENTACIÓN CORRECTA
        # Si la primera dimensión es igual a n_parents (2), entonces viene como (n_parents, n_matings, n_vars)
        # y necesitamos transponerlo a (n_matings, n_parents, n_vars)
        if X.shape[0] == self.n_parents:
            X = np.transpose(X, (1, 0, 2))  # Intercambiar dimensiones 0 y 1
        
        n_matings, n_parents, n_vars = X.shape
        
        # Verificar que tenemos 2 padres
        if n_parents != 2:
            raise ValueError(f"Se esperaban 2 padres pero se recibieron {n_parents}")
        
        # Punto de cruza: mitad de las filas
        nodes_per_row = n_vars // problem.rows
        crossover_point = (problem.rows // 2) * nodes_per_row
        
        # Crear array de salida
        Y = np.zeros((self.n_offsprings, n_matings, n_vars), dtype=int)
        
        # Realizar crossover
        for i in range(n_matings):
            # Hijo 1: primera mitad del padre 1, segunda mitad del padre 2
            Y[0, i, :crossover_point] = X[i, 0, :crossover_point]
            Y[0, i, crossover_point:] = X[i, 1, crossover_point:]
            
            # Hijo 2: primera mitad del padre 2, segunda mitad del padre 1
            Y[1, i, :crossover_point] = X[i, 1, :crossover_point]
            Y[1, i, crossover_point:] = X[i, 0, crossover_point:]
        
        return Y
from pymoo.core.mutation import Mutation

class SpeciesMutation(Mutation):
    """
    Operador de mutación personalizado: cambia la especie en nodos con cierta probabilidad.
    """
    def __init__(self, prob=0.1):
        """
        Parameters:
        -----------
        prob : float
            Probabilidad de mutar cada gen (nodo)
        """
        super().__init__()
        self.prob = prob
    
    def _do(self, problem, X, **kwargs):
        """
        X: array de shape (pop_size, n_variables)
        """
        n_pop, n_vars = X.shape
        Y = X.copy()
        
        for i in range(n_pop):
            for j in range(n_vars):
                if np.random.random() < self.prob:
                    # Cambiar a una especie aleatoria diferente
                    current_species = Y[i, j]
                    available_species = [s for s in range(problem.n_species) if s != current_species]
                    Y[i, j] = np.random.choice(available_species)
        
        return Y

from pymoo.core.sampling import Sampling

class SmartInitialSampling(Sampling):
    """
    Muestreo inicial inteligente que respeta las bandas de objetivos.
    """
    def _do(self, problem, n_samples, **kwargs):
        """
        Genera n_samples soluciones iniciales válidas.
        """
        X = np.zeros((n_samples, problem.n_variables), dtype=int)
        
        for i in range(n_samples):
            # Calcular cuántas plantas de cada especie necesitamos agregar
            species_to_add = []
            
            for s_idx, species_name in enumerate(problem.space.species):
                band = problem.bands_rem[species_name]
                # Agregar entre L_rem y U_rem plantas de esta especie
                n_to_add = np.random.randint(band['L_rem'], band['U_rem'] + 1)
                species_to_add.extend([s_idx] * n_to_add)
            
            # Si no tenemos suficientes plantas, completar con aleatorias
            while len(species_to_add) < problem.n_variables:
                species_to_add.append(np.random.randint(0, problem.n_species))
            
            # Si tenemos demasiadas, recortar
            species_to_add = species_to_add[:problem.n_variables]
            
            # Mezclar aleatoriamente
            np.random.shuffle(species_to_add)
            X[i] = species_to_add
        
        return X
    
# Configurar el algoritmo NSGA-II
from pymoo.operators.selection.tournament import TournamentSelection

# Función de comparación para el torneo
def tournament_compare(pop, P, **kwargs):
    """
    Compara individuos basándose en su rank (dominancia de Pareto).
    Retorna el índice del mejor individuo en cada torneo.
    """
    n_tournaments = P.shape[0]
    S = np.zeros(n_tournaments, dtype=int)
    
    for i in range(n_tournaments):
        # Comparar todos los individuos en este torneo
        tournament_indices = P[i]
        
        # Obtener ranks usando get() - si no existe, usar CV (constraint violation)
        ranks = []
        for idx in tournament_indices:
            ind = pop[idx]
            
            # Intentar obtener rank
            if hasattr(ind, 'rank') and ind.rank is not None:
                rank_val = ind.rank
            elif ind.get("rank") is not None:
                rank_val = ind.get("rank")
            else:
                # Si no hay rank, usar CV (violación de restricciones)
                cv = ind.get("CV")
                if cv is not None:
                    rank_val = 1000 + cv
                else:
                    rank_val = 1000
            
            # Asegurar que sea un escalar
            if isinstance(rank_val, np.ndarray):
                rank_val = float(rank_val.item()) if rank_val.size == 1 else float(rank_val[0])
            else:
                rank_val = float(rank_val)
            
            ranks.append(rank_val)
        
        ranks = np.array(ranks, dtype=float)
        
        # Seleccionar el de menor rank (mejor)
        best_in_tournament = tournament_indices[np.argmin(ranks)]
        S[i] = best_in_tournament
    
    return S


