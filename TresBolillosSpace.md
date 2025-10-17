# Guía rápida de `TresBolillosSpace`

Esta clase modela un arreglo hexagonal de sitios apto para simulaciones de "tres bolillos". Se encarga de generar la geometría, asignar especies y calcular información útil para formulaciones de optimización. A continuación se describe cómo crear un espacio, inicializarlo y aprovechar sus métodos más importantes.

## Crear un espacio

La forma más directa es instanciar la clase con el modo rectangular por defecto:

```python
from TresBolillosSpace import TresBolillosSpace

espacio = TresBolillosSpace()
```

### Parámetros principales

- `mode`: `"rect"` (rejilla rectangular) o `"disk"` (aproximación a disco hexagonal).
- `rows`, `cols`: dimensiones de la rejilla cuando `mode="rect"`. El número total de sitios es `rows * cols`.
- `n_sites`: número de nodos cuando `mode="disk"`.
- `p_occupied`: probabilidad de que un sitio esté ocupado al muestrear un estado inicial.
- `species`, `species_probs`: listas opcionales para personalizar especies y sus probabilidades.
- `seed`: semilla reproducible para los muestreos aleatorios.

También existen constructores de conveniencia:

```python
espacio_rect = TresBolillosSpace.from_rect(rows=10, cols=20)
espacio_disk = TresBolillosSpace.from_disk(n_sites=300)
```

## Inicializar ocupación y especies

Para asignar ocupaciones y especies iniciales con base en las probabilidades configuradas:

```python
espacio.sample_initial()
print(espacio.y_init.sum())        # número de sitios ocupados
print(espacio.counts)              # conteo por especie
```

Puedes proporcionar una semilla específica en cada muestreo (`sample_initial(seed=1234)`) para obtener resultados controlados.

## Integración con PuLP

Si necesitas crear una formulación de programación lineal entera, obtén la información estructural con `to_pulp`:

```python
N, edges, species = espacio.to_pulp()
```

Para fijar en el modelo las especies que ya están sembradas, usa el método estático `add_preplanted_constraints_pulp(model, x, N, species_init, n_species)`. Este método recorre los nodos y agrega igualdades de asignación a la instancia de PuLP recibida.

## Visualización

Si `matplotlib` está disponible, puedes trazar el estado del espacio:

```python
espacio.sample_initial(seed=42)
espacio.plot(spacing=1.2, show_edges=False)
```

Los puntos vacíos y ocupados se diferencian en color y tamaño; las especies ocupadas se muestran con una barra de color.

## Cálculo de cuotas

Cuando se trabajan cuotas por especie, existen utilidades para resumir la información:

- `pre_counts()`: devuelve un diccionario con los conteos actuales por especie.
- `quota_bands_total(targets, tol=0.05)`: calcula bandas inferior/superior alrededor de cada objetivo `T` considerando una tolerancia relativa `tol`.
- `remaining_bands(targets, tol=0.05)`: combina la información anterior y entrega métricas de viabilidad (`free_slots`, `sum_L_rem`, etc.).

Estas funciones permiten validar rápidamente si las metas de reforestación son alcanzables con los espacios libres que quedan.

## Resumen

1. Instancia `TresBolillosSpace` con la geometría deseada.
2. Opcional: ajusta especies, probabilidades y semilla.
3. Llama a `sample_initial()` para obtener una distribución inicial.
4. Usa `to_pulp()` y `add_preplanted_constraints_pulp()` para integrarlo en modelos de optimización.
5. Apóyate en `plot()` y las funciones de cuotas para analizar resultados.

Con esta guía puedes comenzar a experimentar con diferentes configuraciones de plantación y evaluar rápidamente su viabilidad.
