# Período de retorno univariado

## Contexto

El pipeline (`01_copulas.R`) ya calcula, en el PASO 13, el período de retorno
combinado (co-occurrence) de cada cópula final ajustada, análogo a la Figura 7
de Chen et al. 2024: `T(x,y) = N / (n·(1 - F_X(x) - F_Y(y) + C(F_X(x),F_Y(y))))`,
usando el objeto `mvdc` (cópula + ambas marginales) ya ajustado. Esto produce
una grilla 2D con isolíneas de T y un PNG por estación+par de variables.

Falta el análogo univariado: para cada variable individual (`intensidad`,
`magnitud`, `duracion`) ya se elige, en el PASO 6, la distribución que mejor
ajusta por estación (`mejor.ajuste.univariado.x.ubic.var`: `distribucion` +
`parametros`), pero nunca se traduce eso a una curva de período de retorno ni
a un gráfico. Este documento especifica cómo agregar ese cálculo.

## Alcance

- Nuevo PASO en `01_copulas.R` (PASO 14, entre el actual PASO 13 y el cierre
  "Finalizar script"), que calcula y grafica el período de retorno univariado
  para cada combinación estación+variable ya presente en
  `mejor.ajuste.univariado.x.ubic.var` — es decir, las 3 variables usadas en
  `variables_copulas` (`intensidad`, `magnitud`, `duracion`), sin importar en
  qué pares se usen luego.
- No modifica el cálculo bivariado existente (PASO 13) ni requiere que haya
  cópulas ajustadas — solo depende de PASO 6 (mejor ajuste univariado) y de
  `eventos_completos` (ya calculado en PASO 4), ambos disponibles en ese punto
  del script.
- Fuera de alcance: cambiar la fórmula/gráfico del período de retorno
  bivariado existente, o exponer período de retorno univariado para variables
  fuera de `variables_copulas` (`minimo`/`maximo`).

## Fórmula

Igual estructura que el caso bivariado, pero con una sola marginal:

```
T(x) = N / (n · (1 - F(x)))
```

donde:
- `F` es la CDF de la distribución ganadora para esa estación+variable
  (`do.call(paste0("p", distribucion), args = c(list(q = x), parametros))`).
- `N` = extensión del registro en años = `diff(range(fecha_inicio)) / 365.25`
  sobre los eventos observados (`tipo_serie == "observada"`) de esa
  estación+variable en `eventos_completos`.
- `n` = cantidad de eventos observados de esa estación+variable.

Misma fuente de `N`/`n`/valores observados que ya usa `CalcularPeriodoRetornoUC`
para el caso bivariado (filtra `eventos_completos` por `tipo_serie == "observada"`),
para mantener consistencia entre ambos cálculos.

## Grilla y gráfico

- Grilla 1D: `seq(min(x_obs), max(x_obs) + diff(range(x_obs)) * margen_grilla, length.out = resolucion_grilla)`.
- Config reusada de la sección `periodo_retorno` ya existente en el YAML de
  parámetros (`niveles_anios`, `resolucion_grilla`, `margen_grilla`) — no se
  agrega sección nueva.
- Gráfico por estación+variable: eje x = T (años, escala log), eje y = valor
  de la variable, curva teórica calculada sobre la grilla.
- Eventos observados superpuestos con posición de graficación empírica
  (Weibull): para el evento de rango `m` (1 = valor más alto, hasta `n`),
  `T_empirico = N · (n+1) / (n · m)`. Esto permite chequear visualmente si la
  distribución ajustada sigue realmente a los datos (si los puntos se alejan
  mucho de la curva teórica, es una señal de mal ajuste), igual que se usa de
  forma estándar en hidrología — no se calculan con la F teórica porque
  entonces caerían siempre exactos sobre la curva y no aportarían información
  de bondad de ajuste.
- Niveles de años (`niveles_anios`) marcados como líneas de referencia
  verticales en el gráfico, análogo a las curvas de nivel del caso bivariado.

## Archivos generados

- Tabla RDS, una fila por estación+variable, con columnas: `!!id_column`,
  `variable`, `distribucion`, `N`, `n`, `archivo_png`, `grilla` (lista, con
  la grilla 1D calculada). Guardada en la nueva clave de archivo
  `copulas.periodo_retorno_univariado`, agregada tanto a
  `data/configuracion_archivos_utilizados.yml` como a
  `data/configuracion_archivos_utilizados_frio_calor.yml`:
  ```yaml
  periodo_retorno_univariado: "output/periodo_retorno_univariado_<*idc>.rds"
  ```
- Un PNG por estación+variable:
  `output/periodo_retorno_univariado_<id_ubicacion>_<variable>.png`.

## Implementación (funciones nuevas)

- `lib/funciones_periodo_retorno.R`:
  - `CalcularGrillaPeriodoRetornoUV(distribucion, parametros, N, n, grid_x)`
    — análoga a `CalcularGrillaPeriodoRetorno` pero 1D.
  - `GraficarPeriodoRetornoUV(grilla, niveles_anios, x_obs, N, n, nombre_x, titulo, archivo_png)`
    — grafica la curva T-vs-valor + puntos observados en posición Weibull.
- `lib/funciones_worker.R`:
  - `CalcularPeriodoRetornoUV(input.value, script, mejores.ajustes.univariados, eventos.completos, niveles.anios, resolucion.grilla, margen.grilla, dir.salida.png)`
    — worker análogo a `CalcularPeriodoRetornoUC`, una fila por
    estación+variable.
- `01_copulas.R`: nuevo bloque PASO 14, calcado del patrón de tarea
  distribuida del PASO 13 (`Task$new` + `task$run` + logueo de errores +
  `saveRDS`), iterando sobre `mejor.ajuste.univariado.x.ubic.var` en vez de
  sobre cópulas finales.

## Manejo de errores

Mismo patrón defensivo que los fixes de `fc17f05` (no abortar el script ante
un fallo puntual): si el cálculo de `F`/la grilla falla para una
estación+variable específica (parámetros fuera de dominio, distribución sin
función `p*` disponible, etc.), esa fila se degrada a resultado `NA` (sin
`archivo_png` ni `grilla`) y el resto de la tarea distribuida continúa. Un
fallo aislado no debe interrumpir las demás combinaciones estación+variable,
siguiendo el mismo criterio ya aplicado al resto del pipeline.

## Testing

- Correr `pruebas/correr_18_escenarios.sh` (o una corrida single-tipo) y
  verificar que se generan los 3×n_estaciones PNGs univariados esperados por
  corrida, sin afectar los resultados existentes del PASO 13.
- Inspección visual de al menos un par de gráficos (uno cálido, uno frío) para
  confirmar que los puntos Weibull caen razonablemente cerca de la curva
  teórica.
- Confirmar que un fallo forzado en una combinación puntual (p.ej. estación
  con muy pocos eventos) no aborta el resto de la corrida.
