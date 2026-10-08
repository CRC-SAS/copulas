# Arquitectura del pipeline

Documento de referencia técnica de `01_copulas.R`: qué hace cada paso, cómo se
organiza el código y qué significa cada parámetro de configuración. Para
*correr* el pipeline, ver la [guía de ejecución](guia_ejecucion.md). Para el
contexto conceptual (SPI/SPEI, pentadas, generador estocástico) ver
[`guia_video.md`](guia_video.md).

> **Fuera de alcance de este repo:** la identificación de eventos y el
> generador estocástico de series sintéticas son procesos previos, corridos
> por fuera. Acá se **consume** su salida (un CSV de eventos, ver
> [formato del CSV](guia_ejecucion.md#4-el-csv-de-eventos)).

## 1. Flujo de `01_copulas.R`

El script es un único pipeline, dirigido por 3 archivos YAML. Los pasos están
numerados en el código (`# --- PASO n`):

| Paso | Qué hace | Función worker (`lib/funciones_worker.R`) |
|---|---|---|
| 1-3 | Cargar paquetes, leer los 3 YAML, cargar `lib/` e iniciar el log | — |
| 4 | Leer los eventos, filtrarlos (tipo, duración mínima, realizaciones), agregarles el evento ENSO vigente al inicio (`evento_enso`) y generar las series "perturbadas" (con ruido, para propagar incertidumbre) | — |
| 5 | Ajuste univariado de cada (ubicación, variable, distribución) sobre la serie observada | `AjusteUnivariadoUVD` |
| 6 | Mejor distribución por (ubicación, variable) | `MejorAjusteUnivariadoUV` |
| 7 | Aplicar la mejor distribución a cada serie perturbada | `AplicarMejorAjusteASeriesPerturbadas` |
| 8 | Estacionariedad (Box-Pierce/Ljung-Box, punto de cambio, Mann-Kendall, test empírico) | `DeterminarEstacionariedad` |
| 9 | Dependencia entre pares de variables (Kendall, Spearman, cópula empírica) | `DeterminarDependencia` |
| 10 | Ajuste de cópulas por familia (solo si ambas series son estacionarias y dependientes) | `AjustarCopulas` |
| 11 | Mejor familia de cópula por (ubicación, par de variables) | `MejorAjusteMultivariadoUC` |
| 12 | Ajuste multivariado final: cópula + marginales en un objeto `copula::mvdc` | `AplicarMejorAjusteACopulas` |
| 13 | Período de retorno **combinado** y gráfico de isolíneas (un PNG por ubicación y par de variables) | `CalcularPeriodoRetornoUC` |
| 14 | Período de retorno **univariado** y gráfico (un PNG por ubicación y variable) | `CalcularPeriodoRetornoUV` |

Cada paso paralelo se ejecuta con `Task$run()`, que reparte las filas de
`input.values` entre procesos (`doSNOW`/`foreach`), con un máximo de
`max.procesos`. Los errores de las tareas se acumulan y, si hay alguno, el
script aborta después de loguearlos.

**Selección de la mejor distribución (paso 6).** Compara las candidatas por
RMSE/IQR, CCC y test de cuantiles, más los tests KS, Anderson-Darling y
Cramér-von Mises (estos tres se omiten para las variables listadas en
`variables_sin_test_continuidad`). Una variable con menos de
`min_cantidad_valores_ajuste_univariado` eventos (30 por defecto) no se ajusta.

**Selección de la mejor cópula (paso 11).** Por AIC, BIC, MRE, RMSE,
validación cruzada y test de bondad de ajuste Sn.

### Período de retorno

Con `N` = extensión del registro en años (`diff(range(fecha_inicio)) / 365.25`
sobre los eventos observados), `n` = cantidad de eventos observados, `F`/`C` =
marginal y cópula ajustadas:

- Combinado (paso 13): `T(x,y) = N / (n · (1 − F_X(x) − F_Y(y) + C(F_X(x), F_Y(y))))`
  — análogo a la Figura 7 de Chen et al. 2024. Los puntos observados se pintan
  según el evento ENSO vigente al inicio de cada evento (se desactiva con
  `periodo_retorno.colorear_enso: false`).
- Univariado (paso 14): `T(x) = N / (n · (1 − F(x)))`. El gráfico superpone los
  eventos observados en posición de graficación de Weibull,
  `T_emp = N·(n+1) / (n·m)` (`m` = rango, 1 = el mayor), para chequear a ojo si
  la distribución ajustada sigue a los datos. Opcionalmente
  (`graficar_distribucion_ajuste`) se agrega un histograma con la densidad
  ganadora (análogo a la Figura 5 de Chen et al. 2024).

Si el cálculo falla para una combinación puntual (parámetro fuera de dominio,
cópula sin ajustar, etc.), esa fila queda en `NA` y el resto de la corrida
continúa.

## 2. Código

```
01_copulas.R                 → pipeline completo (pasos 1-14)
lib/R/                       → framework de ejecución
  Script.R                     clase R6 de la corrida: logging (futile.logger) a run/CalcCopulas.log
  Task.R                       clase R6 que corre un worker en paralelo, con su log por tarea
  Helpers.R                    utilidades (IdentificarIdColumn: detecta station_id / point_id / *_id)
lib/*.R                      → funciones de negocio
  funciones_worker.R             un wrapper por paso paralelo (tabla de arriba)
  funciones_ajustes_marginales.R ~25 ajustes univariados (L-momentos y máx. verosimilitud, vía lmomco)
  funciones_ajuste_distribuciones.R, ajuste_distribuciones.R  orquestación y CDF de las marginales
  funciones_bondad_ajuste.R      KS, Anderson-Darling, Cramér-von Mises, RMSE/IQR
  funciones_mejor_ajuste_univariado.R  selección de la mejor distribución
  funciones_test_estacionaridad.R / funciones_test_independencia.R
  funciones_ajuste_familias_copulas.R  Gumbel, Frank, AMH, Joe, Clayton, Normal, t
  funciones_bondad_ajuste_copulas.R    Sn, SnB, SnC, AIC, BIC, MRE/RMSE, validación cruzada
  funciones_mejor_ajuste_copula.R      selección de la mejor familia
  funciones_periodo_retorno.R    grillas y gráficos de período de retorno (bivariado y univariado)
  funciones_enso.R               lectura y clasificación del archivo ENSO
  funciones_auxiliares.R         AgregarRuido (series perturbadas), Chi-plot, K-plot
```

> `lib/funciones_bondad_ajuste.R` está en ISO-8859-1 (Latin-1), no UTF-8 como el
> resto. R lo carga sin problema, pero `grep` sin `-a` puede no encontrar
> coincidencias, y un editor que asuma UTF-8 puede corromper los acentos.

## 3. Datos y archivos

```
configuracion_copulas.yml        → rutas de trabajo y procesos paralelos (local, no versionado;
                                   se genera desde configuracion_copulas.yml.tmpl)
parametros_copulas.yml           → qué se analiza y con qué criterios
data/configuracion_archivos_utilizados.yml → nombres de los archivos que lee/escribe cada paso
data/input/                      → insumos: CSV de eventos y enso_roni.txt (no versionado)
data/partial/                    → resultados intermedios (.rds)
data/output/                     → resultados finales (.rds y .png)
run/                             → logs (CalcCopulas.log y uno por tarea) y .pid de la corrida
```

Resultados principales (`<id>` = `identificador_corrida`):

| Archivo | Contenido |
|---|---|
| `output/copulas_<id>.rds` | Resultado final: cópula + marginales elegidas, por ubicación y par de variables |
| `output/mejores_ajustes_univariados_<id>.rds` | Mejor distribución por ubicación y variable |
| `output/mejores_ajustes_multivariados_<id>.rds` | Mejor familia de cópula por ubicación y par |
| `output/periodo_retorno_<id>.rds` + `periodo_retorno_<ubic>_<x>_<y>.png` | Período de retorno combinado |
| `output/periodo_retorno_univariado_<id>.rds` + `periodo_retorno_univariado_<ubic>_<var>.png` | Período de retorno univariado |
| `partial/copulas_<id>.info` | Info de la corrida |

## 4. Referencia de parámetros

### `configuracion_copulas.yml`

| Clave | Qué es |
|---|---|
| `dir.base` | Raíz del repositorio (ruta absoluta) |
| `dir.run` | Carpeta de logs y `.pid` |
| `dir.lib` | Carpeta `lib/R/` (framework de ejecución) |
| `dir.data` | Carpeta `data/` |
| `max.procesos` | Máximo de procesos paralelos |

### `parametros_copulas.yml`

| Clave | Qué define |
|---|---|
| `ubicaciones` | Lista de `{id, nombre}` a procesar. Cada `id` debe existir en la columna `*_id` del CSV (para el tipo de evento elegido) |
| `configuraciones.eventos` | **Una única fila** `conf_id, indice, escala, distribucion, metodo_ajuste`; filtra el CSV por `conf_id`. El script aborta si hay más de una |
| `variables_copulas` | Pares `variable_x`/`variable_y`. Deben ser columnas del CSV: `duracion`, `intensidad`, `magnitud`, `minimo`, `maximo` |
| `eventos.tipo` | Valor de `tipo_evento` a analizar (uno por corrida) |
| `eventos.duracion_minima` | Duración mínima (en pentadas) de los eventos incluidos |
| `valores_minimos` | Valor mínimo de detección por variable, para generar el ruido de las series perturbadas. Debe cubrir cada variable usada en `variables_copulas` |
| `umbral.p.valor` | Nivel de significancia de todos los tests |
| `n.series.perturbadas` | Cantidad de series con ruido por variable |
| `n.realizaciones` | Se usan las realizaciones sintéticas `<= n.realizaciones` |
| `s.series.perturbadas` | Semilla (`set.seed`) para reproducir la corrida |
| `min_cantidad_valores_ajuste_univariado` | (opcional, 30) Mínimo de eventos para ajustar una distribución |
| `variables_sin_test_continuidad` | (opcional) Variables discretas (muchos empates) a las que se omiten KS/AD/CvM |
| `tests.estacionariedad`, `tests.dependencia` | Tests a aplicar (cada uno con su función `Test<nombre>` en `lib/`) |
| `configuracion.ajuste.copula`, `configuracion.ajuste.univariado` | Tablas familia/distribución → función de ajuste + flag `uso_sequias`. Vienen completas; normalmente no se editan |
| `periodo_retorno.niveles_anios` | Isolíneas/líneas de referencia (ej. 2, 5, 10, … 500 años) |
| `periodo_retorno.resolucion_grilla`, `margen_grilla` | Puntos de la grilla y cuánto extenderla más allá del máximo observado (fracción del rango) |
| `periodo_retorno.colorear_enso` | (opcional, `true`) Colorear por ENSO los puntos del gráfico bivariado |
| `periodo_retorno.graficar_distribucion_ajuste` | (opcional, `false`) Histograma + densidad por ubicación y variable (muchos PNG) |

### `data/configuracion_archivos_utilizados.yml`

| Clave | Qué es |
|---|---|
| `identificador_corrida` | Sufijo de todos los archivos generados (`<*idc>` en las rutas). Cambiarlo evita pisar otra corrida |
| `eventos_identificados` | CSV de eventos, relativo a `dir.data` |
| `eventos_enso` | Archivo ENSO (`input/enso_roni.txt`) |
| `copulas.*` | Un `.rds` por artefacto intermedio/final |

## 5. Problema conocido: test `Sn` con variables casi discretas

`TestSn` (`lib/funciones_bondad_ajuste_copulas.R`) usa
`copula::gofCopula(method = "Sn", simulation = "pb")`, que reestima la cópula
en cada réplica con `optim(method = "BFGS")`. Si una de las dos variables tiene
muy pocos valores distintos respecto de la cantidad de eventos (ej. una
duración entera con 7 valores sobre 82 eventos), la log-verosimilitud queda
degenerada y `optim()` falla con `non-finite finite-difference value`. Como
`TestearBondadAjusteCopulas` captura ese fallo pero la línea siguiente
(`dplyr::filter(parametro == "p.value")`) no tolera un resultado sin columnas,
el paso 10 aborta todo el script con `objeto 'parametro' no encontrado`.

**Cómo evitarlo:** elegir pares de variables razonablemente continuas; evitar
una variable muy discretizada como uno de los dos lados de la cópula.

Para ver el error real de un worker (que `doSNOW` oculta), reproducir la
llamada fuera del cluster, en una sesión de R común.
