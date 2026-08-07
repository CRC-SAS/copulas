# Cópulas — Análisis multivariado de eventos secos (CRC-SAS)

## 1. Overview

Este repositorio es un componente del sistema de monitoreo de sequías de **CRC-SAS**
(Comité Regional de Recursos Hídricos / Sistema de Alerta y Sequías). A partir de
eventos secos ya identificados (con sus métricas de duración, intensidad, magnitud
y valor mínimo) en distintas ubicaciones, el pipeline:

1. Ajusta, para cada variable, la distribución univariada que mejor la describe
   (L-momentos y máxima verosimilitud sobre ~20 familias de distribución vía
   [`lmomco`](https://cran.r-project.org/package=lmomco)).
2. Evalúa **estacionariedad** e **independencia** de las series involucradas.
3. Ajusta **cópulas** bivariadas (Gumbel, Frank, Joe, Clayton, Normal, t-Student)
   entre pares de variables (p. ej. duración-intensidad).
4. Selecciona la familia de cópula que mejor ajusta (AIC, BIC, MRE, RMSE,
   validación cruzada, test Sn).
5. Combina la mejor cópula con las mejores distribuciones marginales en una
   distribución multivariada final (`copula::mvdc`), que permite estimar
   **periodos de retorno combinados** (ej. "sequía de 4 meses **y** con
   intensidad extrema").

Para el contexto conceptual completo (SPI/SPEI, pentadas, generador estocástico
de series sintéticas, cópulas) ver [`docs/guia_video.md`](docs/guia_video.md).
El manual metodológico extendido está en
[`docs/manual_copulas.Rmd`](docs/manual_copulas.Rmd)
(renderizado en `docs/manual_copulas.html` / `.pdf`).

> **Fuera de alcance de este repo:** la identificación de eventos secos y el
> generador estocástico de series sintéticas son procesos previos, corridos por
> fuera de este repositorio. Acá se **consume** su salida (ver sección 3).

## 2. Arquitectura

### 2.1 Scripts principales

| Script | Rol | Estado |
|---|---|---|
| `01_copulas.R` | Pipeline de producción. Dirigido por 3 archivos YAML, ejecuta el proceso completo (12 pasos) en paralelo para todas las ubicaciones/variables/cópulas configuradas. | Activo — último fix funcional en 2021-03. |
| `01_script_maestro.R` | Prototipo exploratorio para una sola estación (hardcodea `estacion.usar`, carga un `.rda` desde una ruta local de otro colaborador). Sirvió para iterar antes de construir `01_copulas.R`. | Desactualizado — ver sección "Archivos en desuso" más abajo. |

### 2.2 Flujo del pipeline (`01_copulas.R`)

1. Cargar paquetes y leer los 3 YAML de configuración.
2. Leer los eventos identificados y generar las series "perturbadas" (con ruido,
   según `n.series.perturbadas`) usadas para propagar incertidumbre.
3. Ajuste univariado por (ubicación, variable, distribución) — `AjusteUnivariadoUVD`.
4. Determinar el mejor ajuste univariado por (ubicación, variable) — `MejorAjusteUnivariadoUV`.
5. Aplicar el mejor ajuste a cada serie perturbada — `AplicarMejorAjusteASeriesPerturbadas`.
6. Determinar estacionariedad — `DeterminarEstacionariedad` (Box-Pierce/Ljung-Box,
   punto de cambio uni/autocópula, Mann-Kendall, test empírico).
7. Determinar dependencia entre pares de variables — `DeterminarDependencia`
   (Kendall, Spearman, cópula empírica).
8. Ajustar cópulas por familia — `AjustarCopulas` (solo si la serie es
   estacionaria y las variables son dependientes).
9. Determinar la mejor familia de cópula por (ubicación, par de variables) —
   `MejorAjusteMultivariadoUC`.
10. Aplicar el mejor ajuste multivariado final — `AplicarMejorAjusteACopulas`
    (combina cópula + márgenes en un objeto `copula::mvdc`).

Cada paso se ejecuta mediante `Task$run()`, que paraleliza sobre las filas de
`input.values` usando `doSNOW`/`foreach`/`snow`, respetando `max.procesos`
(definido en YAML). Los errores de cada tarea se acumulan y, si hay alguno, el
script aborta (`script$error(...)`) después de loguearlos.

### 2.3 Framework de ejecución (`lib/R/`)

- `Script.R` — clase R6 que representa la corrida completa; maneja logging
  (`futile.logger`) hacia `run/CalcCopulas.log`, con `start()` / `stop()` /
  `info()` / `warn()` / `error()`.
- `Task.R` — clase R6 que ejecuta una función "worker" en paralelo sobre un
  conjunto de inputs, con su propio log/out por tarea
  (`run/CalcCopulas-<función>.log` / `.out`).
- `Helpers.R` — utilidades chicas (`IdentificarIdColumn`, que detecta si la
  ubicación usa `station_id`, `point_id`, o cualquier columna `*_id`).

### 2.4 Funciones de negocio (`lib/`)

| Archivo | Contenido |
|---|---|
| `funciones_worker.R` | Wrapper de cada paso del pipeline (uno por función paralelizada, ver 2.2) |
| `funciones_ajustes_marginales.R` | ~25 funciones de ajuste univariado por L-momentos y máxima verosimilitud (Gamma, GEV, Log-normal, Weibull, etc., vía `lmomco`) |
| `funciones_ajuste_distribuciones.R`, `ajuste_distribuciones.R` | Orquestación del ajuste univariado (`AjusteUnivariado*`) |
| `funciones_bondad_ajuste.R` | Tests de bondad de ajuste univariado (KS, Anderson-Darling, Cramér-von Mises, RMSE/IQR) |
| `funciones_mejor_ajuste_univariado.R` | Selección de la mejor distribución/método (por RMSE, CCC y test de cuantiles) |
| `funciones_test_estacionaridad.R` | Tests de estacionariedad (Box-Pierce/Ljung-Box, punto de cambio uni/multivariado, Mann-Kendall, empírico) |
| `funciones_test_independencia.R` | Tests de dependencia bivariada (Kendall, Spearman, cópula empírica) |
| `funciones_ajuste_familias_copulas.R` | Ajuste de cada familia de cópula (Gumbel, Frank, AMH, Joe, Clayton, Normal, t) |
| `funciones_bondad_ajuste_copulas.R` | Bondad de ajuste multivariado (Sn, SnB, SnC, AIC, BIC, MRE/RMSE, validación cruzada) |
| `funciones_mejor_ajuste_copula.R` | Selección de la mejor familia de cópula |
| `funciones_auxiliares.R` | Perturbación de series (`AgregarRuido`) y gráficos exploratorios (Chi-plot, K-plot) |

### 2.5 Datos y configuración

```
configuracion_copulas.yml(.tmpl)            → paths de trabajo y cantidad de procesos paralelos
parametros_copulas.yml                      → parámetros del análisis (ubicaciones, variables,
                                               tests, umbrales, distribuciones/cópulas habilitadas)
data/configuracion_archivos_utilizados.yml  → nombres de los archivos de entrada/salida de cada etapa
data/input/                                 → insumos (CSV de eventos, ver sección 3) — no versionado
data/partial/                               → resultados intermedios (semilla, eventos completos,
                                               estacionariedad, dependencia) — no versionado
data/output/                                → resultados finales (ajustes univariados/multivariados,
                                               cópulas) — no versionado
run/                                        → logs y pid de la corrida — no versionado
docs/                                       → guía conceptual y manual metodológico
```

## 3. Archivos necesarios para correr

Se necesitan **3 YAML + 1 CSV**:

**1) `configuracion_copulas.yml`** — ya versionado, define paths absolutos y
cantidad de procesos paralelos:

```yaml
dir:
  base: /devel/CRC-SAS/copulas/
  run: /devel/CRC-SAS/copulas/run/
  lib: /devel/CRC-SAS/copulas/lib/R/
  data: /devel/CRC-SAS/copulas/data/
max.procesos: 8
```

**2) `parametros_copulas.yml`** — ya versionado con valores de ejemplo/prueba;
ajustar a los datos reales:

```yaml
ubicaciones:
  - { id: "87548", nombre: "Junín" }
configuraciones.eventos:                 # solo se soporta UNA fila
  - [ conf_id, indice, escala, distribucion, metodo_ajuste ]
  - [ 3, "SPI", 6, "Gamma", "NoParametrico" ]
variables_copulas:
  - { variable_x: "duracion", variable_y: "minimo" }
  - { variable_x: "duracion", variable_y: "intensidad" }
eventos:
  tipo: "seco"
  duracion_minima: 6          # expresado en pentadas
valores_minimos:
  - { variable: "duracion", valor_minimo_deteccion: 1 }
  - { variable: "intensidad", valor_minimo_deteccion: 0.1 }
n.series.perturbadas: 2
n.realizaciones: 2
s.series.perturbadas: 123
```

**3) `data/configuracion_archivos_utilizados.yml`** — ya versionado, define
nombres de archivos de entrada/salida:

```yaml
identificador_corrida: &idc "id1"
eventos_identificados: "input/eventos_identificados_unconditional.csv"
copulas:
  resultado_final: "output/copulas_<*idc>.rds"
  # ... (ver archivo completo para el resto de los artefactos intermedios)
```

**4) CSV de eventos** (`data/input/eventos_identificados_unconditional.csv` por
defecto, **no versionado**) — una fila por evento, columnas obligatorias:

| Columna | Tipo | Descripción |
|---|---|---|
| `station_id` (o `point_id` / cualquier `*_id`) | texto/num | id de la ubicación |
| `tipo_evento` | texto | debe coincidir con `parametros_copulas.yml: eventos.tipo` |
| `conf_id` | entero | debe coincidir con `configuraciones.eventos` |
| `realizacion` | entero | número de realización sintética (≤ `n.realizaciones`) |
| `numero_evento` | entero | correlativo del evento |
| `fecha_inicio`, `fecha_fin` | fecha (`YYYY-MM-DD`) | límites del evento |
| `referencia_comienzo`, `referencia_fin` | fecha (`YYYY-MM-DD`) | período de referencia |
| `duracion`, `intensidad`, `magnitud`, `minimo`, `maximo` | numérico | métricas del evento (se toman en valor absoluto) |

Ejemplo mínimo:

```csv
station_id,tipo_evento,conf_id,realizacion,numero_evento,fecha_inicio,fecha_fin,referencia_comienzo,referencia_fin,duracion,intensidad,magnitud,minimo,maximo
87548,seco,3,1,1,2001-01-05,2001-03-10,2001-01-01,2001-03-31,13,1.2,15.6,-2.1,-0.8
87548,seco,3,1,2,2005-06-15,2005-08-20,2005-06-01,2005-08-31,10,0.9,9.0,-1.7,-0.5
```

## 4. Ejemplo de ejecución

```bash
# Usando los 3 YAML por defecto: busca configuracion_copulas.yml y
# parametros_copulas.yml en el working directory, y
# configuracion_archivos_utilizados.yml en {dir.data}
Rscript 01_copulas.R

# Indicando explícitamente cada YAML
Rscript 01_copulas.R configuracion_copulas.yml parametros_copulas.yml data/configuracion_archivos_utilizados.yml
```

Salidas de una corrida:

- Log de la corrida: `run/CalcCopulas.log` (+ un `.log`/`.out` por cada tarea paralela)
- Resultados intermedios: `data/partial/*_id1.rds`
- Resultados finales: `data/output/*_id1.rds`
  (`mejores_ajustes_univariados_id1.rds`, `mejores_ajustes_multivariados_id1.rds`,
  `copulas_id1.rds`)
- Info de la corrida: `data/partial/copulas_id1.info`

(`id1` es el `identificador_corrida` definido en
`data/configuracion_archivos_utilizados.yml`.)

## Dependencias de R

Ver `list.of.packages` en `01_copulas.R`: `dplyr`, `purrr`, `lubridate`,
`magrittr`, `lmomco`, `stringr`, `yaml`, `goftest`, `WRS2`, `futile.logger`,
`doSNOW`, `foreach`, `iterators`, `snow`, `yardstick`, `hydroGOF`, `copula`,
`ggplot2`, `R6`, `RPostgres`. Además se usan (vía `::`, sin `require`
explícito en el script): `data.table`, `glue`, `tidyr`, `tibble`, `rlang`,
`xts`, `npcp`, `Kendall`, `caret`.

> **Nota:** las funciones `TestBoxPierceLjungBox`/`TestEstacionaridadEmpirico`
> y los tests de independencia (`lib/funciones_test_estacionaridad.R`,
> `lib/funciones_test_independencia.R`) llaman a una función
> `ParametrosADataFrame` que no está definida en ningún archivo de este
> repositorio. Antes de correr el pipeline conviene confirmar de dónde debe
> provenir (probablemente un paquete interno de CRC-SAS no incluido todavía).
