# Guía de ejecución

Guía paso a paso para instalar el pipeline de cópulas y correrlo sobre una o
varias estaciones, con las variables y el tipo de evento que elijas. No hace
falta saber R: alcanza con editar unos archivos de texto y copiar un comando.
Para entender cómo funciona por dentro, ver [`arquitectura.md`](arquitectura.md).

**Índice**

1. [Qué necesitás](#1-qué-necesitás)
2. [Instalación](#2-instalación)
3. [Configuración inicial (una sola vez)](#3-configuración-inicial-una-sola-vez)
4. [El CSV de eventos](#4-el-csv-de-eventos)
5. [Correr el pipeline](#5-correr-el-pipeline)
6. [ENSO (opcional)](#6-enso-opcional)
7. [Ver los resultados](#7-ver-los-resultados)
8. [Problemas frecuentes](#8-problemas-frecuentes)

## 1. Qué necesitás

- **R** (4.x) y **Git**. En Linux, además, las librerías de sistema que usan
  algunos paquetes de R (`sudo apt-get install libgsl-dev` para `copula`).
- **El CSV de eventos identificados.** El pipeline *no* identifica eventos: parte
  de un CSV ya calculado por fuera (ver [formato](#4-el-csv-de-eventos)).
- **El archivo ENSO** `data/input/enso_roni.txt`. Hoy es obligatorio aunque no uses
  nada de ENSO (ver [sección 6](#6-enso-opcional)).

Ninguno de los dos archivos de datos viene en el repositorio: hay que copiarlos
a `data/input/`.

## 2. Instalación

```bash
git clone https://github.com/CRC-SAS/copulas.git
cd copulas
```

Desde una consola de R (o RStudio), instalar los paquetes una sola vez:

```r
install.packages(c("dplyr", "purrr", "lubridate", "magrittr", "lmomco",
  "stringr", "yaml", "goftest", "WRS2", "futile.logger", "doSNOW", "foreach",
  "iterators", "snow", "yardstick", "hydroGOF", "copula", "ggplot2",
  "ggnewscale", "R6", "RPostgres", "data.table", "glue", "tidyr", "tibble",
  "rlang", "xts", "npcp", "Kendall", "caret"))
```

## 3. Configuración inicial (una sola vez)

El pipeline lee **3 archivos de configuración**. Dos ya vienen listos y uno hay
que crearlo.

### 3.1 `configuracion_copulas.yml` — dónde está el proyecto (hay que crearlo)

Copiá la plantilla y editala:

```bash
cp configuracion_copulas.yml.tmpl configuracion_copulas.yml
```

Reemplazá cada `${base}/Copulas/` por la **ruta absoluta** de la carpeta donde
clonaste el repo, y ajustá `max.procesos` (cuántos núcleos del procesador usar):

```yaml
dir:
  base: /home/usuario/copulas/
  run:  /home/usuario/copulas/run/
  lib:  /home/usuario/copulas/lib/R/
  data: /home/usuario/copulas/data/
max.procesos: 4
```

Este archivo es propio de cada computadora y no se sube al repositorio.

### 3.2 `data/configuracion_archivos_utilizados.yml` — qué archivos usa

Ya viene armado. Solo hay tres líneas que te interesan:

```yaml
identificador_corrida: &idc "id1"                                  # (a)
eventos_identificados: "input/eventos_identificados_unconditional.csv"  # (b)
eventos_enso: "input/enso_roni.txt"                                # (c)
```

- **(a)** es un nombre que se agrega a todos los resultados (`copulas_id1.rds`,
  etc.). Cambialo en cada corrida si no querés pisar los resultados anteriores.
- **(b)** es el nombre de tu CSV de eventos, dentro de `data/`.
- **(c)** es el archivo ENSO, dentro de `data/`.

### 3.3 `parametros_copulas.yml` — qué se analiza

Es el archivo que vas a editar en cada corrida. Ver la sección siguiente.

## 4. El CSV de eventos

Una fila por evento (observado o de una realización sintética). Columnas
obligatorias:

| Columna | Descripción |
|---|---|
| `station_id` (o `point_id`, o cualquier `*_id`) | Id de la ubicación |
| `tipo_evento` | Tipo de evento (ej. `OlaCalorTXyTN`, `seco`) |
| `conf_id` | Entero; debe coincidir con `configuraciones.eventos` |
| `realizacion` | Número de realización sintética (se usan las `<= n.realizaciones`) |
| `numero_evento` | Correlativo del evento |
| `fecha_inicio`, `fecha_fin` | `YYYY-MM-DD` |
| `referencia_comienzo`, `referencia_fin` | `YYYY-MM-DD` |
| `duracion`, `intensidad`, `magnitud`, `minimo`, `maximo` | Métricas del evento (se toman en valor absoluto) |

```csv
station_id,tipo_evento,conf_id,realizacion,numero_evento,fecha_inicio,fecha_fin,referencia_comienzo,referencia_fin,duracion,intensidad,magnitud,minimo,maximo
87548,seco,3,1,1,2001-01-05,2001-03-10,2001-01-01,2001-03-31,13,1.2,15.6,-2.1,-0.8
```

Opcional: una columna `categoria_episodio` (ver [ENSO](#6-enso-opcional)).

## 5. Correr el pipeline

Editá `parametros_copulas.yml`. Para una corrida común se tocan cuatro cosas.
Todo lo demás (umbrales, tests, familias de cópulas) puede quedar como viene.

### 5.1 Estaciones: una o muchas

En `ubicaciones`, una línea por estación. El `id` tiene que existir en el CSV
**para el tipo de evento elegido**:

```yaml
ubicaciones:
  - { id: "87548", nombre: "Junín" }
  - { id: "87624", nombre: "Anguil INTA" }
```

Procesar una o cien estaciones es la misma corrida: alcanza con listarlas. Si
alguna no tiene eventos de ese tipo en el CSV, el script aborta con "No hay datos,
en eventos, para todas las ubicaciones".

### 5.2 Variables

En `variables_copulas`, un par por línea. Cada par genera su propia cópula. Las
variables posibles son `duracion`, `intensidad`, `magnitud`, `minimo` y `maximo`:

```yaml
variables_copulas:
  - { variable_x: "duracion", variable_y: "intensidad" }
  - { variable_x: "magnitud", variable_y: "intensidad" }
```

Cada variable usada tiene que figurar también en `valores_minimos` (el archivo ya
trae `duracion`, `intensidad`, `magnitud` y `minimo`; si usás otra, agregala).
Evitá variables con muy pocos valores distintos (ver
[problemas frecuentes](#8-problemas-frecuentes)).

### 5.3 Tipo de evento

En `eventos.tipo`, **un** tipo por corrida; tiene que coincidir con la columna
`tipo_evento` del CSV. `duracion_minima` descarta los eventos más cortos:

```yaml
eventos:
  tipo: "OlaCalorTXyTN"
  duracion_minima: 6
```

### 5.4 Ejecutar

Desde la carpeta del repo:

```bash
Rscript 01_copulas.R
```

Sin argumentos usa `configuracion_copulas.yml` y `parametros_copulas.yml` de la
carpeta actual y `data/configuracion_archivos_utilizados.yml`. El progreso se ve
en `run/CalcCopulas.log`. Con ~1-5 estaciones y pocos pares tarda minutos; con
~70 estaciones, entre 1 y 4 horas.

### 5.5 Varios tipos de evento

Cada tipo necesita su propia corrida (y su propio `identificador_corrida`, para no
pisar resultados). Para cada tipo:

1. Copiá los archivos con otro nombre:
   ```bash
   cp parametros_copulas.yml parametros_olafrio.yml
   cp data/configuracion_archivos_utilizados.yml data/archivos_olafrio.yml
   ```
2. En `parametros_olafrio.yml` cambiá `eventos.tipo` (y `ubicaciones`, porque no
   todas las estaciones tienen todos los tipos).
3. En `data/archivos_olafrio.yml` cambiá `identificador_corrida` (ej. `"olafrio"`).
4. Corré pasando **los tres** archivos, en este orden:
   ```bash
   Rscript 01_copulas.R configuracion_copulas.yml parametros_olafrio.yml data/archivos_olafrio.yml
   ```
5. Los `.rds` ya llevan el identificador en el nombre, pero los **PNG no**
   (`data/output/periodo_retorno_*.png`): movelos a otra carpeta antes de la
   próxima corrida, o se pisan.

Corré los tipos de a uno y revisá el log antes de lanzar el siguiente.

## 6. ENSO (opcional)

El pipeline relaciona los eventos con El Niño / La Niña a partir de
`data/input/enso_roni.txt` (índice RONI por trimestre móvil, tabulado, columnas
`Year, Season, Start_Month, End_Month, RONI, Tipo_evento`). Hay tres niveles de
uso; solo el primero es automático.

**a) Clasificación ENSO de cada evento (automática).** Cada evento recibe una
columna `evento_enso` (Niño/Niña/Neutro, con su intensidad) según el mes en que
empezó. Es lo que usan los gráficos. No requiere configuración, pero **el archivo
`enso_roni.txt` tiene que existir** aunque no te interese ENSO: hoy el script falla
al leerlo si falta.

**b) Gráficos con o sin colores ENSO.** Por defecto los puntos observados del
gráfico bivariado se pintan según el ENSO vigente. Para un gráfico sin colores
agregá a `parametros_copulas.yml`:

```yaml
periodo_retorno:
  colorear_enso: false
```

**c) Cópulas sobre los episodios ENSO en sí.** Si lo que querés analizar son los
propios episodios Niño/Niña (duración e intensidad de cada uno), armá un CSV de
eventos en el que cada episodio es una fila, con `station_id = ENSO`,
`tipo_evento = Nina` (o `Nino`) y una columna `categoria_episodio` (ej.
`Niña Moderada`), que se usa para colorear los puntos. Luego, en
`parametros_copulas.yml`:

```yaml
ubicaciones:
  - { id: "ENSO", nombre: "ENSO" }
eventos:
  tipo: "Nina"          # una corrida por Nina y otra por Nino
  duracion_minima: 1
variables_copulas:
  - { variable_x: "duracion", variable_y: "intensidad" }
min_cantidad_valores_ajuste_univariado: 20   # son pocos episodios (~25); el mínimo por defecto es 30
```

Y apuntá `eventos_identificados` al CSV de episodios en
`data/configuracion_archivos_utilizados.yml` (usando una copia, como en 5.5).

## 7. Ver los resultados

Todo queda en `data/output/` (`<id>` = `identificador_corrida`):

- `copulas_<id>.rds` — resultado final: cópula y distribuciones elegidas, por
  estación y par de variables.
- `periodo_retorno_<estación>_<x>_<y>.png` — gráfico de período de retorno
  combinado: para cada combinación de valores, cada cuántos años se espera un
  evento así de extremo.
- `periodo_retorno_univariado_<estación>_<variable>.png` — lo mismo para cada
  variable por separado.
- `periodo_retorno_<id>.rds` y `periodo_retorno_univariado_<id>.rds` — las tablas
  detrás de los gráficos.

Los `.rds` se abren desde R con `readRDS("data/output/copulas_id1.rds")`. La lista
completa de archivos está en [`arquitectura.md`](arquitectura.md#3-datos-y-archivos).

## 8. Problemas frecuentes

**"El archivo ... no existe"** al arrancar. Falta alguno de los 3 YAML, el CSV de
eventos o `data/input/enso_roni.txt`, o la ruta en `configuracion_copulas.yml` está
mal. Revisá que las rutas sean absolutas y terminen en la carpeta correcta.

**"Paquete no encontrado: …"**. Falta instalar ese paquete (sección 2).

**"No hay datos, en eventos, para todas las ubicaciones/variables…"**. Alguna
estación de `ubicaciones` no tiene eventos del `eventos.tipo` elegido (con la
`duracion_minima` y `n.realizaciones` dados), o alguna variable de
`variables_copulas` no está en el CSV o en `valores_minimos`.

**El script aborta a mitad de camino con `objeto 'parametro' no encontrado`.** Es el
test `Sn` fallando con una variable de pocos valores distintos (ej. duraciones
enteras con rango chico). Elegí otro par de variables; el detalle está en
[`arquitectura.md`](arquitectura.md#5-problema-conocido-test-sn-con-variables-casi-discretas).

**Una variable no se ajusta (resultado `NA`).** Tiene menos eventos que
`min_cantidad_valores_ajuste_univariado` (30 por defecto). Bajalo en
`parametros_copulas.yml` si corresponde.

**Quiero ver qué está pasando.** `run/CalcCopulas.log` tiene el resumen; cada paso
paralelo deja además su log en `run/CalcCopulas-<NombreDelPaso>.log`.
