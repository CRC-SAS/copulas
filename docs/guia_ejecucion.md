# Guía de ejecución (paso a paso)

Esta guía es para alguien que se baja este repositorio por primera vez y
necesita hacerlo correr, sin conocer todavía los detalles internos. Si
después querés entender *cómo* funciona el pipeline por dentro (qué hace
cada paso, qué significa cada test estadístico), esa información está en el
[`README.md`](../README.md) principal.

Acá vemos **dos formas de correr el pipeline**:

1. **Un escenario simple** — una sola corrida, sobre las estaciones y
   variables que vos elijas.
2. **Múltiples escenarios (batch)** — encadenar varias corridas
   automáticamente, típicamente una por cada "tipo de evento" (ola de calor,
   ola de frío, período cálido, etc.).

Antes de eso, una aclaración corta que evita mucha confusión más adelante.

## ¿Qué es "un escenario"?

Una corrida del pipeline (`01_copulas.R`) procesa, **todas juntas, en la
misma ejecución**:

- Todas las **estaciones** que le indiques (pueden ser 1 o pueden ser 100).
- Todas las **combinaciones de variables** que le indiques (ej.
  intensidad-duración, intensidad-magnitud).

Es decir: **procesar muchas estaciones, o probar varios pares de
variables a la vez, ya es "gratis" dentro de un único escenario** — no hace
falta ningún mecanismo de batch para eso, alcanza con listarlas en la
configuración (sección siguiente).

Entonces, ¿para qué existe el modo batch? Porque hay una cosa que **sí**
obliga a hacer una corrida aparte por cada valor: el **tipo de evento**
(ej. "ola de calor" vs. "período frío"). Cada tipo de evento:

- tiene su propio identificador de corrida (para no pisar los resultados de
  otro tipo),
- y no todas las estaciones tienen datos para todos los tipos (una estación
  puede tener olas de calor identificadas pero no períodos fríos).

El modo batch (caso 2 de esta guía) automatiza exactamente eso: correr el
mismo pipeline una vez por cada tipo de evento, con la lista de estaciones
correcta para cada uno.

## Qué necesitás antes de arrancar

**1) R y los paquetes del pipeline.** La lista completa está en el
`README.md` (sección "Dependencias de R"). Instalación rápida desde R:

```r
install.packages(c("dplyr", "purrr", "lubridate", "magrittr", "lmomco",
  "stringr", "yaml", "goftest", "WRS2", "futile.logger", "doSNOW", "foreach",
  "iterators", "snow", "yardstick", "hydroGOF", "copula", "ggplot2",
  "ggnewscale", "R6", "RPostgres", "data.table", "glue", "tidyr", "tibble",
  "rlang", "xts", "npcp", "Kendall", "caret"))
```

> Si `install.packages("copula")` falla, probablemente falte la librería de
> sistema `libgsl-dev` (Linux: `sudo apt-get install libgsl-dev`) — `copula`
> depende del paquete de R `gsl`, que a su vez necesita esa librería
> compilada en el sistema.

**2) Un archivo de eventos ya identificados.** Este pipeline **no**
identifica eventos secos por sí solo — parte de un CSV con eventos ya
calculados por fuera del repo (duración, intensidad, magnitud, mínimo,
máximo de cada evento, por estación). Sin este archivo no hay nada para
analizar. El formato exacto está documentado en el `README.md`, sección 3.

**3) El archivo de referencia ENSO** (`data/input/enso_roni.txt`). El
pipeline enriquece cada evento con el fenómeno ENSO (Niño/Niña/Neutro)
vigente al momento de su inicio, y **esto es obligatorio**: si este archivo
no está, el paso 4 del pipeline falla para cualquier corrida (simple o
batch). Es un archivo de datos, no de código — no viene versionado en el
repo (igual que el CSV de eventos), así que hay que copiarlo a mano a esa
ruta antes de la primera corrida.

## Caso 1: Un escenario simple

Pensalo como "quiero correr el análisis para tal estación (o tales
estaciones), con tales pares de variables".

### Paso 1 — Configurar el escenario

Editá `parametros_copulas.yml` (ya viene un ejemplo armado, hay que
ajustarlo a tu caso):

```yaml
ubicaciones:
  - { id: "87548", nombre: "Junín" }
  # agregá una línea por cada estación que quieras incluir

variables_copulas:
  - { variable_x: "duracion", variable_y: "intensidad" }
  # agregá un par por cada combinación de variables que quieras probar

eventos:
  tipo: "seco"              # debe coincidir con la columna tipo_evento del CSV
  duracion_minima: 6
```

(Hay más parámetros en ese archivo — umbrales estadísticos, cantidad de
series con ruido, etc. — el README los explica todos en detalle. Para
arrancar, con tocar `ubicaciones`, `variables_copulas` y `eventos.tipo`
alcanza.)

### Paso 2 — Decirle al pipeline dónde está tu CSV de eventos

En `data/configuracion_archivos_utilizados.yml`:

```yaml
eventos_identificados: "input/mi_archivo_de_eventos.csv"
```

Ese archivo tiene que existir en `data/input/`.

### Paso 3 — Correr

```bash
Rscript 01_copulas.R
```

(Usa por defecto `configuracion_copulas.yml`, `parametros_copulas.yml` y
`data/configuracion_archivos_utilizados.yml`. Ver README sección 4 si
querés pasar rutas distintas por línea de comandos.)

### Paso 4 — Ver los resultados

- **Mientras corre:** el progreso se va escribiendo en `run/CalcCopulas.log`.
- **Al terminar**, en `data/output/`:
  - `copulas_<identificador>.rds` — el resultado final (cópula + ajustes
    elegidos, por estación y par de variables).
  - `periodo_retorno_<identificador>.rds` + un `.png` por estación y par de
    variables — el gráfico de período de retorno combinado (el que muestra,
    para cada combinación de valores, cada cuántos años se espera un evento
    así de extremo).
  - `periodo_retorno_univariado_<identificador>.rds` + PNGs — lo mismo pero
    mirando cada variable por separado.

`<identificador>` es el valor de `identificador_corrida` en
`data/configuracion_archivos_utilizados.yml` (por defecto `id1`).

**¿Cuánto tarda?** Para pocas estaciones (1 a 5) y pocos pares de
variables, minutos. El tiempo crece con la cantidad de estaciones × pares de
variables — para tener una referencia, corridas de producción con ~70
estaciones tardaron entre 1 y 4 horas cada una (ver caso 2).

## Caso 2: Múltiples escenarios (batch, por tipo de evento)

Usalo cuando tenés un CSV con **varios tipos de evento** (por ejemplo,
`OlaCalorTXyTN`, `OlaFrioTXyTN`, `PeriodoCalidoTN`, `PeriodoCalidoTX`,
`PeriodoFrioTN`, `PeriodoFrioTX`) y querés correr el pipeline para todos,
sin editar el YAML a mano cada vez.

Esta tooling vive en `pruebas/` y **no está versionada** (es infraestructura
de trabajo, no parte del pipeline en sí) — si no la encontrás en tu
checkout, hay que recrearla o pedírsela a alguien del equipo.

### Paso 1 — Preparar el archivo de configuración "por tipo"

A diferencia del escenario simple, acá se usa una variante de
`data/configuracion_archivos_utilizados.yml` pensada para elegir el tipo de
evento en cada corrida: `data/configuracion_archivos_utilizados_frio_calor.yml`.
Tampoco está versionada — si no existe, se crea copiando la base y
cambiando una sola línea:

```bash
cp data/configuracion_archivos_utilizados.yml data/configuracion_archivos_utilizados_frio_calor.yml
```

y en la copia, cambiar:

```yaml
eventos_identificados: "input/eventos_identificados_frio_calor.csv"
```

(apuntando al CSV que tiene los 6 tipos de evento juntos, con la columna
`tipo_evento` distinguiéndolos).

### Paso 2 — Definir los pares de variables una sola vez

`pruebas/parametros_copulas_test.yml` (sí está versionado) es la base de
parámetros que usa el modo batch. Los pares de variables
(`variables_copulas`) y el resto de los umbrales se definen ahí **una sola
vez** — valen para los 6 tipos por igual, porque (como se explicó arriba)
elegir variables no requiere una corrida por separado.

### Paso 3 — Correr un tipo de evento

```bash
bash pruebas/correr_un_tipo_evento.sh OlaCalorTXyTN
```

Esto, automáticamente:

1. Genera la lista de estaciones correcta para ese tipo (filtra el CSV:
   no todas las estaciones tienen todos los tipos de evento).
2. Le pone a esa corrida un identificador propio (ej. `olacalor`), para que
   no pise los resultados de otro tipo.
3. Corre `01_copulas.R` con esa configuración.
4. Mueve los gráficos generados a su propia carpeta:
   `data/output/png_olacalor/`.

Tipos válidos: `OlaCalorTXyTN`, `OlaFrioTXyTN`, `PeriodoCalidoTN`,
`PeriodoCalidoTX`, `PeriodoFrioTN`, `PeriodoFrioTX` (tienen que coincidir
con los valores de la columna `tipo_evento` del CSV).

**¿Cuánto tarda?** Con ~70 estaciones, entre 1 y 4 horas por tipo,
dependiendo de cuántas estaciones tiene ese tipo específico (algunos tipos
tienen menos estaciones con datos que otros).

### Paso 4 — Correr los 6 tipos

**Para una corrida de producción real (muchas estaciones), correr los 6
tipos de a uno**, repitiendo el paso 3 con cada nombre. No conviene
encadenarlos automáticamente: cada uno puede tardar horas, y encadenar los
6 sin pausas puede saturar la máquina por casi un día entero sin que nadie
esté controlando que vaya bien.

```bash
bash pruebas/correr_un_tipo_evento.sh OlaCalorTXyTN
# esperar a que termine, revisar el log, recién ahí seguir con el próximo
bash pruebas/correr_un_tipo_evento.sh OlaFrioTXyTN
# ... y así con los 6
```

Para corridas **chicas** (pocas estaciones, de prueba/diagnóstico), sí existe
un atajo que encadena los 6 automáticamente:

```bash
bash pruebas/correr_18_escenarios.sh
```

(el nombre viene de 6 tipos × 3 pares de variables = 18 combinaciones — pero
ojo, como se explicó arriba, eso no son 18 corridas: son 6 corridas, cada
una resolviendo sus 3 pares de variables de una sola vez).

### Paso 5 — Ver los resultados

Cada tipo deja sus resultados con su propio identificador, todos dentro de
`data/output/`:

- `copulas_<tipo>.rds`, `periodo_retorno_<tipo>.rds`,
  `periodo_retorno_univariado_<tipo>.rds`
- `png_<tipo>/` — todos los gráficos de ese tipo (bivariados + univariados)

(`<tipo>` acá es la versión corta que usa `correr_un_tipo_evento.sh`, ej.
`olacalor`, `periodocalidotn` — no el nombre largo `OlaCalorTXyTN`.)

### Cómo agregar/cambiar algo

- **¿Querés agregar una estación?** No hace falta tocar nada del modo
  batch: la lista de estaciones por tipo se genera sola a partir del CSV de
  eventos (`pruebas/generar_parametros_tipo.py`, filtra por `tipo_evento`).
  Si la estación tiene eventos de ese tipo en el CSV, ya va a aparecer.
- **¿Querés agregar/cambiar un par de variables?** Editá
  `variables_copulas` en `pruebas/parametros_copulas_test.yml` (una sola
  vez, vale para los 6 tipos). Si agregás una variable nueva que no esté ya
  en `valores_minimos` en ese mismo archivo, agregala ahí también.
- **¿Querés agregar un tipo de evento nuevo?** Tiene que existir esa
  combinación en la columna `tipo_evento` del CSV. Después alcanza con
  correr `bash pruebas/correr_un_tipo_evento.sh <TipoNuevo>` — la tooling no
  tiene una lista fija de tipos "permitidos", lee lo que haya en el CSV.

## Problemas frecuentes

**"El archivo ... no existe"** al arrancar — alguno de los 3 YAML, el CSV de
eventos, o `data/input/enso_roni.txt` no está en la ruta esperada. Revisá
la sección "Qué necesitás antes de arrancar" de esta guía.

**El script aborta a mitad de camino con un error de R** — buscá primero en
la sección "Problemas conocidos" del `README.md` principal, ahí está
documentado el problema más común (variables con muy pocos valores
distintos rompiendo el test `Sn`).

**Quiero ver qué está pasando mientras corre** — `run/CalcCopulas.log` tiene
el resumen; si un paso específico falla, cada paso paralelo deja su propio
log más detallado en `run/CalcCopulas-<NombreDelPaso>.log`.
