# Cópulas — Análisis multivariado de eventos extremos (CRC-SAS)

## Overview

Componente del sistema de monitoreo de sequías y eventos extremos de **CRC-SAS**.
A partir de eventos ya identificados (sequías, olas de calor, olas de frío,
períodos cálidos/fríos…) con sus métricas de **duración, intensidad, magnitud y
valor mínimo/máximo**, el pipeline estima **cuán raro es un evento** teniendo en
cuenta varias de esas variables a la vez:

1. Ajusta, para cada variable, la distribución que mejor la describe.
2. Verifica estacionariedad e independencia de las series.
3. Ajusta **cópulas** (Gumbel, Frank, Joe, Clayton, Normal, t) entre pares de
   variables (ej. duración-intensidad) y elige la mejor familia.
4. Combina cópula y distribuciones en una distribución multivariada y calcula
   **períodos de retorno combinados** (ej. "una sequía de 4 meses *y* de
   intensidad extrema ocurre una vez cada N años"), con gráficos de isolíneas.

Opcionalmente, los eventos se clasifican según la fase ENSO (Niño/Niña/Neutro)
vigente al inicio, y los gráficos las distinguen por color.

> La identificación de eventos y el generador estocástico de series sintéticas
> son procesos previos, **fuera de este repositorio**. Acá se consume su salida
> (un CSV de eventos). Contexto conceptual: [`docs/guia_video.md`](docs/guia_video.md).

## Arquitectura

Un único script, [`01_copulas.R`](01_copulas.R), ejecuta el pipeline completo
(14 pasos) en paralelo para todas las ubicaciones, variables y cópulas
configuradas. Se dirige con 3 archivos YAML (rutas y procesos, parámetros del
análisis, nombres de archivos) y se apoya en un pequeño framework de ejecución
(`lib/R/`: logging y tareas paralelas) y en funciones de negocio (`lib/`).

```
CSV de eventos ─► ajuste univariado ─► estacionariedad / dependencia ─► ajuste de cópulas
                                                                         │
            gráficos y períodos de retorno ◄─ distribución multivariada ◄┘
```

Detalle de pasos, código y parámetros: [`docs/arquitectura.md`](docs/arquitectura.md).

## Instalación y configuración

Necesitás **R** (4.x) y **Git**.

```bash
git clone https://github.com/CRC-SAS/copulas.git
cd copulas
cp configuracion_copulas.yml.tmpl configuracion_copulas.yml   # y editar las rutas
```

1. Instalar los paquetes de R (lista en la guía).
2. En `configuracion_copulas.yml`, poner la ruta absoluta del proyecto y la
   cantidad de procesos.
3. Copiar a `data/input/` el **CSV de eventos** y el archivo ENSO `enso_roni.txt`
   (no vienen en el repo).
4. En `parametros_copulas.yml`, elegir **estaciones**, **variables** y **tipo de
   evento**.
5. Correr:

```bash
Rscript 01_copulas.R
```

Los resultados quedan en `data/output/` y el avance en `run/CalcCopulas.log`.

**Guía completa, para quien nunca lo corrió:**
[`docs/guia_ejecucion.md`](docs/guia_ejecucion.md) — una o varias estaciones,
dónde indicar variables y tipos de evento, y uso opcional de ENSO.
