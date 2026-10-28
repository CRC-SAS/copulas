# Figuras y números de la presentación de cópulas

Scripts que generan las figuras y las cifras de la presentación "Cópulas para el análisis de olas de calor y de frío" (deck de slides, octubre de 2026). Cada script se corre desde la raíz del repo.

## Entradas

- Salidas de la corrida con `Olas_20.rds` en `data/output/` (`copulas_*.rds`, `ajustes_multivariados_*.rds`, `mejores_ajustes_multivariados_*.rds`).
- `data/input/eventos_identificados_frio_calor.csv` y `data/input/enso_roni.txt`.
- `data/input/regiones.csv`: regiones de las estaciones (archivo de Natalia Herrera; Latin-1, separado por `;`). Se puede indicar otra ruta con la variable `REGIONES_CSV`. Cuatro estaciones no están en ese archivo y se asignaron por coordenadas (87345 a CENTRO, 87420 a CUYO, 87569 a BUENOS AIRES, 87934 a PATAGONIA).

## Orden de ejecución

Todos se corren desde la raíz del repo (`Rscript scripts/presentacion/<script>.R`). Los de la columna "Requiere" van después de los que allí se listan.

| Script | Requiere | Qué genera |
|---|---|---|
| `resumen_ajustes.R` | corrida | `resumen_ajustes.csv`: cola de la cópula ganadora, región, n de eventos y tau por estación y tipo |
| `resumen_ajustes_filtrados.R` | `resumen_ajustes.R` | `id_filtrado.csv`: pares con cópula ganadora tras los filtros de dependencia y estacionariedad (Bonferroni); % de cola |
| `aic_por_familia.R` | corrida | `aic_olafrio.rds`: AIC por familia (mediana sobre series perturbadas) |
| `dependencia_tests.R` | `data/partial/dependencia_<tipo>.rds` | `dependencia_tests.csv/.txt`: pares que pasan Kendall, Spearman y cópula empírica |
| `figuras_basicas.R` | corrida | Scatter de 87585, de los datos a la cópula, familias y colas, cópula ajustada, scatter por fase ENSO |
| `figura_evento.R` | corrida | Isolíneas de 87585 con la ola de diciembre de 2013 |
| `figura_tipos_copulas.R` | - | Los cuatro tipos de cópulas por cola (simulaciones) |
| `figuras_parsimonia_limites.R` | `aic_por_familia.R` | Parsimonioso contra sobreajustado, AIC con n chico, dependencia débil, ajustes sin ganadora |
| `figuras_regiones_n.R` | `resumen_ajustes.R` | Cola por grupo, mapa de estaciones, cola según n de eventos |
| `figuras_caso_87585.R` | corrida | Caso Buenos Aires 87585 (períodos cálidos y fríos de TX): isolíneas rotuladas, período de retorno, `caso_87585.txt` con las cifras citadas |
| `figuras_limites.R` | `figuras_caso_87585.R`, `resumen_ajustes_filtrados.R` | AIC con n chico y período de retorno con cópula contra independencia |
| `figuras_filtros.R` | `resumen_ajustes_filtrados.R`, `dependencia_tests.R`, `est_todos_tests.rds` (*) | Bonferroni, embudo de pares y tests de dependencia |
| `figura_quiaca.R` | `resumen_ajustes_filtrados.R` | Tendencia en la frecuencia de eventos (La Quiaca) y `tendencia_frecuencia.csv` |
| `periodo_retorno_via_magnitud.R` | corrida | Período de retorno (duración, intensidad) vía la cópula duración-magnitud y su validación (`via_magnitud_validacion.csv`) |
| `via_magnitud_resumen.R` | `periodo_retorno_via_magnitud.R` | Resumen de la validación |
| `figuras_via_magnitud.R` | `periodo_retorno_via_magnitud.R` | Isolíneas vía D-M (87585) y figuras de validación |
| `salvadas_via_magnitud.R` | `resumen_ajustes_filtrados.R`, `periodo_retorno_via_magnitud.R` | Estaciones con cópula D-M y sin I-D (`salvadas_via_magnitud.csv/.txt`) |
| `mapa_olas_copula.R` | `resumen_ajustes_filtrados.R` | Mapa de estaciones con cópula I-D en olas de calor y de frío |
| `mapa_olas_salvadas.R` | `salvadas_via_magnitud.R` | Mapas con las estaciones agregadas por D-M |
| `figuras_comentarios_nati.R`, `figuras_comentarios_nati_2.R` | `resumen_ajustes.R`, `regiones.csv`, serie diaria de 87585 (`SERIE_87585`) | Figuras y tablas ENSO/regiones de las rondas de comentarios (tablas en `tablas_enso/`) |

I-D: intensidad-duración; D-M: duración-magnitud.

(*) `est_todos_tests.rds` (p-valores de los 5 tests de estacionariedad sobre las series observadas) se calculó fuera del pipeline, porque el pipeline corta al primer test que rechaza; no hay un script versionado que lo genere.

## Salidas

Van a `DIR_PRES` (por defecto `data/output/presentacion`, ignorado por git); las figuras, a `DIR_PRES/figuras`.

## Notas

- Las figuras de ejemplo (familias, tipos de cópulas, parsimonia) son simulaciones ilustrativas, no resultados del pipeline.
- La cópula ajustada de 87585 usa los parámetros del pipeline (Gumbel con theta 1,1656; marginales GLO para la intensidad y PE3 para la duración).
- El período de retorno combinado es T = N / (n · P), con P = 1 − F_I(i) − F_D(d) + C(F_I(i), F_D(d)); es el mismo cálculo que `CalcularGrillaPeriodoRetorno` en `lib/funciones_periodo_retorno.R`.
