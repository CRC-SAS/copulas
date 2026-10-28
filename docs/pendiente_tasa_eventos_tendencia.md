# Pendiente: tasa de eventos con tendencia en la frecuencia (no incluido en el pipeline)

Estado (2026-10-28): decidido NO incluirlo en el pipeline por ahora. Solo se usa en la presentacion.
Este documento guarda lo necesario para incorporarlo mas adelante.

## Problema
El periodo de retorno combinado es T = N / (n * P), con N anios de registro, n eventos y P la
probabilidad de excedencia de la copula (`lib/funciones_periodo_retorno.R`, ~linea 317). Supone una
tasa de eventos constante n/N. Si la frecuencia de eventos tiene tendencia (ej. La Quiaca, olas de
calor: 25 olas en 65 anios, Mann-Kendall sobre conteos anuales p < 0,001), esa tasa subestima la
actual y T queda sobreestimado. P no depende de la tasa.

## Que se hizo en la presentacion
- Test de Mann-Kendall (paquete Kendall) sobre el conteo anual de eventos (1961-2025, anios sin eventos = 0).
- Si p < 0,01: tasa = (eventos entre 1991 y 2020) / 30, T = 1 / (tasa * P). Si no: n/N (como el pipeline).
- Casos: Buenos Aires 87585, PeriodoCalidoTX (p = 0,0016; 66/30 = 2,20 por anio vs n/N = 2,21) y olas de calor
  via copula duracion-magnitud (p < 0,001; 35/30). TX frio (p = 0,017) no la usa. Efecto: dic 2013 pasa de
  107 a 108 anios (via D-M: 120 a 123).
- La Quiaca: tasa 1991-2020 = 0,67 por anio (20 olas en 30 anios), T = 7,5 anios con P = 0,2 (13 con todo el registro).
- Codigo: `scripts/presentacion/figuras_caso_87585.R` (~lineas 33-38), `figuras_via_magnitud.R` (~lineas 26-29),
  `figura_quiaca.R` (figura y conteo de pares con tendencia).

## Decisiones abiertas (a definir con Natalia)
- Umbral del test: 0,01 (usado) o 0,05. Con 0,01 hay pocos casos; ver `tendencia_frecuencia.csv`
  (salida de `figura_quiaca.R`) para cuantos de los 196 pares con copula tienen tendencia.
- Ventana: la normal climatologica 1991-2020 (usada) u otra (ultimos 10 o 30 anios, ventana movil).
- Como marcarlo: columna en los resultados con la tasa usada y el p-valor.
- Efecto en otros resultados: cambia T en las estaciones con tendencia, hay que regenerar los zips.

## Como implementarlo
1. Parametros en `parametros_copulas.yml` (propuesta): `tasa.eventos.tendencia.alfa` (0,01),
   `tasa.eventos.periodo.inicio` / `.fin` (1991 / 2020), y un interruptor (por defecto desactivado
   para no cambiar resultados existentes).
2. Funcion en `lib/funciones_periodo_retorno.R` que calcule el test y la tasa; usarla donde hoy se calcula N/(n*P).
   Ojo: los archivos de `lib/*.R` tienen comentarios mezclados Latin-1/UTF-8: editar a nivel de bytes, solo ASCII.
3. Columnas nuevas en los resultados: tasa usada, origen (n/N o 1991-2020) y p-valor de tendencia.
4. Test en `tests/` y un commit aparte.
