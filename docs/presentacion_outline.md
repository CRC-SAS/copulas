# Outline — Presentación "Cópulas para el análisis de olas de calor y de frío"

> **Estado:** v3, con las respuestas de la experta a las preguntas abiertas (2026-10-15). Todavía sin contenido final.
> **Formato final:** artifact HTML con estructura de charla (una sección = un bloque de slides).
> **Audiencia:** meteorólogos. Charla relajada, no es una clase de cópulas.
> **Duración estimada:** ~35 min + 10-15 min de preguntas.

**Cambios respecto de la v1**
- El hilo conductor son olas de calor y de frío (no sequía). `guia_video.md` queda como esqueleto: concepto → métricas → probabilidad → uso.
- Sklar solo se menciona en una frase, sin slide propio.
- Matemática "un poco": una slide de familias y colas, lenguaje de meteorólogo.
- Resultados: 1-2 casos representativos, no los 6 tipos de evento completos.
- Se incluyen los resultados ENSO.
- Se agrega un bloque nuevo: comparación por región y frío vs. calor (hay un primer resumen abajo).

---

## Bloque 0 — Apertura (2 min)

1. Título, autoría, financiamiento
2. Qué vamos a ver
   - Por qué mirar duración e intensidad por separado no alcanza
   - Qué es una cópula, en simple
   - Qué obtenemos y qué encontramos en frío/calor y en ENSO

## Bloque 1 — El problema (7 min)

3. **Una ola de calor (o de frío) no es un número**
   - Tiene duración, intensidad y magnitud (acumulado), definidas con umbrales de TX y TN
   - Un ejemplo concreto de una estación: dos eventos con la misma duración y muy distinta intensidad
4. **La pregunta que queremos responder**
   - "¿Cada cuánto ocurre una ola de **al menos 6 días y** de intensidad **al menos X**?"
   - El período de retorno de cada variable por separado no la responde
5. **Duración e intensidad no son independientes**
   - Un scatter real (ola de calor) con su tau de Kendall
   - Si las trato como independientes, el riesgo conjunto sale mal
6. **Pocos eventos por estación**
   - Decenas de olas por estación (mínimo 20 casos): hace falta un modelo parsimonioso

## Bloque 2 — Qué es una cópula (7 min)

7. **Idea central:** separar *cómo se distribuye cada variable* de *cómo se relacionan*
   - Marginales (GPA, Weibull, Normal, ...) + una función de dependencia
   - Mención de una línea al teorema de Sklar (garantiza que esta separación siempre existe)
8. **Intuición visual**
   - Scatter original → rangos (0,1) → la cópula
9. **Familias y colas: lo que importa para extremos**
   - Gaussiana, Frank: sin dependencia en las colas
   - Gumbel, Joe: los extremos altos "viajan juntos"
   - Clayton: dependencia en el extremo bajo
   - Una sola figura con 4 paneles; sin fórmulas

## Bloque 3 — Para qué sirve (5 min)

10. **Período de retorno combinado**
    - Probabilidad de superar duración *y* intensidad a la vez
11. **Isolíneas: cómo se lee el gráfico**
    - Un ejemplo: ¿cada cuánto una ola de ≥ d días con intensidad ≥ x?
    - Contraste con el período de retorno univariado del mismo evento
12. **Para qué le sirve a un meteorólogo**
    - Caracterizar olas históricas (¿qué tan raro fue este evento?)
    - Insumo para alertas y para estudios de impacto (salud, energía)

## Bloque 4 — Cómo lo hacemos (4 min)

13. **Pipeline en un diagrama:** eventos → marginales → estacionariedad → dependencia → cópula → período de retorno
14. **Qué decide el pipeline solo:** mejor marginal, mejor cópula (AIC/BIC, validación), solo si las series son estacionarias y dependientes
15. **Limitaciones:** pocos eventos, estaciones sin ajuste ("parámetro NA"), el test de bondad Sn

## Bloque 5 — Resultados: casos representativos (6 min)

16. **Caso 1 — Ola de calor en Buenos Aires (87585)** (elegida con la experta; región BUENOS AIRES / Pampeana)
    - 73 olas; intensidad-duración: Gumbel (cola superior), tau de Kendall 0,20
    - Scatter + mejor marginal + mejor cópula + isolíneas, coloreado por ENSO
17. **Caso 2 — Ola de frío, la misma estación (87585)**
    - 54 olas; intensidad-duración: Gumbel, tau de Kendall 0,05 (dependencia débil)
    - Misma lectura, para comparar calor vs. frío lado a lado
18. **Qué nos dicen las isolíneas**
    - Un evento histórico notable ubicado en el gráfico, con su período de retorno
18b. **Intensidad-magnitud, rápido** (87585, calor y frío)
    - Dependencia mucho más fuerte (tau 0,65-0,71) porque la magnitud es intensidad × duración; solo para mostrar que el pipeline también lo ajusta. El énfasis sigue en intensidad-duración

## Bloque 6 — ENSO (4 min)

19. **Qué agregamos:** cada evento se etiqueta con la fase ENSO vigente al inicio
    - Los puntos de los gráficos se colorean por Niño / Niña / Neutro
20. **Cópula sobre los episodios ENSO en sí** (duración-intensidad)
    - 27 Niño / 24 Niña; Niña: Joe, marginales Weibull/GPA; Niño: Normal, GPA/GPA
    - Limitante: n chico; clasificación por mes de inicio vs. categoría del episodio

## Bloque 7 — ¿Hay patrón por región o entre frío y calor? (5 min)

21. **Qué miramos:** familia de cópula ganadora (par intensidad-duración) por estación, agrupada por las 7 regiones de la experta (`regiones.csv`) y por tipo de evento
22. **Frío vs. calor**
    - Casi sin diferencias en la mezcla de familias (ver resumen abajo)
23. **Por región**
    - Mapa de estaciones coloreado por familia/tipo de cola; ¿se agrupan?
24. **Limitación explícita: pocos eventos por estación** (slide propia, pedida por la experta)
    - Olas: mediana de 39 eventos por estación (mínimo 20); períodos: mediana de 142
    - Con n chico, que gane Frank o Normal puede ser "no se distinguen"; no sobreinterpretar

## Bloque 8 — Cierre (2 min)

25. Síntesis: problema → cópula → utilidad → qué encontramos
26. Próximos pasos y puntos para discutir
27. Preguntas

---

## Resumen de región y frío/calor (exploratorio)

Calculado sobre `data/output/mejores_ajustes_multivariados_*.rds` (6 corridas con Olas_20), par **intensidad-duración**, 73 estaciones × 6 tipos = 422 ajustes (no todas las estaciones tienen los 6 tipos); 416 con familia ganadora (6 sin ajuste, 5 de ellos olas de frío). Regiones: las 7 de `regiones.csv` de la experta (NOA, NORTE, NEA, CENTRO, CUYO, BUENOS AIRES, PATAGONIA). Cuatro estaciones no estaban en su archivo y se asignaron por coordenadas (a confirmar con ella): 87345 → CENTRO, 87420 → CUYO, 87569 → BUENOS AIRES, 87934 → PATAGONIA. Los porcentajes son sobre ajustes (estación × tipo).

**Por tipo de cola de la cópula ganadora (% de ajustes)**

| | Ajustes | Sin cola (Normal/Frank) | Cola superior (Gumbel/Joe) | Cola inferior (Clayton) | Ambas (t) |
|---|---|---|---|---|---|
| Calor (olas + períodos cálidos) | 215 | 60 | 23 | 12 | 6 |
| Frío (olas + períodos fríos) | 201 | 63 | 21 | 9 | 6 |
| Olas (calor + frío) | 125 | 43 | 22 | 20 | 14 |
| Períodos (cálidos + fríos) | 291 | 69 | 22 | 7 | 2 |

**Por región (todos los tipos)**

| Región | Ajustes | Sin cola | Cola superior | Cola inferior | Ambas | n eventos (mediana) | tau (mediana) |
|---|---|---|---|---|---|---|---|
| NOA | 38 | 87 | 8 | 5 | 0 | 141 | 0,23 |
| NORTE | 20 | 65 | 25 | 5 | 5 | 144 | 0,17 |
| NEA | 80 | 72 | 11 | 10 | 6 | 150 | 0,20 |
| CUYO | 42 | 57 | 24 | 10 | 10 | 136 | 0,24 |
| CENTRO | 91 | 59 | 23 | 13 | 4 | 123 | 0,20 |
| BUENOS AIRES | 82 | 49 | 32 | 13 | 6 | 128 | 0,18 |
| PATAGONIA | 63 | 52 | 29 | 10 | 10 | 98 | 0,22 |

**Frío vs. calor por región (% sin cola / % con cola)**

| Región | Calor | Frío |
|---|---|---|
| NOA | 85 / 15 | 89 / 11 |
| NORTE | 58 / 42 | 75 / 25 |
| NEA | 71 / 29 | 74 / 26 |
| CUYO | 43 / 57 | 71 / 29 |
| CENTRO | 58 / 42 | 60 / 40 |
| BUENOS AIRES | 56 / 44 | 41 / 59 |
| PATAGONIA | 48 / 52 | 56 / 44 |

**Lectura**
- **Frío vs. calor:** prácticamente iguales en el total (2-3 puntos). Por región hay diferencias puntuales (CUYO más cola en calor, BUENOS AIRES más cola en frío), pero con 8-48 ajustes por celda no son concluyentes.
- **Región:** gradiente NOA/NEA (más Normal/Frank) → BUENOS AIRES/PATAGONIA/CUYO (más cola superior). NORTE son solo 20 ajustes.
- **Diferencia más marcada:** entre *olas* (57% con cola) y *períodos* (31%), más que entre frío y calor o entre regiones. Es también donde hay menos eventos (olas: mediana 39 por estación; períodos: 142), así que puede reflejar poder estadístico y no solo el fenómeno.
- **Tau de Kendall** (intensidad-duración): mediana ~0,2 en todas las regiones; dependencia débil a moderada, por lo que Normal/Frank pueden ser indistinguibles de las colas con n chico.
- **Cautela:** el ganador por AIC con pocos eventos puede no ser estadísticamente distinguible del segundo. Entre los 36 ajustes con n < 30 (todos olas), 19 ganan Normal/Frank.

---

## Decisiones acordadas con la experta (2026-10-15)

1. **Caso representativo:** estación 87585 (Buenos Aires, región BUENOS AIRES/Pampeana), calor y frío.
2. **Regionalización:** las de su `regiones.csv`, sin agrupar más (BUENOS AIRES y CENTRO se mantienen separadas).
3. **Magnitud:** énfasis en intensidad-duración; intensidad-magnitud solo como slide rápida.
4. **Limitación de n chico:** explícita, en slide propia.
5. **Estaciones fuera de su archivo:** asignación confirmada (87345 → CENTRO, 87420 → CUYO, 87569 → BUENOS AIRES, 87934 → PATAGONIA).
