# Metodología Integral para el Monitoreo y Análisis de Sequías

## 1. Fundamentos Conceptuales
La sequía se define como una manifestación extrema y transitoria del ciclo hidrológico. Es fundamental distinguirla de la **aridez**, que es una condición climática permanente. Una sequía es un proceso con un inicio y un fin definidos que puede ocurrir en cualquier región, independientemente de su régimen hídrico habitual.

Para su análisis, el sistema se apoya en datos *in situ* provenientes de estaciones meteorológicas convencionales, complementados con sensores remotos que miden índices de vegetación (NDVI) y evapotranspiración (ALEXI).

---

## 2. Índices Estandarizados de Sequía
El pilar del monitoreo son los índices estandarizados, que permiten la comparación objetiva entre regiones con regímenes de precipitación muy distintos.

### Índice de Precipitación Estandarizado (SPI)
Se basa exclusivamente en datos de precipitación. Su cálculo implica ajustar los totales acumulados a una distribución de probabilidad (usualmente Gamma) para luego transformarlos en una distribución normal con media 0 y desviación estándar 1.

### Índice de Precipitación-Evapotranspiración Estandarizado (SPEI)
A diferencia del SPI, este índice considera el balance hídrico ($P - ETP$).
* **ETP (Evapotranspiración Potencial):** Se calcula mediante la fórmula de **Hargreaves (2002)**, que utiliza temperaturas máximas y mínimas. Este método es robusto frente a la falta de otros datos meteorológicos y captura la variabilidad de la demanda atmosférica.
* **Ajuste Estadístico:** Debido a que el balance hídrico puede arrojar valores negativos, se utiliza una distribución **Log-logística** para el ajuste.

---

## 3. Innovación en Resolución Temporal: El Uso de Pentadas
Para superar la limitación de los índices mensuales tradicionales, se emplea una agregación por **pentadas** (periodos de 5 días).
* **Estructura:** El mes se divide en 6 pentadas. Las primeras cinco son de 5 días; la última varía entre 3 y 6 días según el mes y el año.
* **Ventajas:** Proporciona **72 valores anuales** en lugar de 12. Esta alta frecuencia permite identificar el comienzo y fin de un evento seco con precisión de días, capturando variaciones extremas que el promedio mensual suele suavizar.

---

## 4. Flujo de Cálculo y Validación
El proceso de generación de índices sigue un protocolo riguroso:
1.  **Agregación Temporal:** Los datos de lluvia y temperatura se acumulan en la escala deseada (1, 3, 6 meses, etc.) usando ventanas móviles de pentadas.
2.  **Ajuste de Distribuciones:** Se aplican métodos de **Máxima Verosimilitud** o métodos **No Paramétricos** (basados en *splines*). Estos últimos son más flexibles cuando los datos no se ajustan bien a una curva teórica predefinida.
3.  **Manejo de Ceros:** En zonas semiáridas, la alta frecuencia de registros sin lluvia (ceros) puede sesgar el índice. Se aplica una corrección matemática para asegurar que la distribución se mantenga centrada en cero y conserve su interpretabilidad.
4.  **Tests de Bondad de Ajuste:** Se realizan pruebas estadísticas para validar que la función matemática elegida representa fielmente el comportamiento histórico de la estación meteorológica.

---

## 5. Identificación y Métricas de Eventos Secos
Un evento seco se identifica mediante un **umbral de detección** (ej. -0.5, -1.0 o -1.5) definido por el usuario según su objetivo. El evento comienza cuando el índice perfora dicho umbral y termina cuando lo supera nuevamente.

Cada evento se caracteriza mediante cuatro métricas clave:
* **Duración:** Número de periodos (pentadas o meses) que dura el evento.
* **Intensidad:** El valor promedio del índice durante el transcurso del evento.
* **Magnitud:** El déficit acumulado (suma de los valores del índice).
* **Valor Mínimo:** El punto de mayor severidad alcanzado durante el evento.

---

## 6. Análisis Probabilístico Avanzado
Para determinar el riesgo y el **periodo de retorno** (recurrencia) de una sequía, el sistema emplea dos herramientas de vanguardia:

### Generador Estocástico de Series Sintéticas
Dado que los registros históricos son a menudo insuficientes para capturar eventos de "una vez cada 100 años", se utiliza un modelo para generar **series sintéticas**. Esto permite simular miles de años de clima que mantienen las propiedades estadísticas de los datos reales, creando una base de datos masiva para el análisis de eventos extremos que aún no han ocurrido pero que son físicamente probables.

### Cúpulas Matemáticas (Análisis Multivariado)
Las métricas de sequía (como duración e intensidad) son dependientes entre sí. Las **Cúpulas** permiten unir estas variables en una sola distribución multivariada.
* **Periodo de Retorno Combinado:** Gracias a esto, es posible calcular la probabilidad de que ocurra una combinación específica de factores (ej. una sequía de 4 meses que además tenga una intensidad extrema).
* **Isolíneas de Probabilidad:** El resultado son gráficos donde se visualizan curvas de igual riesgo para distintas combinaciones de duración y severidad.

---

## 7. Aplicaciones Transversales
Este marco metodológico permite transformar datos abstractos en herramientas de gestión:
* **Agricultura:** Filtrar eventos por ventanas críticas (ej. floración de cultivos estivales) para calcular primas de seguros o riesgos de pérdida de rendimiento.
* **Recursos Hídricos:** Estimar la probabilidad de que un embalse no alcance su nivel crítico basándose en la recurrencia de la intensidad y duración de la sequía.
* **Sistemas de Alerta:** Utilizar la alta resolución de las pentadas para activar protocolos de emergencia de manera temprana.
