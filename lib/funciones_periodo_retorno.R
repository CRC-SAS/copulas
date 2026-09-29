# ------------------------------------------------------------------------------#
# ---- Funciones para el periodo de retorno combinado (co-occurrence) de una ----
# ---- copula ya ajustada, analogo a la Figura 7 de Chen et al. 2024         ----
# ------------------------------------------------------------------------------#

# Paleta divergente rojo (Nino) <-> azul (Nina), usada para colorear el
# relleno de los puntos observados en GraficarPeriodoRetorno segun el evento
# ENSO vigente al inicio de cada evento. Validada con el script de validacion
# de paletas (separacion CVD/vision normal >= umbral entre las ramas roja y
# azul; rampas monotonas en luminosidad dentro de cada rama) - ver
# docs/superpowers si se necesita regenerar. Neutro va sin relleno (solo el
# borde negro del punto) para que los eventos ENSO resalten sobre el heatmap.
COLORES_ENSO <- c(
  "Niño Muy Fuerte" = "#a50f15",
  "Niño Fuerte"     = "#de2d26",
  "Niño Moderado"   = "#fb6a4a",
  "Niño Débil" = "#fcae91",
  "Neutro"               = "transparent",
  "Niña Débil" = "#bdd7e7",
  "Niña Moderada"   = "#6baed6",
  "Niña Fuerte"     = "#08519c"
)
NIVELES_ENSO <- names(COLORES_ENSO)
# Neutro se dibuja como cuadrado (22) y el resto como circulo (21): asi los
# neutros (sin relleno) no se confunden con los NA (circulo blanco)
FORMAS_ENSO <- stats::setNames(ifelse(NIVELES_ENSO == "Neutro", 22, 21), NIVELES_ENSO)
#
# Formula (caso AND / co-occurrence, ambas variables superan simultaneamente
# el umbral x,y):
#   T(x,y) = N / (n * (1 - F_X(x) - F_Y(y) + C(F_X(x), F_Y(y))))
# donde N es la extension del registro en anios, n la cantidad de eventos,
# F_X/F_Y las marginales ajustadas y C la copula ajustada. Se reutiliza el
# objeto mvdc (copula + ambas marginales) ya calculado por el pipeline.

# Primera letra en mayuscula, para nombres de variables/familias en titulos y ejes
Capitalizar <- function(texto) {
  paste0(toupper(substr(texto, 1, 1)), substring(texto, 2))
}

# Etiqueta legible de una variable para titulos/ejes (los nombres internos van sin tilde)
EtiquetaVariable <- function(variable) {
  Capitalizar(dplyr::recode(variable, duracion = "duración"))
}

CalcularGrillaPeriodoRetorno <- function(mvdc, N, n, grid_x, grid_y) {
  # La familia ganadora puede tener parametro NA (ver AplicarMejorAjusteACopulas
  # / fc17f05 / 373e1df). Para clayton/normal/t, pCopula/pMvdc abortan con
  # error ante un parametro NA (queda atrapado por el tryCatch del caller,
  # se degrada a archivo_png=NA como corresponde). Para joe/gumbel, en cambio,
  # pMvdc devuelve NaN en silencio en vez de abortar: sin este chequeo
  # explicito, la grilla completa queda en NA pero NO se detecta como error,
  # y se termina generando un PNG "valido" en apariencia (con titulo y
  # familia) pero con el heatmap/isolineas completamente vacios. Se fuerza
  # el mismo tratamiento (stop(), atrapado por el caller) para todas las
  # familias por igual.
  if (any(is.na(mvdc@copula@parameters))) {
    stop("parameter is NA")
  }

  grilla <- tidyr::crossing(x = grid_x, y = grid_y)

  F_X  <- do.call(what = paste0("p", mvdc@margins[1]), args = c(list(q = grilla$x), mvdc@paramMargins[[1]]))
  F_Y  <- do.call(what = paste0("p", mvdc@margins[2]), args = c(list(q = grilla$y), mvdc@paramMargins[[2]]))
  F_XY <- copula::pMvdc(cbind(grilla$x, grilla$y), mvdc)

  # P(X>=x, Y>=y), acotada para evitar T infinito/negativo por errores de redondeo
  prob_conjunta <- pmax(1 - F_X - F_Y + F_XY, 1e-6)
  grilla$T <- N / (n * prob_conjunta)

  # Red de seguridad adicional: si por cualquier otra razon (no solo
  # parametro NA) la grilla entera queda sin valores validos, tratarlo
  # tambien como fallo en vez de dejar pasar un grafico vacio.
  if (all(is.na(grilla$T))) {
    stop("todos los valores de T resultaron NA/NaN")
  }

  return(grilla)
}

GraficarPeriodoRetorno <- function(grilla, niveles_anios, x_obs, y_obs, enso_obs,
                                    nombre_x, nombre_y, titulo, archivo_png) {
  grid_x <- sort(unique(grilla$x))
  grid_y <- sort(unique(grilla$y))

  # Construccion de T_mat indexando por posicion (match), no por orden de
  # filas de grilla: asumir que grilla esta ordenada como un crossing(x, y)
  # con x variando mas lento (lo que matrix(grilla$T, nrow=, ncol=) requeria)
  # es fragil, y de hecho estaba mal (quedaba transpuesta) cuando la funcion
  # de graficacion recibe una tabla con las columnas x/y ya intercambiadas
  # (ver CalcularPeriodoRetornoUC, inversion de eje para "duracion"). Con
  # match() la matriz queda correcta sin importar el orden de filas de grilla.
  T_mat <- matrix(NA_real_, nrow = length(grid_x), ncol = length(grid_y))
  T_mat[cbind(match(grilla$x, grid_x), match(grilla$y, grid_y))] <- grilla$T

  # Posicion de las etiquetas de cada curva de nivel: punto medio de cada linea
  lineas <- grDevices::contourLines(grid_x, grid_y, T_mat, levels = niveles_anios)
  etiquetas <- purrr::map_dfr(lineas, function(l) {
    medio <- ceiling(length(l$x) / 2)
    tibble::tibble(nivel = l$level, x = l$x[medio], y = l$y[medio])
  })

  # evento_enso puede venir NA (indeterminado, fecha_inicio fuera del rango
  # del archivo ENSO) - se factoriza con todos los niveles conocidos para que
  # la leyenda muestre siempre el mismo orden/colores entre graficos, y solo
  # aparezca la entrada "NA" si realmente hay algun punto sin clasificar.
  observados <- tibble::tibble(x = x_obs, y = y_obs,
                               evento_enso = factor(enso_obs, levels = NIVELES_ENSO))

  p <- ggplot2::ggplot(grilla, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = T), interpolate = TRUE) +
    ggplot2::geom_contour(ggplot2::aes(z = T), breaks = niveles_anios,
                          color = "white", linewidth = 0.4) +
    ggplot2::scale_fill_viridis_c(trans = "log10", breaks = niveles_anios,
                                  limits = range(niveles_anios), oob = scales::squish,
                                  name = "Período de\nretorno (años)",
                                  guide = ggplot2::guide_colourbar(
                                    direction = "vertical", order = 1,
                                    theme = ggplot2::theme(legend.key.height = ggplot2::unit(10, "lines")))) +
    # El heatmap de arriba ya usa la estetica "fill" (escala continua viridis);
    # ggnewscale permite una segunda escala "fill" independiente para los
    # puntos observados (escala discreta ENSO), en vez de compartir una sola
    # escala de fill para todo el grafico
    ggnewscale::new_scale_fill() +
    # fill y shape comparten name/breaks, por lo que ggplot fusiona ambas
    # escalas en una unica leyenda "Intensidad ENSO"
    ggplot2::geom_point(data = observados, ggplot2::aes(x = x, y = y, fill = evento_enso, shape = evento_enso),
                        color = "black", size = 1.8, stroke = 0.4) +
    ggplot2::scale_fill_manual(values = COLORES_ENSO, breaks = NIVELES_ENSO,
                               name = "Intensidad ENSO",
                               na.value = "white", drop = TRUE,
                               guide = ggplot2::guide_legend(ncol = 1, order = 2)) +
    ggplot2::scale_shape_manual(values = FORMAS_ENSO, breaks = NIVELES_ENSO,
                                name = "Intensidad ENSO",
                                na.value = 21, drop = TRUE,
                                guide = ggplot2::guide_legend(ncol = 1, order = 2)) +
    ggplot2::labs(x = EtiquetaVariable(nombre_x), y = EtiquetaVariable(nombre_y), title = titulo) +
    ggplot2::theme_minimal(base_size = 12) +
    # Ambas leyendas apiladas a la derecha (en vez de abajo) y margenes
    # minimos, para dejarle al panel la mayor area posible
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   plot.title.position = "plot",
                   legend.position = "right", legend.box = "vertical",
                   legend.justification = "center",
                   legend.box.spacing = ggplot2::unit(4, "pt"),
                   legend.spacing.y = ggplot2::unit(8, "pt"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   plot.margin = ggplot2::margin(4, 4, 4, 4))

  # La duracion es conceptualmente entera (dias), pero se trata como continua
  # en la grilla de evaluacion: sin esto, los breaks automaticos de ggplot
  # eligen incrementos "redondos" (2.5, 7.5, ...) que no tienen sentido para
  # una duracion. Se fuerza a mostrar todos los enteros del rango graficado.
  # expand = 0 elimina el padding de ejes: el heatmap ocupa todo el panel.
  breaks_x <- if (identical(nombre_x, "duracion")) seq(floor(min(grid_x)), ceiling(max(grid_x)), by = 1) else ggplot2::waiver()
  breaks_y <- if (identical(nombre_y, "duracion")) seq(floor(min(grid_y)), ceiling(max(grid_y)), by = 1) else ggplot2::waiver()
  p <- p + ggplot2::scale_x_continuous(breaks = breaks_x, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = breaks_y, expand = c(0, 0)) +
    # clip = "off": con expand = 0 los puntos observados sobre el borde del
    # panel (la grilla arranca en el minimo observado) quedarian cortados
    ggplot2::coord_cartesian(clip = "off")

  if (nrow(etiquetas) > 0) {
    p <- p + ggplot2::geom_label(data = etiquetas, ggplot2::aes(x = x, y = y, label = nivel),
                                 size = 3, linewidth = 0, label.padding = ggplot2::unit(0.12, "lines"),
                                 fill = grDevices::rgb(1, 1, 1, 0.75), color = "gray20")
  }

  ggplot2::ggsave(filename = archivo_png, plot = p, width = 7.5, height = 5.5, dpi = 150)
}

# ------------------------------------------------------------------------------#
# ---- Funciones para el periodo de retorno univariado de una variable       ----
# ---- individual (intensidad/magnitud/duracion), analogo al caso bivariado  ----
# ---- de arriba pero con una sola marginal                                  ----
# ------------------------------------------------------------------------------#
#
# Formula: T(x) = N / (n * (1 - F(x))), donde F es la distribucion ganadora
# del mejor ajuste univariado (PASO 6) para esa estacion+variable.

CalcularGrillaPeriodoRetornoUV <- function(distribucion, parametros, N, n, grid_x) {
  F_X <- do.call(what = paste0("p", distribucion), args = c(list(q = grid_x), parametros))

  # P(X>=x), acotada para evitar T infinito/negativo por errores de redondeo
  prob_excedencia <- pmax(1 - F_X, 1e-6)
  grilla <- tibble::tibble(x = grid_x, T = N / (n * prob_excedencia))

  return(grilla)
}

GraficarPeriodoRetornoUV <- function(grilla, niveles_anios, x_obs, N, n,
                                     nombre_x, titulo, archivo_png) {
  # Posicion de graficacion empirica (Weibull) de los eventos observados: para
  # el evento de rango m (1 = valor mas alto, hasta n), T_empirico = N*(n+1)/(n*m).
  # No se usa la F teorica para no forzar que caigan siempre exactos sobre la
  # curva: la distancia entre puntos y curva es la senal visual de bondad de ajuste.
  m <- rank(-x_obs)
  observados <- tibble::tibble(T = N * (n + 1) / (n * m), x = x_obs)

  p <- ggplot2::ggplot(grilla, ggplot2::aes(x = T, y = x)) +
    ggplot2::geom_vline(xintercept = niveles_anios, color = "gray80", linewidth = 0.3) +
    ggplot2::geom_line(color = "steelblue4", linewidth = 0.8) +
    ggplot2::geom_point(data = observados, ggplot2::aes(x = T, y = x),
                        shape = 21, fill = "white", color = "black", size = 1.8) +
    ggplot2::scale_x_log10(breaks = niveles_anios, labels = niveles_anios) +
    ggplot2::labs(x = "Período de retorno (años)", y = nombre_x, title = titulo) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  # Misma razon que en GraficarPeriodoRetorno: duracion es conceptualmente
  # entera (dias), forzar breaks enteros en el eje de valores si corresponde.
  if (identical(nombre_x, "duracion")) {
    p <- p + ggplot2::scale_y_continuous(breaks = seq(floor(min(grilla$x)), ceiling(max(grilla$x)), by = 1))
  }

  ggplot2::ggsave(filename = archivo_png, plot = p, width = 7.5, height = 5.5, dpi = 150)
}

# ------------------------------------------------------------------------------#
# ---- Grafico opcional (parametro periodo_retorno.graficar_distribucion_    ----
# ---- ajuste): histograma de los eventos observados + densidad de la        ----
# ---- distribucion ganadora superpuesta, analogo a la Figura 5 de           ----
# ---- Chen et al. 2024                                                      ----
# ------------------------------------------------------------------------------#

GraficarDistribucionUnivariada <- function(x_obs, distribucion, parametros,
                                           nombre_x, titulo, archivo_png) {
  observados <- tibble::tibble(x = x_obs)
  densidad <- function(x) do.call(what = paste0("d", distribucion), args = c(list(x = x), parametros))

  # Duracion es conceptualmente entera (dias): un bins=N generico calcula el
  # ancho de bin como rango/(N-1), que solo por casualidad coincide con un
  # ancho que alinee los bins a los enteros (depende de que el rango
  # observado en esa estacion sea multiplo de ese ancho). Si no coincide, las
  # barras quedan visualmente corridas respecto de los enteros del eje (breaks
  # forzados mas abajo), aunque el eje si este alineado. Se fuerza en cambio
  # binwidth=1 con boundary=0.5, que garantiza un bin por cada entero sin
  # importar el rango de cada estacion. Para intensidad/magnitud (continuas)
  # se mantiene el bins=15 generico.
  histograma <- if (identical(nombre_x, "duracion")) {
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            binwidth = 1, boundary = 0.5,
                            fill = "steelblue3", color = "white", alpha = 0.7)
  } else {
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = 15, fill = "steelblue3", color = "white", alpha = 0.7)
  }

  p <- ggplot2::ggplot(observados, ggplot2::aes(x = x)) +
    histograma +
    ggplot2::stat_function(fun = densidad, color = "steelblue4", linewidth = 0.9) +
    ggplot2::labs(x = nombre_x, y = "Densidad", title = titulo) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  # Misma razon que en las demas funciones de graficacion: forzar breaks
  # enteros en el eje si corresponde.
  if (identical(nombre_x, "duracion")) {
    p <- p + ggplot2::scale_x_continuous(breaks = seq(floor(min(x_obs)), ceiling(max(x_obs)), by = 1))
  }

  ggplot2::ggsave(filename = archivo_png, plot = p, width = 7.5, height = 5.5, dpi = 150)
}
