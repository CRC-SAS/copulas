# ------------------------------------------------------------------------------#
# ---- Funciones para el periodo de retorno combinado (co-occurrence) de una ----
# ---- copula ya ajustada, analogo a la Figura 7 de Chen et al. 2024         ----
# ------------------------------------------------------------------------------#
#
# Formula (caso AND / co-occurrence, ambas variables superan simultaneamente
# el umbral x,y):
#   T(x,y) = N / (n * (1 - F_X(x) - F_Y(y) + C(F_X(x), F_Y(y))))
# donde N es la extension del registro en anios, n la cantidad de eventos,
# F_X/F_Y las marginales ajustadas y C la copula ajustada. Se reutiliza el
# objeto mvdc (copula + ambas marginales) ya calculado por el pipeline.

CalcularGrillaPeriodoRetorno <- function(mvdc, N, n, grid_x, grid_y) {
  grilla <- tidyr::crossing(x = grid_x, y = grid_y)

  F_X  <- do.call(what = paste0("p", mvdc@margins[1]), args = c(list(q = grilla$x), mvdc@paramMargins[[1]]))
  F_Y  <- do.call(what = paste0("p", mvdc@margins[2]), args = c(list(q = grilla$y), mvdc@paramMargins[[2]]))
  F_XY <- copula::pMvdc(cbind(grilla$x, grilla$y), mvdc)

  # P(X>=x, Y>=y), acotada para evitar T infinito/negativo por errores de redondeo
  prob_conjunta <- pmax(1 - F_X - F_Y + F_XY, 1e-6)
  grilla$T <- N / (n * prob_conjunta)

  return(grilla)
}

GraficarPeriodoRetorno <- function(grilla, niveles_anios, x_obs, y_obs,
                                    nombre_x, nombre_y, titulo, archivo_png) {
  grid_x <- sort(unique(grilla$x))
  grid_y <- sort(unique(grilla$y))
  T_mat  <- matrix(grilla$T, nrow = length(grid_x), ncol = length(grid_y))

  # Posicion de las etiquetas de cada curva de nivel: punto medio de cada linea
  lineas <- grDevices::contourLines(grid_x, grid_y, T_mat, levels = niveles_anios)
  etiquetas <- purrr::map_dfr(lineas, function(l) {
    medio <- ceiling(length(l$x) / 2)
    tibble::tibble(nivel = l$level, x = l$x[medio], y = l$y[medio])
  })

  observados <- tibble::tibble(x = x_obs, y = y_obs)

  p <- ggplot2::ggplot(grilla, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = T), interpolate = TRUE) +
    ggplot2::geom_contour(ggplot2::aes(z = T), breaks = niveles_anios,
                          color = "white", linewidth = 0.4) +
    ggplot2::geom_point(data = observados, ggplot2::aes(x = x, y = y),
                        shape = 21, fill = "white", color = "black", size = 1.8) +
    ggplot2::scale_fill_viridis_c(trans = "log10", breaks = niveles_anios,
                                  limits = range(niveles_anios), oob = scales::squish,
                                  name = "Período de\nretorno (años)") +
    ggplot2::labs(x = nombre_x, y = nombre_y, title = titulo) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  if (nrow(etiquetas) > 0) {
    p <- p + ggplot2::geom_label(data = etiquetas, ggplot2::aes(x = x, y = y, label = nivel),
                                 size = 3, linewidth = 0, label.padding = ggplot2::unit(0.12, "lines"),
                                 fill = grDevices::rgb(1, 1, 1, 0.75), color = "gray20")
  }

  ggplot2::ggsave(filename = archivo_png, plot = p, width = 7.5, height = 5.5, dpi = 150)
}
