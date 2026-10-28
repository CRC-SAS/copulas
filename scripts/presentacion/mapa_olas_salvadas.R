# Mapas de olas de calor y de frio: estaciones con copula intensidad-duracion (I-D), estaciones que se agregan con la copula duracion-magnitud (D-M)
# y estaciones sin ninguna de las dos (circulo vacio). I-D: intensidad-duracion; D-M: duracion-magnitud.
# Se corre desde la raiz del repo, despues de salvadas_via_magnitud.R
# Entradas: DIR_PRES/salvadas_via_magnitud.csv, data/input/regiones.csv (Latin-1). Salida: DIR_PRES/figuras/mapa_olas_salvadas.png
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
REGIONES_CSV <- Sys.getenv("REGIONES_CSV", "data/input/regiones.csv")
suppressMessages({library(dplyr); library(ggplot2); library(patchwork)})
tema <- theme_minimal(base_size = 15) + theme(panel.grid = element_blank(), axis.text = element_blank(), plot.title = element_text(face = "bold", size = 17),
  plot.subtitle = element_text(size = 12.5, colour = "#56606B"), plot.title.position = "plot", plot.background = element_rect(fill = "white", colour = NA), legend.position = "bottom")
rg <- read.table(REGIONES_CSV, sep = ";", header = TRUE, dec = ".", quote = "", fileEncoding = "latin1", stringsAsFactors = FALSE)
co <- data.frame(station_id = as.character(rg[[1]]), lon = rg[[4]], lat = rg[[5]])
co <- bind_rows(co, data.frame(station_id = c("87345", "87420", "87569", "87934"), lon = c(-64.198, -68.873, -58.731, -67.759), lat = c(-31.421, -32.895, -34.555, -53.780)))
s <- read.csv(file.path(DIR_PRES, "salvadas_via_magnitud.csv")) %>% filter(slug %in% c("olacalor", "olafrio")) %>% mutate(station_id = as.character(station_id)) %>% left_join(co, by = "station_id")
stopifnot(!anyNA(s$lon))
arg <- ggplot2::map_data("world", region = "Argentina")
mapa <- function(sl, titulo, col_id, col_dm) {
  d <- s %>% filter(slug == sl) %>% mutate(estado = factor(case_when(di ~ "Con I–D", solo_dm ~ "Se agrega con D–M", TRUE ~ "Sin cópula"),
                                                            levels = c("Con I–D", "Se agrega con D–M", "Sin cópula")))
  n <- table(d$estado); cat(titulo, ":", paste(names(n), n, sep = " = ", collapse = "; "), "| total", nrow(d), "\n")
  ggplot() + geom_polygon(data = arg, aes(long, lat, group = group), fill = "#EEF2F5", colour = "grey60", linewidth = .3) +
    geom_point(data = filter(d, estado == "Sin cópula"), aes(lon, lat, colour = estado), shape = 1, size = 3, stroke = .8) +
    geom_point(data = filter(d, estado != "Sin cópula"), aes(lon, lat, colour = estado), shape = 16, size = 4.8) +
    scale_colour_manual(values = setNames(c(col_id, col_dm, "grey45"), levels(d$estado)), name = NULL, drop = FALSE) +
    guides(colour = guide_legend(nrow = 3, override.aes = list(shape = c(16, 16, 1), size = c(4.8, 4.8, 3)))) +
    coord_quickmap(xlim = c(-74, -53), ylim = c(-55.5, -21.5)) +
    labs(x = NULL, y = NULL, title = titulo, subtitle = sprintf("%d de %d estaciones con I–D;\n%d más con D–M", sum(d$di), nrow(d), sum(d$solo_dm))) + tema
}
p <- mapa("olacalor", "Olas de calor", "#C0392B", "#F29BBE") + mapa("olafrio", "Olas de frío", "#1F5F8B", "#8CCBF0")
ggsave(file.path(DIR_PRES, "figuras", "mapa_olas_salvadas.png"), p, width = 7.4, height = 7.9, dpi = 130)
