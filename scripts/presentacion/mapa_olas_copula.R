# Mapa de las estaciones con copula ganadora duracion-intensidad en olas de calor (rojo) y de frio (azul);
# circulo vacio donde hay olas pero la copula no paso los filtros. Ninguna estacion tiene copula en ambos tipos.
# Se corre desde la raiz del repo, despues de resumen_ajustes.R y resumen_ajustes_filtrados.R
# Entradas: DIR_PRES/id_filtrado.csv, data/input/regiones.csv (Latin-1). Salida: DIR_PRES/figuras/mapa_olas_copula.png
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
REGIONES_CSV <- Sys.getenv("REGIONES_CSV", "data/input/regiones.csv")
suppressMessages({library(dplyr); library(ggplot2)})
tema <- theme_minimal(base_size = 15) + theme(panel.grid = element_blank(), axis.text = element_blank(), plot.title = element_text(face = "bold", size = 16),
  plot.subtitle = element_text(size = 12, colour = "#56606B"), plot.title.position = "plot", plot.background = element_rect(fill = "white", colour = NA), legend.position = "bottom")
id <- read.csv(file.path(DIR_PRES, "id_filtrado.csv")) %>% filter(clase == "Ola") %>% mutate(station_id = as.character(station_id))
rg <- read.table(REGIONES_CSV, sep = ";", header = TRUE, dec = ".", quote = "", fileEncoding = "latin1", stringsAsFactors = FALSE)
co <- data.frame(station_id = as.character(rg[[1]]), lon = rg[[4]], lat = rg[[5]])
co <- bind_rows(co, data.frame(station_id = c("87345", "87420", "87569", "87934"), lon = c(-64.198, -68.873, -58.731, -67.759), lat = c(-31.421, -32.895, -34.555, -53.780)))
d <- id %>% left_join(co, by = "station_id") %>% group_by(station_id, lon, lat) %>%
  summarise(estado = case_when(any(gana & estacion_tipo == "Calor") ~ "Cópula en olas de calor", any(gana & estacion_tipo == "Frio") ~ "Cópula en olas de frío", TRUE ~ "Sin cópula en olas"), .groups = "drop")
stopifnot(!anyNA(d$lon))
cat("estaciones por estado:\n"); print(table(d$estado))
d$estado <- factor(d$estado, levels = c("Cópula en olas de calor", "Cópula en olas de frío", "Sin cópula en olas"))
arg <- ggplot2::map_data("world", region = "Argentina")
p <- ggplot() + geom_polygon(data = arg, aes(long, lat, group = group), fill = "#EEF2F5", colour = "grey60", linewidth = .3) +
  geom_point(data = filter(d, estado == "Sin cópula en olas"), aes(lon, lat, colour = estado), shape = 1, size = 3, stroke = .8) +
  geom_point(data = filter(d, estado != "Sin cópula en olas"), aes(lon, lat, colour = estado), shape = 16, size = 4.6) +
  scale_colour_manual(values = c("Cópula en olas de calor" = "#C0392B", "Cópula en olas de frío" = "#1F5F8B", "Sin cópula en olas" = "grey45"), name = NULL, drop = FALSE) +
  guides(colour = guide_legend(nrow = 3, override.aes = list(shape = c(16, 16, 1), size = c(4.6, 4.6, 3)))) +
  coord_quickmap(xlim = c(-74, -53), ylim = c(-55.5, -21.5)) +
  labs(x = NULL, y = NULL, title = "Cópula en las olas", subtitle = "Duración–intensidad; ninguna estación la tiene en las dos") + tema
ggsave(file.path(DIR_PRES, "figuras", "mapa_olas_copula.png"), p, width = 5.6, height = 7.6, dpi = 130)
