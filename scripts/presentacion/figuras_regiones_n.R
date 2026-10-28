# Figuras: cola por grupo, mapa de estaciones y cola segun n de eventos (requiere resumen_ajustes.R)
# Se corre desde la raiz del repo: Rscript scripts/presentacion/figuras_regiones_n.R
# Entradas: data/output/*.rds (corrida Olas_20), data/input/eventos_identificados_frio_calor.csv
# Salidas: DIR_PRES (por defecto data/output/presentacion); las figuras van en DIR_PRES/figuras
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
dir.create(file.path(DIR_PRES, "figuras"), recursive = TRUE, showWarnings = FALSE)
REGIONES_CSV <- Sys.getenv("REGIONES_CSV", "data/input/regiones.csv")
suppressMessages({library(dplyr); library(ggplot2)})
out <- file.path(DIR_PRES, "figuras")
azul <- "#1F5F8B"; calor <- "#B4491F"; gris <- "#56606B"
tema <- theme_minimal(base_size = 15) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15),
  plot.subtitle = element_text(size = 12, colour = gris), plot.title.position = "plot", plot.background = element_rect(fill = "white", colour = NA))
ev <- readr::read_csv("data/input/eventos_identificados_frio_calor.csv", show_col_types = FALSE) %>% filter(duracion >= 3) %>%
  mutate(station_id = as.character(station_id)) %>% count(station_id, tipo = tipo_evento, name = "n")
r <- readr::read_csv(file.path(DIR_PRES, "resumen_ajustes.csv"), show_col_types = FALSE) %>%
  filter(variable_x == "intensidad", variable_y == "duracion", !is.na(familia)) %>% mutate(station_id = as.character(station_id)) %>% select(-n) %>%
  left_join(ev, by = c("station_id", "tipo")) %>% mutate(con_cola = cola != "Sin cola")
cat("ajustes:", nrow(r), "\n")

# 23: apiladas por grupo
niv <- c("Sin cola", "Superior", "Inferior", "Ambas")
lab <- c("Sin cola", "Cola superior", "Cola inferior", "Ambas")
g <- bind_rows(r %>% mutate(grupo = ifelse(estacion_tipo == "Calor", "Calor", "Frío")), r %>% mutate(grupo = ifelse(clase == "Ola", "Olas", "Períodos"))) %>%
  count(grupo, cola) %>% group_by(grupo) %>% mutate(p = 100 * n / sum(n), tot = sum(n)) %>% ungroup() %>%
  mutate(cola = factor(cola, levels = rev(niv), labels = rev(lab)), grupo = factor(grupo, levels = c("Calor", "Frío", "Olas", "Períodos")),
         etq = paste0(grupo, "\n(", tot, " ajustes)"))
g$etq <- factor(g$etq, levels = unique(g$etq[order(g$grupo)]))
cols <- c("Sin cola" = "#C5CFD8", "Cola superior" = calor, "Cola inferior" = azul, "Ambas" = "#6B4E9B")
p23 <- ggplot(g, aes(etq, p, fill = cola)) + geom_col(width = .7, colour = "white") +
  geom_text(aes(label = ifelse(p >= 5, paste0(round(p), " %"), "")), position = position_stack(vjust = .5), size = 4.4, colour = "black") +
  scale_fill_manual(values = cols, name = NULL, breaks = lab) + labs(x = NULL, y = "% de ajustes", title = "Cola de la cópula ganadora", subtitle = "Par intensidad-duración, todas las estaciones") + tema
ggsave(file.path(out, "cola_por_grupo.png"), p23 + theme(axis.text.x = element_text(size = 11)), width = 8.8, height = 5.4, dpi = 130)

# 24: mapa
rg <- read.table(REGIONES_CSV, sep = ";", header = TRUE, dec = ".", quote = "", fileEncoding = "latin1", stringsAsFactors = FALSE)
co <- data.frame(station_id = as.character(rg[[1]]), lon = rg[[4]], lat = rg[[5]])
co <- bind_rows(co, data.frame(station_id = c("87345", "87420", "87569", "87934"), lon = c(-64.198, -68.873, -58.731, -67.759), lat = c(-31.421, -32.895, -34.555, -53.780)))
s <- r %>% group_by(station_id) %>% summarise(k = n(), cc = sum(con_cola), p = cc / k, .groups = "drop") %>% left_join(co, by = "station_id") %>%
  mutate(grupo = factor(ifelse(cc == 0, "Ningún tipo con cola", "Al menos un tipo con cola"), levels = c("Ningún tipo con cola", "Al menos un tipo con cola")))
cat("sin coord:", sum(is.na(s$lon)), " | estaciones:", nrow(s), "\n"); print(table(s$grupo))
arg <- ggplot2::map_data("world", region = "Argentina")
p24 <- ggplot() + geom_polygon(data = arg, aes(long, lat, group = group), fill = "#EEF2F5", colour = "grey60", linewidth = .3) +
  geom_point(data = s, aes(lon, lat, colour = grupo), size = 3.8, alpha = .95) +
  scale_colour_manual(values = c("Ningún tipo con cola" = azul, "Al menos un tipo con cola" = calor), name = NULL) +
  coord_quickmap(xlim = c(-74, -53), ylim = c(-55.5, -21.5)) +
  labs(x = NULL, y = NULL, title = "Estaciones por cola de la cópula", subtitle = "Según los tipos de evento con cópula ganadora") +
  tema + theme(panel.grid = element_blank(), axis.text = element_blank(), legend.position = "bottom")
ggsave(file.path(out, "mapa_estaciones.png"), p24, width = 5.4, height = 7.0, dpi = 130)

# 25: de mayor a menor numero de eventos
h <- r %>% mutate(grupo = case_when(clase == "Ola" & n < 30 ~ "Olas, menos de 30", clase == "Ola" ~ "Olas, 30 o más", TRUE ~ "Períodos")) %>%
  group_by(grupo) %>% summarise(aj = n(), p = 100 * mean(con_cola), nm = round(median(n)), .groups = "drop")
h$grupo <- factor(h$grupo, levels = c("Períodos", "Olas, 30 o más", "Olas, menos de 30"))
h$etq <- factor(paste0(h$grupo, "\n", h$nm, " eventos"), levels = paste0(levels(h$grupo), "\n", h$nm[order(h$grupo)], " eventos"))
print(as.data.frame(h))
p25 <- ggplot(h, aes(etq, p, fill = grupo)) + geom_col(width = .6) +
  geom_text(aes(label = paste0(round(p), " %")), vjust = -.4, size = 5.6, fontface = "bold") +
  scale_fill_manual(values = c("Períodos" = azul, "Olas, 30 o más" = "#D98A66", "Olas, menos de 30" = calor), guide = "none") +
  coord_cartesian(ylim = c(0, 75)) +
  labs(x = NULL, y = "% de ajustes con cola", title = "Cópulas ganadoras con cola", subtitle = sprintf("%d ajustes con cópula ganadora", sum(h$aj))) + tema
ggsave(file.path(out, "n_chico.png"), p25, width = 6.2, height = 5.0, dpi = 130)
