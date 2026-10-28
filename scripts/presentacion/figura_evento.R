# Figura: isolineas de 87585 (olas de calor) con la ola de diciembre de 2013 marcada
# Se corre desde la raiz del repo: Rscript scripts/presentacion/figura_evento.R
# Entradas: data/output/*.rds (corrida Olas_20), data/input/eventos_identificados_frio_calor.csv
# Salidas: DIR_PRES (por defecto data/output/presentacion); las figuras van en DIR_PRES/figuras
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
dir.create(file.path(DIR_PRES, "figuras"), recursive = TRUE, showWarnings = FALSE)
REGIONES_CSV <- Sys.getenv("REGIONES_CSV", "data/input/regiones.csv")
suppressMessages({library(dplyr); library(ggplot2); library(copula)})
source("lib/funciones_ajuste_distribuciones.R", encoding="UTF-8")
x <- readRDS("data/output/copulas_olacalor.rds") %>% filter(station_id=="87585", variable_x=="intensidad", variable_y=="duracion", tipo_serie=="observada")
mv <- x$ajuste_multivariado[[1]]
ev <- readr::read_csv("data/input/eventos_identificados_frio_calor.csv", show_col_types=FALSE) %>% filter(station_id=="87585", tipo_evento=="OlaCalorTXyTN", duracion>=3)
f <- as.Date(ev$fecha_inicio); n <- length(f); N <- as.numeric(diff(range(f)))/365.25
FI <- function(q) do.call("pglo", c(list(q=q), mv@paramMargins[[1]])); FD <- function(q) do.call("ppe3", c(list(q=q), mv@paramMargins[[2]]))
g <- expand.grid(d = seq(2.9, 11.5, length.out = 220), i = seq(0.8, 5.8, length.out = 220))
u <- FI(g$i); v <- FD(g$d); g$T <- N/(n*pmax(1-u-v+pCopula(cbind(u,v), mv@copula), 1e-6))
ev$ev <- ev$fecha_inicio == as.Date("2013-12-22")
p <- ggplot(g, aes(d, i)) + geom_raster(aes(fill = T)) +
  scale_fill_viridis_c(trans = "log10", limits = c(2, 500), oob = scales::squish, breaks = c(2,5,10,20,50,100,200,500), name = "Período de\nretorno (años)") +
  geom_contour(aes(z = T), breaks = c(2,5,10,20,50,100,200,500), colour = "white", linewidth = .5) +
  geom_point(data = ev %>% filter(!ev), aes(duracion, intensidad), shape = 21, fill = "white", size = 2.4) +
  geom_point(data = ev %>% filter(ev), aes(duracion, intensidad), shape = 8, size = 6, stroke = 1.5, colour = "#B4491F") +
  annotate("label", x = 9, y = 4.45, label = "Dic. 2013: 9 días, 3,6\nT combinado ≈ 138 años", size = 4.6, fontface = "bold", colour = "#B4491F") +
  labs(x = "Duración (días)", y = "Intensidad", title = "Período de retorno combinado, olas de calor, Buenos Aires") +
  theme_minimal(base_size = 16) + theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 16), plot.background = element_rect(fill = "white", colour = NA)) +
  coord_cartesian(expand = FALSE)
ggsave(file.path(DIR_PRES, "figuras", "isolineas_evento.png"), p, width = 8.1, height = 5.94, dpi = 140)
