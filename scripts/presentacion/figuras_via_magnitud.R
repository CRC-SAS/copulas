# Figuras del periodo de retorno via duracion-magnitud: isolineas de las olas de calor de Buenos Aires (87585, sin copula intensidad-duracion)
# y validacion (acuerdo con la copula intensidad-duracion y calibracion contra los eventos observados).
# Se corre desde la raiz del repo, despues de periodo_retorno_via_magnitud.R. Entradas: DIR_PRES/via_magnitud_validacion.csv, data/output/copulas_olacalor.rds
# Salidas: DIR_PRES/figuras/via_magnitud_isolineas.png, via_magnitud_validacion.png y DIR_PRES/via_magnitud_caso.txt
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
out <- file.path(DIR_PRES, "figuras")
suppressMessages({library(dplyr); library(ggplot2); library(patchwork); library(copula)})
source("lib/funciones_ajuste_distribuciones.R", encoding = "UTF-8")
src <- readLines("lib/funciones_periodo_retorno.R", encoding = "latin1")
ini <- grep("^ProbabilidadConjuntaViaMagnitud <- function", src); eval(parse(text = src[ini:length(src)]))
azul <- "#1F5F8B"; calor <- "#B4491F"; gris <- "#56606B"
tema <- theme_minimal(base_size = 16) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 16),
  plot.subtitle = element_text(size = 12, colour = gris), plot.title.position = "plot", plot.background = element_rect(fill = "white", colour = NA))
breaks_dias <- function(l) seq(ceiling(l[1]), floor(l[2]), by = 1)
# Etiquetas de las curvas de nivel: punto medio de cada linea (igual criterio que el grafico del pipeline)
etiquetas_iso <- function(g, x, y, niveles) {
  l <- grDevices::contourLines(sort(unique(g[[x]])), sort(unique(g[[y]])), matrix(g$T[order(g[[y]], g[[x]])], nrow = length(unique(g[[x]]))), levels = niveles)
  bind_rows(lapply(l, function(z) { m <- ceiling(length(z$x) * .3); data.frame(d = z$x[m], i = z$y[m], nivel = z$level) }))
}
sink(file.path(DIR_PRES, "via_magnitud_caso.txt"), split = TRUE)

# ---- 1. isolineas de las olas de calor de Buenos Aires via magnitud
ev <- readr::read_csv("data/input/eventos_identificados_frio_calor.csv", show_col_types = FALSE) %>%
  filter(station_id == "87585", tipo_evento == "OlaCalorTXyTN", duracion >= 3) %>% mutate(fecha_inicio = as.Date(fecha_inicio))
n <- nrow(ev); N <- as.numeric(diff(range(ev$fecha_inicio))) / 365.25
# Tasa de eventos del periodo de retorno (solo para la presentacion): si el conteo anual tiene tendencia (Mann-Kendall p < 0,01) se usa la de 1991-2020 (nT eventos en NT anios)
anio <- as.integer(format(ev$fecha_inicio, "%Y")); p_tend <- Kendall::MannKendall(tabulate(anio - 1960, nbins = 65))$sl
if (p_tend < .01) { NT <- 30; nT <- sum(anio >= 1991 & anio <= 2020) } else { NT <- N; nT <- n }
cat(sprintf("Tendencia en la frecuencia: p = %.5f; tasa usada %.3f por anio (%.0f eventos en %.0f anios)\n", p_tend, nT / NT, nT, NT))
x <- readRDS("data/output/copulas_olacalor.rds") %>% filter(station_id == "87585", variable_x == "duracion", variable_y == "magnitud", tipo_serie == "observada")
mv <- x$ajuste_multivariado[[1]]
cat(sprintf("Olas de calor 87585: n=%d, N=%.1f anios; copula %s, marginales %s (duracion) y %s (magnitud)\n", n, N, x$mejor_copula[[1]]@class, x$mejor_ajuste_x_dist, x$mejor_ajuste_y_dist))
g <- expand.grid(d = seq(3, 9.5, length.out = 40), i = seq(.8, 4.6, length.out = 40))
g$T <- PeriodoRetornoViaMagnitud(mv, g$d - .5, g$i, NT, nT, nsim = 400000)   # x = k  <->  D >= k  <->  D > k - 0.5
niveles <- c(2, 5, 10, 20, 50, 100, 200, 500)
dic13 <- ev %>% filter(fecha_inicio == as.Date("2013-12-22"))
Tdic <- PeriodoRetornoViaMagnitud(mv, dic13$duracion - .5, dic13$intensidad, NT, nT, nsim = 2000000)
cat(sprintf("Ola del 22/12/2013 (%d dias, intensidad %.2f): T via magnitud = %.0f anios\n", dic13$duracion, dic13$intensidad, Tdic))
ej <- data.frame(d = c(5, 6, 7), i = c(2, 2.5, 2.5)); ej$T <- PeriodoRetornoViaMagnitud(mv, ej$d - .5, ej$i, NT, nT, nsim = 2000000)
ej$obs <- mapply(function(a, b) sum(ev$duracion >= a & ev$intensidad > b), ej$d, ej$i); ej$esp <- n / ((nT / NT) * ej$T)   # eventos esperados en el registro: n * P, con P = 1 / (tasa * T)
print(ej)
p1 <- ggplot(g, aes(d, i)) + geom_raster(aes(fill = T)) +
  scale_fill_viridis_c(trans = "log10", limits = c(2, 500), oob = scales::squish, breaks = niveles, name = "Período de\nretorno (años)") +
  geom_contour(aes(z = T), breaks = niveles, colour = "white", linewidth = .5) +
  geom_label(data = etiquetas_iso(g, "d", "i", niveles), aes(d, i, label = nivel), size = 4.2, linewidth = 0, label.padding = unit(.14, "lines"), fill = rgb(1, 1, 1, .8), colour = "grey20") +
  geom_point(data = ev, aes(duracion, intensidad), shape = 21, fill = "white", colour = "grey20", size = 2.4, stroke = .4) +
  geom_point(data = dic13, aes(duracion, intensidad), shape = 8, size = 6, stroke = 1.5, colour = "#E0501F") +
  annotate("label", x = dic13$duracion - 1.3, y = dic13$intensidad - 1.0, label = sprintf("Dic. 2013\nT ≈ %.0f años", Tdic), size = 4.4, fontface = "bold", colour = "#B4491F") +
  scale_x_continuous(breaks = breaks_dias) + coord_cartesian(expand = FALSE) +
  labs(x = "Duración (días): ola de al menos esos días", y = "Intensidad", title = "Período de retorno combinado, vía duración–magnitud",
       subtitle = "Olas de calor, Buenos Aires (87585): sin cópula intensidad–duración") + tema + theme(panel.grid = element_blank())
ggsave(file.path(out, "via_magnitud_isolineas.png"), p1, width = 8.1, height = 5.94, dpi = 140)

# ---- 2. validacion
r <- readr::read_csv(file.path(DIR_PRES, "via_magnitud_validacion.csv"), show_col_types = FALSE) %>% filter(corr == "d + 0,5") %>% mutate(E_dm = n * p_dm, E_di = n * p_di, E_ind = n * p_ind_emp)
b <- r %>% filter(ambas, E_dm >= 3, E_di >= 3)
pa <- ggplot(b, aes(1 / p_di, 1 / p_dm)) + geom_abline(slope = 1, intercept = 0, colour = "grey50") + geom_point(alpha = .25, size = 1.2, colour = azul) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Con la cópula intensidad–duración", y = "Vía duración–magnitud", title = "Los dos cálculos coinciden",
       subtitle = sprintf("1 / P en %d umbrales, %d estaciones-tipo\ncon las dos cópulas (mediana del cociente: %.2f)", nrow(b), n_distinct(paste(b$slug, b$station_id)), median(b$p_di / b$p_dm))) + tema
cal <- r %>% mutate(tramo = cut(E_dm, c(1, 3, 10, 30, Inf), right = FALSE, labels = c("1 a 3", "3 a 10", "10 a 30", "más de 30"))) %>% filter(!is.na(tramo)) %>% group_by(tramo) %>%
  summarise(`Vía duración–magnitud` = sum(obs) / sum(E_dm), `Independencia` = sum(obs) / sum(E_ind), .groups = "drop") %>% tidyr::pivot_longer(-tramo, names_to = "metodo", values_to = "cociente")
print(as.data.frame(cal))
pb <- ggplot(cal, aes(tramo, cociente, fill = metodo)) + geom_hline(yintercept = 1, colour = "grey40") + geom_col(position = position_dodge(.75), width = .7) +
  geom_text(aes(label = sprintf("%.2f", cociente)), position = position_dodge(.75), vjust = -.4, size = 4.2) +
  scale_fill_manual(values = c("Vía duración–magnitud" = azul, "Independencia" = "grey70"), name = NULL) + coord_cartesian(ylim = c(0, 1.7)) +
  labs(x = "Eventos esperados por el modelo", y = "Observados / esperados", title = "Calibración con los eventos observados", subtitle = "1 = el modelo acierta la frecuencia\n321 estaciones-tipo con cópula duración–magnitud") + tema + theme(legend.position = "top")
ggsave(file.path(out, "via_magnitud_validacion.png"), pa + pb, width = 13, height = 5.4, dpi = 130)
sink()
