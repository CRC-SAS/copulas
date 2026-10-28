# Figuras nuevas y corregidas tras los comentarios de Natalia sobre el deck (octubre de 2026)
# Se corre desde la raiz del repo: Rscript scripts/presentacion/figuras_comentarios_nati.R
# Entradas: data/output/*.rds (corrida Olas_20), data/input/eventos_identificados_frio_calor.csv, data/input/enso_roni.txt,
#   data/input/regiones.csv y la serie diaria de 87585 (variable SERIE_87585)
# Requiere resumen_ajustes.R (para el mapa de estaciones)
# Salidas: DIR_PRES (por defecto data/output/presentacion); las figuras van en DIR_PRES/figuras
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
dir.create(file.path(DIR_PRES, "figuras"), recursive = TRUE, showWarnings = FALSE)
REGIONES_CSV <- Sys.getenv("REGIONES_CSV", "data/input/regiones.csv")
SERIE_87585 <- Sys.getenv("SERIE_87585", "/data/devel/crc-sas/climindvis-webinar-2025/data/Registros_87585_Buenos_Aires.csv")
suppressMessages({library(dplyr); library(ggplot2); library(patchwork); library(copula); library(lubridate)})
source("lib/funciones_enso.R", encoding = "UTF-8")
source("lib/funciones_ajuste_distribuciones.R", encoding = "UTF-8")
txt <- readLines("lib/funciones_periodo_retorno.R", encoding = "UTF-8")
# El eval solo evalua la definicion de ConstantesEnso tomada de lib/funciones_periodo_retorno.R del propio repo (mismo patron que los otros scripts)
i <- grep("^ConstantesEnso <- function", txt); j <- i + which(txt[(i+1):length(txt)] == "}")[1]
eval(parse(text = txt[i:j])); enso <- ConstantesEnso()
out <- file.path(DIR_PRES, "figuras")
tema <- theme_minimal(base_size = 16) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 17),
  plot.background = element_rect(fill = "white", colour = NA))
col_calor <- "#B4491F"; col_frio <- "#1F6F9F"; azul <- "#1F5F8B"; gris <- "#56606B"
coma <- function(x, d = 1) format(round(x, d), decimal.mark = ",", nsmall = d)

ev <- readr::read_csv("data/input/eventos_identificados_frio_calor.csv", show_col_types = FALSE) %>%
  mutate(station_id = as.character(station_id), fecha_inicio = as.Date(fecha_inicio)) %>% filter(duracion >= 3)
ba <- ev %>% filter(station_id == "87585")
calor <- ba %>% filter(tipo_evento == "OlaCalorTXyTN"); frio <- ba %>% filter(tipo_evento == "OlaFrioTXyTN")

# ---- Slide 3: serie diaria con umbrales (datos reales de 87585, ola del 22 al 30 de diciembre de 2013)
d <- read.delim(SERIE_87585, stringsAsFactors = FALSE) %>%
  mutate(fecha = as.Date(fecha), tmax = suppressWarnings(as.numeric(tmax)), tmin = suppressWarnings(as.numeric(tmin)))
base <- d %>% filter(year(fecha) >= 1961, year(fecha) <= 2010, month(fecha) %in% c(10:12, 1:3))
u_tx <- unname(quantile(base$tmax, .9, na.rm = TRUE)); u_tn <- unname(quantile(base$tmin, .9, na.rm = TRUE))
cat(sprintf("umbrales P90 oct-mar 1961-2010: TX %.2f, TN %.2f\n", u_tx, u_tn))
w <- d %>% filter(fecha >= as.Date("2013-12-12"), fecha <= as.Date("2014-01-03"))
o <- w %>% filter(fecha >= as.Date("2013-12-22"), fecha <= as.Date("2013-12-30"))
cat(sprintf("intensidad recalculada: %.3f\n", ((mean(o$tmax) - u_tx) + (mean(o$tmin) - u_tn)) / 2))
etq <- function(x) format(x, "%d/%m")
p3 <- ggplot(w, aes(fecha)) +
  geom_ribbon(data = o, aes(ymin = u_tx, ymax = tmax), fill = col_calor, alpha = .35) +
  geom_ribbon(data = o, aes(ymin = u_tn, ymax = tmin), fill = "#7FC4E8", alpha = .55) +
  geom_hline(yintercept = c(u_tx, u_tn), linetype = "dashed", colour = c(col_calor, col_frio), linewidth = .6) +
  geom_line(aes(y = tmax), colour = col_calor, linewidth = 1) + geom_point(aes(y = tmax), colour = col_calor, size = 1.8) +
  geom_line(aes(y = tmin), colour = col_frio, linewidth = 1) + geom_point(aes(y = tmin), colour = col_frio, size = 1.8) +
  annotate("text", x = as.Date("2013-12-26"), y = 14.5, size = 4.8, fontface = "bold", colour = "grey20", label = "Ola: 22 al 30 de diciembre, 9 días") +
  annotate("text", x = as.Date("2013-12-20"), y = 38.6, size = 4.6, colour = col_calor, label = "TX (máxima)", fontface = "bold") +
  annotate("text", x = as.Date("2013-12-17"), y = 18.2, size = 4.6, colour = col_frio, label = "TN (mínima)", fontface = "bold") +
  scale_x_date(breaks = seq(as.Date("2013-12-14"), as.Date("2014-01-03"), by = "4 days"), labels = etq) +
  coord_cartesian(ylim = c(13.5, 40)) +
  labs(x = "Fecha (dd/mm)", y = "Temperatura (°C)", title = "Buenos Aires (87585), diciembre de 2013",
       subtitle = sprintf("Líneas punteadas: umbrales de TX (%s °C) y de TN (%s °C),\npercentil 90 de octubre a marzo, 1961-2010. Sombreado: exceso sobre el umbral", coma(u_tx), coma(u_tn))) + tema
ggsave(file.path(out, "umbrales_ola.png"), p3, width = 10.5, height = 6.3, dpi = 130)

# ---- Slide 4: scatter sin jitter (la duracion es entera)
tau_c <- cor(calor$intensidad, calor$duracion, method = "kendall")
p5 <- ggplot(calor, aes(duracion, intensidad)) + geom_count(alpha = .7, colour = col_calor) +
  scale_size_area(max_size = 9, breaks = c(1, 5, 15, 30), name = "Olas") +
  scale_x_continuous(breaks = 3:9) +
  labs(title = "Olas de calor, Buenos Aires (87585)", x = "Duración (días)", y = "Intensidad",
       subtitle = sprintf("n = %d olas · tau de Kendall = %s · el tamaño indica cuántas olas coinciden", nrow(calor), coma(tau_c, 2))) + tema +
  theme(plot.subtitle = element_text(size = 11, colour = gris))
ggsave(file.path(out, "scatter_calor.png"), p5, width = 8.1, height = 5.94, dpi = 140)

# ---- Slide 7: de los datos a la copula (duracion en x, sin jitter)
u <- pobs(cbind(calor$duracion, calor$intensidad)); colnames(u) <- c("u", "v")
cat("rangos distintos de la duracion:", paste(round(sort(unique(u[, 1])), 3), collapse = " "), "\n")
th <- 1.165608
a <- ggplot(calor, aes(duracion, intensidad)) + geom_count(alpha = .75, colour = col_calor) + scale_size_area(max_size = 7, guide = "none") +
  scale_x_continuous(breaks = 3:9) + labs(title = "1. Los datos (73 olas)", x = "Duración (días)", y = "Intensidad") + tema
b <- ggplot(as.data.frame(u), aes(u, v)) + geom_count(alpha = .75, colour = azul) + scale_size_area(max_size = 7, guide = "none") +
  coord_fixed() + labs(title = "2. Sus rangos", x = "Rango de duración", y = "Rango de intensidad") + tema
set.seed(11); sm <- as.data.frame(rCopula(500, gumbelCopula(th))); names(sm) <- c("u", "v")
c3 <- ggplot(sm, aes(u, v)) + geom_point(size = 1.2, alpha = .4, colour = azul) +
  coord_fixed() + labs(title = "3. Cópula ajustada", x = "Rango de duración", y = "Rango de intensidad") + tema
ggsave(file.path(out, "de_datos_a_copula.png"), a + b + c3 + plot_layout(ncol = 3), width = 13, height = 4.6, dpi = 130)

# ---- Slide 20: scatter por fase ENSO (sin jitter)
cl0 <- ClasificarMesesEnso(LeerEventosEnso("data/input/enso_roni.txt"))
pl <- bind_rows(calor %>% mutate(tipo = "Olas de calor"), frio %>% mutate(tipo = "Olas de frío")) %>%
  mutate(evento_enso = factor(AsignarEventoEnso(fecha_inicio, cl0), levels = enso$niveles))
pen <- ggplot(pl, aes(duracion, intensidad, fill = evento_enso, shape = evento_enso)) +
  geom_point(size = 3.6, colour = "black", stroke = .4, alpha = .85) +
  scale_fill_manual(values = enso$colores, drop = FALSE, name = "Fase ENSO") + scale_shape_manual(values = enso$formas, drop = FALSE, name = "Fase ENSO") +
  scale_x_continuous(breaks = 3:9) + facet_wrap(~tipo, scales = "free") +
  labs(x = "Duración (días)", y = "Intensidad", subtitle = "Buenos Aires (87585), fase ENSO al inicio de cada ola") +
  tema + theme(strip.text = element_text(face = "bold", size = 16), legend.position = "right")
ggsave(file.path(out, "scatter_enso.png"), pen, width = 11.5, height = 5.4, dpi = 130)

# ---- Slide 8: familias, muestras con las esquinas marcadas, y curvas de nivel de la densidad
set.seed(7); n <- 600
sim <- function(cop, nom) { x <- rCopula(n, cop); data.frame(u = x[,1], v = x[,2], fam = nom) }
nf <- c("Normal: sin colas", "Gumbel: cola superior", "Clayton: cola inferior", "t: ambas colas")
dd <- bind_rows(sim(normalCopula(iTau(normalCopula(), .5)), nf[1]), sim(gumbelCopula(iTau(gumbelCopula(), .5)), nf[2]),
                sim(claytonCopula(iTau(claytonCopula(), .5)), nf[3]),
                sim(tCopula(iTau(tCopula(df = 3, df.fixed = TRUE), .5), df = 3, df.fixed = TRUE), nf[4]))
dd$fam <- factor(dd$fam, levels = nf)
cols <- setNames(c("#7C8A97", col_calor, azul, "#6B4E9B"), nf)
rect <- bind_rows(data.frame(fam = nf[2], x0 = .8, x1 = 1, y0 = .8, y1 = 1), data.frame(fam = nf[3], x0 = 0, x1 = .2, y0 = 0, y1 = .2),
                  data.frame(fam = nf[4], x0 = .8, x1 = 1, y0 = .8, y1 = 1), data.frame(fam = nf[4], x0 = 0, x1 = .2, y0 = 0, y1 = .2))
rect$fam <- factor(rect$fam, levels = nf)
p9 <- ggplot(dd, aes(u, v, colour = fam)) +
  geom_rect(data = rect, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fam), inherit.aes = FALSE, alpha = .2) +
  geom_point(size = .9, alpha = .5) + facet_wrap(~fam, nrow = 2) + coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = cols, guide = "none") + scale_fill_manual(values = cols, guide = "none") +
  labs(x = "Rango de la primera variable", y = "Rango de la segunda variable", subtitle = "Mismo tau de Kendall (0,5). Sombreado: la esquina donde hay cola") +
  tema + theme(strip.text = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 12, colour = gris))
ggsave(file.path(out, "familias_colas.png"), p9, width = 8.1, height = 7.4, dpi = 130)

tau <- .5; g <- expand.grid(x = seq(-2.6, 2.6, length.out = 160), y = seq(-2.6, 2.6, length.out = 160))
fams <- list("Normal" = normalCopula(iTau(normalCopula(), tau)),
             "t de Student" = tCopula(iTau(tCopula(df = 4, df.fixed = TRUE), tau), df = 4, df.fixed = TRUE),
             "Frank" = frankCopula(iTau(frankCopula(), tau)),
             "Clayton" = claytonCopula(iTau(claytonCopula(), tau)),
             "Joe" = joeCopula(iTau(joeCopula(), tau)), "Gumbel" = gumbelCopula(iTau(gumbelCopula(), tau)))
dens <- bind_rows(lapply(names(fams), function(k) {
  m <- mvdc(fams[[k]], c("norm", "norm"), list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))
  g$z <- dMvdc(as.matrix(g[, c("x", "y")]), m); g$fam <- k; g }))
dens$fam <- factor(dens$fam, levels = names(fams))
pd <- ggplot(dens, aes(x, y, z = z)) + geom_contour(colour = azul, bins = 10, linewidth = .6) +
  facet_wrap(~fam, nrow = 2) + coord_fixed() +
  labs(x = "Primera variable (normalizada)", y = "Segunda variable (normalizada)",
       subtitle = "Curvas de nivel de la densidad conjunta con marginales normales y el mismo tau de Kendall (0,5)") + tema +
  theme(strip.text = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 12, colour = gris))
ggsave(file.path(out, "familias_densidad.png"), pd, width = 11, height = 7.4, dpi = 130)

# ---- Slide 9: marginales (u y v) y copula ajustada con duracion en x
x_ <- readRDS("data/output/copulas_olacalor.rds") %>% filter(station_id == "87585", variable_x == "intensidad", variable_y == "duracion", tipo_serie == "observada")
mv <- x_$ajuste_multivariado[[1]]
FI <- function(q) do.call("pglo", c(list(q = q), mv@paramMargins[[1]]))
FD <- function(q) do.call("ppe3", c(list(q = q), mv@paramMargins[[2]]))
i0 <- 3; d0 <- 5; u0 <- FD(d0); v0 <- FI(i0); C0 <- pCopula(c(FI(i0), FD(d0)), mv@copula)
cat(sprintf("u (duracion <= 5) = %.4f, v (intensidad <= 3) = %.4f, C = %.4f, producto = %.4f\n", u0, v0, C0, u0 * v0))
xs <- seq(2.9, 9.6, length.out = 300)
pm1 <- ggplot() + stat_ecdf(data = calor, aes(duracion), geom = "step", colour = "grey55", linewidth = .8) +
  geom_line(data = data.frame(x = xs, F = FD(xs)), aes(x, F), colour = azul, linewidth = 1.1) +
  annotate("segment", x = d0, xend = d0, y = 0, yend = u0, linetype = "dashed", colour = "black") +
  annotate("segment", x = 2.9, xend = d0, y = u0, yend = u0, linetype = "dashed", colour = "black") +
  annotate("point", x = d0, y = u0, shape = 8, size = 4.5, stroke = 1.3) +
  annotate("label", x = 6.2, y = .55, label = sprintf("F(5 días) = %s", coma(u0, 2)), size = 4.6, fontface = "bold") +
  scale_x_continuous(breaks = 3:9) + coord_cartesian(xlim = c(2.9, 9.6), ylim = c(0, 1)) +
  labs(title = "Duración", x = "Duración (días)", y = "Probabilidad acumulada", subtitle = "Gris: datos. Azul: marginal ajustada (PE3)") + tema +
  theme(plot.subtitle = element_text(size = 11, colour = gris))
xi <- seq(.8, 5.8, length.out = 300)
pm2 <- ggplot() + stat_ecdf(data = calor, aes(intensidad), geom = "step", colour = "grey55", linewidth = .8) +
  geom_line(data = data.frame(x = xi, F = FI(xi)), aes(x, F), colour = col_calor, linewidth = 1.1) +
  annotate("segment", x = i0, xend = i0, y = 0, yend = v0, linetype = "dashed", colour = "black") +
  annotate("segment", x = .8, xend = i0, y = v0, yend = v0, linetype = "dashed", colour = "black") +
  annotate("point", x = i0, y = v0, shape = 8, size = 4.5, stroke = 1.3) +
  annotate("label", x = 4.3, y = .55, label = sprintf("F(intensidad 3) = %s", coma(v0, 2)), size = 4.6, fontface = "bold") +
  coord_cartesian(xlim = c(.8, 5.8), ylim = c(0, 1)) +
  labs(title = "Intensidad", x = "Intensidad", y = NULL, subtitle = "Gris: datos. Naranja: marginal ajustada (GLO)") + tema +
  theme(plot.subtitle = element_text(size = 11, colour = gris))
ggsave(file.path(out, "marginales_uv.png"), pm1 + pm2, width = 12, height = 5, dpi = 130)

pts <- data.frame(u = FD(calor$duracion), v = FI(calor$intensidad))
gr <- expand.grid(u = seq(.005, .995, length.out = 200), v = seq(.005, .995, length.out = 200))
gr$C <- pCopula(cbind(gr$v, gr$u), mv@copula); gr$ind <- gr$u * gr$v
br <- c(.05, .1, .2, .3, .5, .7)
pc <- ggplot() +
  geom_contour(data = gr, aes(u, v, z = ind), breaks = br, colour = "grey55", linetype = "dashed", linewidth = .6) +
  geom_contour(data = gr, aes(u, v, z = C), breaks = br, colour = azul, linewidth = 1) +
  geom_point(data = pts, aes(u, v), size = 2.2, alpha = .55, colour = col_calor) +
  annotate("point", x = u0, y = v0, shape = 8, size = 5, stroke = 1.4) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "F de la duración (u)", y = "F de la intensidad (v)", title = "Cópula ajustada: Gumbel, θ = 1,17",
       subtitle = "Continuas: cópula ajustada. Punteadas: independencia") + tema + theme(plot.subtitle = element_text(size = 11))
ggsave(file.path(out, "copula_ajustada.png"), pc, width = 8.1, height = 7.4, dpi = 130)

# ---- Slide 20 (nueva): estaciones donde la fase ENSO se nota mas
cl <- ClasificarMesesEnso(LeerEventosEnso("data/input/enso_roni.txt"))
colapsar <- function(f) factor(ifelse(startsWith(f, "Niño"), "Niño", ifelse(startsWith(f, "Niña"), "Niña", f)), levels = c("Niña", "Neutro", "Niño"))
olas <- ev %>% filter(grepl("^Ola", tipo_evento)) %>% mutate(fase = colapsar(as.character(AsignarEventoEnso(fecha_inicio, cl)))) %>% filter(!is.na(fase))
prueba <- olas %>% group_by(station_id, tipo_evento) %>% filter(n() >= 30, n_distinct(fase) == 3) %>%
  summarise(n = n(), p_int = kruskal.test(intensidad ~ fase)$p.value, p_dur = kruskal.test(duracion ~ fase)$p.value, .groups = "drop")
cat(sprintf("pruebas de Kruskal-Wallis: %d estaciones-tipo (%d pruebas); p < 0,05: intensidad %d, duracion %d (esperables por azar: %.1f)\n",
            nrow(prueba), 2 * nrow(prueba), sum(prueba$p_int < .05), sum(prueba$p_dur < .05), .05 * 2 * nrow(prueba)))
print(as.data.frame(prueba %>% arrange(pmin(p_int, p_dur)) %>% head(8)))
colf <- c("Niña" = unname(enso$colores["Niña"]), "Neutro" = unname(enso$colores["Neutro"]), "Niño" = unname(enso$colores["Niño"]))
caso <- function(st, tp, var, titulo, ylab) {
  z <- olas %>% filter(station_id == st, tipo_evento == tp)
  pv <- prueba %>% filter(station_id == st, tipo_evento == tp) %>% pull(if (var == "intensidad") p_int else p_dur)
  ggplot(z, aes(fase, .data[[var]], fill = fase)) + geom_boxplot(alpha = .35, outlier.shape = NA, width = .55, colour = "grey30") +
    geom_jitter(width = .12, height = 0, shape = 21, size = 2.8, colour = "black", stroke = .3) +
    scale_fill_manual(values = colf, guide = "none") +
    labs(title = titulo, subtitle = sprintf("%d olas · prueba de Kruskal-Wallis p = %s", nrow(z), coma(pv, 3)), x = "Fase ENSO al inicio de la ola", y = ylab) + tema +
    theme(plot.subtitle = element_text(size = 11, colour = gris))
}
pe <- caso("87393", "OlaCalorTXyTN", "duracion", "Monte Caseros: olas de calor", "Duración (días)") +
  caso("87934", "OlaFrioTXyTN", "intensidad", "Estación 87934: olas de frío", "Intensidad")
ggsave(file.path(out, "enso_estaciones.png"), pe, width = 12, height = 5.2, dpi = 130)

# ---- Slide 22 (nueva): mapa de las regiones
rg <- read.table(REGIONES_CSV, sep = ";", header = TRUE, dec = ".", quote = "", fileEncoding = "latin1", stringsAsFactors = FALSE)
co <- data.frame(station_id = as.character(rg[[1]]), lon = rg[[4]], lat = rg[[5]], region = toupper(rg[[6]]))
co <- bind_rows(co, data.frame(station_id = c("87345", "87420", "87569", "87934"), lon = c(-64.198, -68.873, -58.731, -67.759),
                               lat = c(-31.421, -32.895, -34.555, -53.780), region = c("CENTRO", "CUYO", "BUENOS AIRES", "PATAGONIA")))
co <- co %>% filter(station_id %in% unique(ev$station_id), station_id != "87420")  # 87420 excluida (decision de Natalia, 2026-10-23)
regcol <- c("NOA" = "#8C564B", "NORTE" = "#E39B2D", "NEA" = "#2CA02C", "CENTRO" = "#6B4E9B", "CUYO" = "#D6477F", "BUENOS AIRES" = "#1F5F8B", "PATAGONIA" = "#17A2B8")
nreg <- co %>% count(region); co$regetq <- factor(paste0(co$region, " (", nreg$n[match(co$region, nreg$region)], ")"), levels = paste0(names(regcol), " (", nreg$n[match(names(regcol), nreg$region)], ")"))
arg <- ggplot2::map_data("world", region = "Argentina")
pr <- ggplot() + geom_polygon(data = arg, aes(long, lat, group = group), fill = "#F1F4F7", colour = "grey60", linewidth = .3) +
  geom_point(data = co, aes(lon, lat, fill = regetq), shape = 21, colour = "grey20", size = 5, stroke = .4) +
  scale_fill_manual(values = setNames(unname(regcol), levels(co$regetq)), name = "Región (estaciones)") +
  coord_quickmap(xlim = c(-74, -53), ylim = c(-55.5, -21.5)) + labs(x = NULL, y = NULL) + tema +
  theme(panel.grid = element_blank(), axis.text = element_blank(), legend.position = c(.77, .22), legend.background = element_rect(fill = "white", colour = "grey80"),
        legend.title = element_text(size = 13, face = "bold"), legend.text = element_text(size = 13))
ggsave(file.path(out, "mapa_regiones.png"), pr, width = 6.6, height = 8.6, dpi = 130)

# ---- Slide 24: mapa de estaciones, con clases de proporcion y mas grande
r <- readr::read_csv(file.path(DIR_PRES, "resumen_ajustes.csv"), show_col_types = FALSE) %>%
  filter(variable_x == "intensidad", variable_y == "duracion", !is.na(familia)) %>% mutate(station_id = as.character(station_id), con_cola = cola != "Sin cola")
s <- r %>% group_by(station_id) %>% summarise(k = n(), cc = sum(con_cola), p = cc / k, .groups = "drop") %>% left_join(co, by = "station_id") %>%
  mutate(clase = cut(p, c(-.01, 0, .34, .67, 1), labels = c("Ninguno", "Hasta un tercio", "Entre un tercio y dos tercios", "Más de dos tercios")))
cat("estaciones del mapa:", nrow(s), "\n"); print(table(s$clase))
pm <- ggplot() + geom_polygon(data = arg, aes(long, lat, group = group), fill = "#EEF2F5", colour = "grey60", linewidth = .3) +
  geom_point(data = s, aes(lon, lat, fill = clase), shape = 21, colour = "grey20", size = 4.8, stroke = .4) +
  scale_fill_manual(values = c("Ninguno" = "#DCE6EE", "Hasta un tercio" = "#F4C7A8", "Entre un tercio y dos tercios" = "#E0804F", "Más de dos tercios" = "#9E3A14"),
                    name = "Tipos de evento de la\nestación con cópula con cola") +
  coord_quickmap(xlim = c(-74, -53), ylim = c(-55.5, -21.5)) + labs(x = NULL, y = NULL) + tema +
  theme(panel.grid = element_blank(), axis.text = element_blank(), legend.position = c(.72, .2), legend.background = element_rect(fill = "white", colour = "grey80"),
        legend.title = element_text(size = 11), legend.text = element_text(size = 11))
ggsave(file.path(out, "mapa_estaciones.png"), pm, width = 6.0, height = 8.6, dpi = 130)
