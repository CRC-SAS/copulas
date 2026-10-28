# Figuras y tablas nuevas tras la segunda ronda de comentarios de Natalia (octubre de 2026)
# Se corre desde la raiz del repo: Rscript scripts/presentacion/figuras_comentarios_nati_2.R
# Entradas: data/input/eventos_identificados_frio_calor.csv, data/input/enso_roni.txt, data/input/regiones.csv (Latin-1),
#   y la serie diaria de 87585 (variable SERIE_87585)
# Salidas: DIR_PRES/figuras (png) y DIR_PRES/tablas_enso (csv)
# La estacion 87420 se excluye de todas las tablas: Natalia la va a analizar por separado (tiene muchos mas eventos que el resto)
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
out <- file.path(DIR_PRES, "figuras"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DIR_PRES, "tablas_enso"), recursive = TRUE, showWarnings = FALSE)
SERIE_87585 <- Sys.getenv("SERIE_87585", "/data/devel/crc-sas/climindvis-webinar-2025/data/Registros_87585_Buenos_Aires.csv")
suppressMessages({library(dplyr); library(ggplot2); library(patchwork); library(lubridate)})
source("lib/funciones_enso.R", encoding = "UTF-8")
tema <- theme_minimal(base_size = 16) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 17),
  plot.background = element_rect(fill = "white", colour = NA))
col_calor <- "#B4491F"; col_frio <- "#1F6F9F"; azul <- "#1F5F8B"; gris <- "#56606B"
coma <- function(x, d = 1) format(round(x, d), decimal.mark = ",", nsmall = d)
EXCLUIDA <- "87420"

ev <- readr::read_csv("data/input/eventos_identificados_frio_calor.csv", show_col_types = FALSE) %>%
  mutate(station_id = as.character(station_id), fecha_inicio = as.Date(fecha_inicio)) %>% filter(duracion >= 3)
calor <- ev %>% filter(station_id == "87585", tipo_evento == "OlaCalorTXyTN")

# ---- Slide 3: la serie empieza el 18/12 para que la ola (22 al 30/12) quede centrada
d <- read.delim(SERIE_87585, stringsAsFactors = FALSE) %>%
  mutate(fecha = as.Date(fecha), tmax = suppressWarnings(as.numeric(tmax)), tmin = suppressWarnings(as.numeric(tmin)))
base <- d %>% filter(year(fecha) >= 1961, year(fecha) <= 2010, month(fecha) %in% c(10:12, 1:3))
u_tx <- unname(quantile(base$tmax, .9, na.rm = TRUE)); u_tn <- unname(quantile(base$tmin, .9, na.rm = TRUE))
w <- d %>% filter(fecha >= as.Date("2013-12-18"), fecha <= as.Date("2014-01-03"))
o <- w %>% filter(fecha >= as.Date("2013-12-22"), fecha <= as.Date("2013-12-30"))
p3 <- ggplot(w, aes(fecha)) +
  geom_ribbon(data = o, aes(ymin = u_tx, ymax = tmax), fill = col_calor, alpha = .35) +
  geom_ribbon(data = o, aes(ymin = u_tn, ymax = tmin), fill = "#7FC4E8", alpha = .55) +
  geom_hline(yintercept = c(u_tx, u_tn), linetype = "dashed", colour = c(col_calor, col_frio), linewidth = .6) +
  geom_line(aes(y = tmax), colour = col_calor, linewidth = 1) + geom_point(aes(y = tmax), colour = col_calor, size = 1.8) +
  geom_line(aes(y = tmin), colour = col_frio, linewidth = 1) + geom_point(aes(y = tmin), colour = col_frio, size = 1.8) +
  annotate("text", x = as.Date("2013-12-26"), y = 14.5, size = 4.8, fontface = "bold", colour = "grey20", label = "Ola: 22 al 30 de diciembre, 9 días") +
  annotate("text", x = as.Date("2013-12-19"), y = 38.6, size = 4.6, colour = col_calor, label = "TX (máxima)", fontface = "bold") +
  annotate("text", x = as.Date("2013-12-19"), y = 18.2, size = 4.6, colour = col_frio, label = "TN (mínima)", fontface = "bold") +
  scale_x_date(breaks = seq(as.Date("2013-12-18"), as.Date("2014-01-03"), by = "2 days"), labels = function(x) format(x, "%d/%m")) +
  coord_cartesian(ylim = c(13.5, 40)) +
  labs(x = "Fecha (dd/mm)", y = "Temperatura (°C)", title = "Buenos Aires (87585), diciembre de 2013",
       subtitle = sprintf("Líneas punteadas: umbrales de TX (%s °C) y de TN (%s °C),\npercentil 90 de octubre a marzo, 1961-2010. Sombreado: exceso sobre el umbral", coma(u_tx), coma(u_tn))) + tema
ggsave(file.path(out, "umbrales_ola.png"), p3, width = 10.5, height = 6.3, dpi = 130)

# ---- Slide 9 (nueva): de los datos a los rangos, con barras (la duracion es entera)
n <- nrow(calor)
tb <- calor %>% count(duracion) %>% arrange(duracion) %>% mutate(nn = n, fin = cumsum(nn) / sum(nn), ini = fin - nn / sum(nn),
  rango = (cumsum(nn) - (nn - 1) / 2) / (sum(nn) + 1))
cat("rangos por duracion:\n"); print(as.data.frame(tb))
b1 <- ggplot(tb, aes(factor(duracion), nn)) + geom_col(fill = col_calor, alpha = .8, width = .7) +
  geom_text(aes(label = nn), vjust = -.4, size = 5.2, fontface = "bold") + scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "1. Cuántas olas hay de cada duración", x = "Duración (días)", y = "Número de olas") + tema
tb$etq <- ifelse(tb$duracion <= 5, sprintf("%d días: %d olas\nrango %s", tb$duracion, tb$nn, coma(tb$rango, 2)), NA)
etq_alta <- sprintf("6 a 9 días: %d olas\nrangos de %s a %s", sum(tb$nn[tb$duracion >= 6]), coma(min(tb$rango[tb$duracion >= 6]), 2), coma(max(tb$rango[tb$duracion >= 6]), 2))
b2 <- ggplot(tb) + geom_rect(aes(xmin = 0, xmax = 1, ymin = ini, ymax = fin, fill = factor(duracion)), colour = "white", linewidth = 1) +
  geom_point(aes(x = .5, y = rango), shape = 21, fill = "white", colour = "black", size = 3.5) +
  geom_text(data = tb[tb$duracion <= 5, ], aes(x = 1.08, y = (ini + fin) / 2, label = etq), hjust = 0, size = 4.2, lineheight = .9) +
  annotate("text", x = 1.08, y = .94, label = etq_alta, hjust = 0, size = 4.2, lineheight = .9) +
  scale_fill_manual(values = colorRampPalette(c("#F2C9B5", col_calor, "#6B2A10"))(nrow(tb)), guide = "none") +
  scale_x_continuous(limits = c(0, 2.2)) + scale_y_continuous(breaks = seq(0, 1, .25)) +
  labs(title = "2. Ordenadas de menor a mayor: el rango", x = NULL, y = "Rango de duración (0 a 1)") +
  tema + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
ggsave(file.path(out, "duracion_barras_rangos.png"), b1 + b2 + plot_layout(widths = c(1, 1.15)), width = 13, height = 5.3, dpi = 130)

# ---- Slide 7: maximos por tipo de evento, sin la estacion excluida y con nombres
rg <- read.csv2("data/input/regiones.csv", fileEncoding = "latin1")
rg <- data.frame(station_id = as.character(rg[[1]]), loc = trimws(rg[[2]]), reg = rg[[6]])
rg <- rbind(rg, data.frame(station_id = c("87345", "87420", "87569", "87934"), loc = c("Córdoba Obs.", "Mendoza Obs.", "San Miguel", "Río Grande Aero"), reg = c("CENTRO", "CUYO", "BUENOS AIRES", "PATAGONIA")))
rg <- rg[!duplicated(rg$station_id), ]
cnt <- ev %>% count(station_id, tipo_evento) %>% filter(n >= 20, station_id != EXCLUIDA) %>% left_join(rg, by = "station_id")
cat("\nextremos sin", EXCLUIDA, ":\n")
print(as.data.frame(cnt %>% group_by(tipo_evento) %>% summarise(est = n(), mediana = median(n), max = max(n), quien = paste(loc[n == max(n)], collapse = " y "))))

# ---- Slides nuevas: porcentaje de eventos por fase ENSO (al inicio del evento), por region y por estacion
cl <- ClasificarMesesEnso(LeerEventosEnso("data/input/enso_roni.txt"))
colapsar <- function(f) factor(ifelse(startsWith(f, "Niño"), "Niño", ifelse(startsWith(f, "Niña"), "Niña", f)), levels = c("Niña", "Neutro", "Niño"))
cl$fase <- colapsar(cl$evento_enso); cl$anio <- (cl$ym - 1) %/% 12; cl$mes <- (cl$ym - 1) %% 12 + 1
ev$fase <- colapsar(as.character(AsignarEventoEnso(ev$fecha_inicio, cl)))
ee <- ev %>% filter(!is.na(fase), station_id != EXCLUIDA) %>% group_by(station_id, tipo_evento) %>% filter(n() >= 20) %>% ungroup() %>% left_join(rg, by = "station_id")
stopifnot(!anyNA(ee$reg))
TIPOS <- c(OlaCalorTXyTN = "Olas de calor", OlaFrioTXyTN = "Olas de frío", PeriodoCalidoTX = "Períodos cálidos de TX",
           PeriodoCalidoTN = "Períodos cálidos de TN", PeriodoFrioTX = "Períodos fríos de TX", PeriodoFrioTN = "Períodos fríos de TN")
# referencia: % de meses de cada fase entre 1961 y 2025 en los meses del semestre del evento (calor: oct-mar, frio: abr-sep)
ref <- function(tp) { m <- if (grepl("Calor|Calido", tp)) c(10:12, 1:3) else 4:9
  b <- cl %>% filter(anio >= 1961, anio <= 2025, mes %in% m); 100 * as.numeric(prop.table(table(b$fase))) }
refs <- t(sapply(names(TIPOS), ref)); colnames(refs) <- c("Niña", "Neutro", "Niño"); print(round(refs, 1))
pct <- function(x) 100 * c(mean(x == "Niña"), mean(x == "Neutro"), mean(x == "Niño"))
tab_reg <- ee %>% group_by(tipo_evento, reg) %>% summarise(est = n_distinct(station_id), n = n(), Nina = pct(fase)[1], Neutro = pct(fase)[2], Nino = pct(fase)[3], .groups = "drop")
tab_est <- ee %>% group_by(tipo_evento, station_id, loc, reg) %>% summarise(n = n(), Nina = pct(fase)[1], Neutro = pct(fase)[2], Nino = pct(fase)[3], .groups = "drop")
write.csv(tab_reg, file.path(DIR_PRES, "tablas_enso", "enso_por_region.csv"), row.names = FALSE)
write.csv(tab_est, file.path(DIR_PRES, "tablas_enso", "enso_por_estacion.csv"), row.names = FALSE)
cat("\nolas de calor por region:\n"); print(as.data.frame(tab_reg %>% filter(tipo_evento == "OlaCalorTXyTN")), digits = 3)

orden_reg <- c("NOA", "NORTE", "NEA", "CENTRO", "CUYO", "BUENOS AIRES", "PATAGONIA")
abrev <- c(NOA = "NOA", NORTE = "NORTE", NEA = "NEA", CENTRO = "CENTRO", CUYO = "CUYO", `BUENOS AIRES` = "BA", PATAGONIA = "PAT")
# tabla simple: % de eventos que empezaron en cada fase (en negrita la fase mas frecuente de cada fila)
dibujar <- function(df, tp, etiqueta, sz = 5, pie = NULL) {
  # df: etq, n, Nina, Neutro, Nino (ya ordenado). Devuelve un ggplot con la tabla
  f <- data.frame(y = rep(rev(seq_len(nrow(df))), 3), x = rep(1:3, each = nrow(df)), v = c(df$Nina, df$Neutro, df$Nino))
  f$mx <- ave(f$v, f$y, FUN = function(z) z == max(z))
  ggplot() +
    geom_text(data = f, aes(x + 3.0, y, label = sprintf("%.0f", v), fontface = ifelse(mx == 1, "bold", "plain")), size = sz) +
    geom_text(data = data.frame(y = rev(seq_len(nrow(df))), l = df$etq, n = df$n), aes(.2, y, label = l), hjust = 0, size = sz) +
    geom_text(data = data.frame(y = rev(seq_len(nrow(df))), n = df$n), aes(7.1, y, label = n), hjust = 1, size = sz, colour = gris) +
    annotate("text", x = c(4.0, 5.0, 6.0), y = nrow(df) + 1, label = c("% Niña", "% Neut.", "% Niño"), fontface = "bold", size = sz) +
    annotate("text", x = .2, y = nrow(df) + 1, label = etiqueta, hjust = 0, fontface = "bold", size = sz) +
    annotate("text", x = 7.1, y = nrow(df) + 1, label = "N", hjust = 1, fontface = "bold", size = sz, colour = gris) +
    annotate("text", x = .2, y = -.1, hjust = 0, size = sz, fontface = "italic", colour = gris,
             label = "% de eventos por fase de inicio (N: eventos)") +
    scale_x_continuous(limits = c(0, 7.15)) + scale_y_continuous(limits = c(-1.1, nrow(df) + 1.6)) + theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = NA))
}
# a) por region: dos tipos por imagen
reg_df <- function(tp) { z <- tab_reg %>% filter(tipo_evento == tp) %>% mutate(reg = factor(reg, orden_reg)) %>% arrange(reg)
  tot <- ee %>% filter(tipo_evento == tp)
  z <- z %>% mutate(etq = sprintf("%s (%d est.)", reg, est))
  tt <- data.frame(reg = "Total", est = n_distinct(tot$station_id), n = nrow(tot), Nina = pct(tot$fase)[1], Neutro = pct(tot$fase)[2], Nino = pct(tot$fase)[3])
  tt$etq <- sprintf("Total (%d est.)", tt$est); bind_rows(z %>% select(etq, n, Nina, Neutro, Nino), tt %>% select(etq, n, Nina, Neutro, Nino)) }
grupos <- list(regional_calor_frio = c("OlaCalorTXyTN", "OlaFrioTXyTN"), regional_calido = c("PeriodoCalidoTX", "PeriodoCalidoTN"),
               regional_frio = c("PeriodoFrioTX", "PeriodoFrioTN"))
for (g in names(grupos)) {
  tp <- grupos[[g]]
  pl <- lapply(tp, function(t) dibujar(reg_df(t), t, TIPOS[[t]], sz = 6))
  ggsave(file.path(out, sprintf("enso_%s.png", g)), pl[[1]] + pl[[2]], width = 16.64, height = 5.6, dpi = 100)
}
# b) por estacion: un png por tipo, en tres bloques
for (t in names(TIPOS)) {
  z <- tab_est %>% filter(tipo_evento == t) %>% mutate(reg = factor(reg, orden_reg), loc = ifelse(is.na(loc), station_id, sub("^Pcia\\. ", "", sub(" Aero$", "", loc)))) %>%
    arrange(reg, loc) %>% mutate(etq = sprintf("%s · %s", abrev[as.character(reg)], loc))
  k <- ceiling(nrow(z) / 3); idx <- split(seq_len(nrow(z)), ceiling(seq_len(nrow(z)) / k))
  pl <- lapply(idx, function(i) dibujar(z[i, ], t, "Región · estación", sz = 4.0))
  while (length(pl) < 3) pl[[length(pl) + 1]] <- plot_spacer()
  ggsave(file.path(out, sprintf("enso_estacion_%s.png", t)), wrap_plots(pl, ncol = 3), width = 16.64, height = 6.5, dpi = 100)
}
cat("listo\n")
