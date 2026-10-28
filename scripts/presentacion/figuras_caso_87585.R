# Figuras del caso representativo: Buenos Aires (87585), periodos calidos y frios de TX, par intensidad-duracion
# (con los filtros corregidos de dependencia y estacionariedad, 87585 solo tiene copula intensidad-duracion en estos dos tipos)
# Se corre desde la raiz del repo: Rscript scripts/presentacion/figuras_caso_87585.R
# Entradas: data/output/copulas_periodocalidotx.rds y copulas_periodofriotx.rds, data/input/eventos_identificados_frio_calor.csv, data/input/enso_roni.txt
# Salidas: DIR_PRES/figuras (png) y DIR_PRES/caso_87585.txt (numeros que se citan en las slides)
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
out <- file.path(DIR_PRES, "figuras"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
suppressMessages({library(dplyr); library(ggplot2); library(patchwork); library(copula); library(Kendall)})
source("lib/funciones_ajuste_distribuciones.R", encoding = "UTF-8")
source("lib/funciones_enso.R", encoding = "UTF-8")
txt <- readLines("lib/funciones_periodo_retorno.R", encoding = "UTF-8")
i <- grep("^ConstantesEnso <- function", txt); j <- i + which(txt[(i+1):length(txt)] == "}")[1]
eval(parse(text = txt[i:j])); enso <- ConstantesEnso()
tema <- theme_minimal(base_size = 16) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 17),
  plot.subtitle = element_text(size = 12, colour = "#56606B"), plot.background = element_rect(fill = "white", colour = NA))
col_calor <- "#B4491F"; col_frio <- "#1F6F9F"; azul <- "#1F5F8B"; gris <- "#56606B"
# Duraciones enteras: sin jitter en el eje x y marcas en dias enteros (cada 1 o cada 2 segun el rango)
breaks_dias <- function(l) seq(ceiling(l[1]), floor(l[2]), by = ifelse(diff(l) > 9, 2, 1))
coma <- function(x, d = 1) format(round(x, d), decimal.mark = ",", nsmall = d)
sink_file <- file.path(DIR_PRES, "caso_87585.txt"); sink(sink_file, split = TRUE)

ev_all <- readr::read_csv("data/input/eventos_identificados_frio_calor.csv", show_col_types = FALSE) %>%
  filter(station_id == "87585", duracion >= 3) %>% mutate(fecha_inicio = as.Date(fecha_inicio))

caso <- function(tp, slug, col) {
  e <- ev_all %>% filter(tipo_evento == tp)
  x <- readRDS(sprintf("data/output/copulas_%s.rds", slug)) %>%
    filter(station_id == "87585", variable_x == "intensidad", variable_y == "duracion", tipo_serie == "observada")
  mv <- x$ajuste_multivariado[[1]]
  fx <- get(paste0("p", x$mejor_ajuste_x_dist)); fy <- get(paste0("p", x$mejor_ajuste_y_dist))
  FI <- function(q) do.call(fx, c(list(q = q), mv@paramMargins[[1]])); FD <- function(q) do.call(fy, c(list(q = q), mv@paramMargins[[2]]))
  n <- nrow(e); N <- as.numeric(diff(range(e$fecha_inicio))) / 365.25
  # Tasa de eventos por anio para el periodo de retorno: n / N, salvo que el conteo anual tenga tendencia (Mann-Kendall p < 0,01):
  # en ese caso (solo para la presentacion) se usa la tasa de 1991-2020. P no depende de la tasa.
  anio <- as.integer(format(e$fecha_inicio, "%Y")); conteo <- tabulate(anio - 1960, nbins = 65)
  p_tend <- MannKendall(conteo)$sl; n9120 <- sum(anio >= 1991 & anio <= 2020)
  lam <- if (p_tend < .01) n9120 / 30 else n / N
  cat(sprintf("%s: tendencia en la frecuencia p = %.4f; tasa usada %.3f por anio (%s; n/N = %.3f, 1991-2020 = %.3f)\n", tp, p_tend, lam, ifelse(p_tend < .01, "1991-2020", "n/N"), n / N, n9120 / 30))
  tau <- cor(e$intensidad, e$duracion, method = "kendall"); rho <- mv@copula@parameters
  # Correccion de continuidad: la duracion es entera; "al menos d dias" se evalua en d - 0.5
  Tc <- function(d, i) { u <- FI(i); v <- FD(d - .5); 1 / (lam * pmax(1 - u - v + pCopula(cbind(u, v), mv@copula), 1e-9)) }
  cat(sprintf("\n== %s: n=%d, N=%.1f anios, tau=%.3f, %s rho=%.3f, marginales %s (intensidad) y %s (duracion)\n", tp, n, N, tau, x$mejor_copula, rho, x$mejor_ajuste_x_dist, x$mejor_ajuste_y_dist))
  list(e = e, x = x, mv = mv, FI = FI, FD = FD, n = n, N = N, lam = lam, tau = tau, rho = rho, Tc = Tc, col = col, tp = tp)
}
C <- caso("PeriodoCalidoTX", "periodocalidotx", col_calor)
F <- caso("PeriodoFrioTX", "periodofriotx", col_frio)

# ---- punto de ejemplo del texto: mas de 5 dias e intensidad mayor a 2 (duracion entera; el pipeline usa P(D > d))
ejemplo <- function(K, d0, i0) {
  u0 <- K$FD(d0 - .5); v0 <- K$FI(i0); C0 <- pCopula(c(u0, v0), K$mv@copula); P <- 1 - u0 - v0 + C0   # al menos d0 dias <-> D > d0 - 0.5
  obs <- sum(K$e$duracion >= d0 & K$e$intensidad > i0)
  cat(sprintf("%s: D>=%g (evaluada en %g), I>%g: u=%.3f v=%.3f C=%.3f ind=%.3f P=%.4f T=%.1f anios (indep %.1f); observados %d en %.0f anios (esperados %.1f)\n",
              K$tp, d0, d0 - .5, i0, u0, v0, C0, u0 * v0, P, 1 / (K$lam * P), 1 / (K$lam * (1 - u0) * (1 - v0)), obs, K$N, K$n * P))
  c(u0 = u0, v0 = v0, C0 = C0)
}
for (K in list(C, F)) for (p in list(c(5, 2), c(5, 3), c(7, 2), c(4, 2))) ejemplo(K, p[1], p[2])
D0 <- 5; I0 <- 2; pc0 <- ejemplo(C, D0, I0)

# ---- 1. scatter
set.seed(1)
sc <- function(K, ttl) ggplot(K$e, aes(duracion, intensidad)) + geom_point(size = 3, alpha = .5, colour = K$col) +
  scale_x_continuous(breaks = breaks_dias) +
  labs(title = ttl, x = "Duración (días)", y = "Intensidad", subtitle = sprintf("n = %d períodos · tau de Kendall = %s", K$n, coma(K$tau, 2))) + tema
ggsave(file.path(out, "caso_scatter_calor.png"), sc(C, "Períodos cálidos de TX, Buenos Aires (87585)"), width = 8.1, height = 5.94, dpi = 140)
ggsave(file.path(out, "caso_scatter_frio.png"), sc(F, "Períodos fríos de TX, Buenos Aires (87585)"), width = 8.1, height = 5.94, dpi = 140)

# ---- 2. de los datos a la copula
u <- pobs(cbind(C$e$duracion, C$e$intensidad)); colnames(u) <- c("u", "v")
a <- ggplot(C$e, aes(duracion, intensidad)) + geom_point(size = 2.2, alpha = .5, colour = col_calor) + scale_x_continuous(breaks = breaks_dias) +
  labs(title = "1. Los datos", x = "Duración", y = "Intensidad") + tema
b <- ggplot(as.data.frame(u), aes(u, v)) + geom_point(size = 2.2, alpha = .6, colour = azul) + coord_fixed() +
  labs(title = "2. Sus rangos", x = "Rango de duración", y = "Rango de intensidad") + tema
set.seed(11); sm <- as.data.frame(rCopula(500, C$mv@copula)); names(sm) <- c("u", "v")
c3 <- ggplot(sm, aes(u, v)) + geom_point(size = 1.4, alpha = .45, colour = azul) + coord_fixed() +
  labs(title = "3. La cópula ajustada", x = "Rango de duración", y = "Rango de intensidad") + tema
ggsave(file.path(out, "caso_de_datos_a_copula.png"), a + b + c3 + plot_layout(ncol = 3), width = 13, height = 4.6, dpi = 130)
tb <- C$e %>% count(duracion) %>% arrange(duracion) %>% mutate(fin = cumsum(n) / sum(n), ini = fin - n / sum(n), rango = (cumsum(n) - (n - 1) / 2) / (sum(n) + 1))
cat("\nrangos por duracion:\n"); print(as.data.frame(tb))

# ---- 3. barras de rangos
b1 <- ggplot(tb, aes(factor(duracion), n)) + geom_col(fill = col_calor, alpha = .8, width = .7) +
  geom_text(aes(label = n), vjust = -.4, size = 5.2, fontface = "bold") + scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "1. Períodos por duración", x = "Duración (días)", y = "Número de períodos") + tema
tb$etq <- ifelse(tb$duracion <= 5, sprintf("%d días: %d períodos\nrango %s", tb$duracion, tb$n, coma(tb$rango, 2)), NA)
etq_alta <- sprintf("6 a %d días: %d períodos\nrangos de %s a %s", max(tb$duracion), sum(tb$n[tb$duracion >= 6]), coma(min(tb$rango[tb$duracion >= 6]), 2), coma(max(tb$rango[tb$duracion >= 6]), 2))
b2 <- ggplot(tb) + geom_rect(aes(xmin = 0, xmax = 1, ymin = ini, ymax = fin, fill = factor(duracion)), colour = "white", linewidth = 1) +
  geom_point(aes(x = .5, y = rango), shape = 21, fill = "white", colour = "black", size = 3.5) +
  geom_text(data = tb[tb$duracion <= 5, ], aes(x = 1.08, y = (ini + fin) / 2, label = etq), hjust = 0, size = 4.2, lineheight = .9) +
  annotate("text", x = 1.08, y = .94, label = etq_alta, hjust = 0, size = 4.2, lineheight = .9) +
  scale_fill_manual(values = colorRampPalette(c("#F2C9B5", col_calor, "#6B2A10"))(nrow(tb)), guide = "none") +
  scale_x_continuous(limits = c(0, 2.2)) + scale_y_continuous(breaks = seq(0, 1, .25)) +
  labs(title = "2. Ordenados de menor a mayor: el rango", x = NULL, y = "Rango de duración (0 a 1)") +
  tema + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
ggsave(file.path(out, "caso_duracion_barras_rangos.png"), b1 + b2 + plot_layout(widths = c(1, 1.15)), width = 13, height = 5.3, dpi = 130)

# ---- 4. marginales y 5. copula ajustada
xs <- seq(2.9, max(C$e$duracion) + .5, length.out = 300); xi <- seq(min(C$e$intensidad) - .1, max(C$e$intensidad) + .3, length.out = 300)
pm1 <- ggplot() + stat_ecdf(data = C$e, aes(duracion), geom = "step", colour = "grey55", linewidth = .8) +
  geom_line(data = data.frame(x = xs, F = C$FD(xs)), aes(x, F), colour = azul, linewidth = 1.1) +
  annotate("segment", x = D0 - .5, xend = D0 - .5, y = 0, yend = pc0["u0"], linetype = "dashed") + annotate("segment", x = 2.9, xend = D0 - .5, y = pc0["u0"], yend = pc0["u0"], linetype = "dashed") +
  annotate("point", x = D0 - .5, y = pc0["u0"], shape = 8, size = 4.5, stroke = 1.3) +
  annotate("label", x = 9.5, y = .55, label = sprintf("Menos de %d días: u = %s\nAl menos %d días: 1 − u = %s", D0, coma(pc0["u0"], 2), D0, coma(1 - pc0["u0"], 2)), size = 4.6, fontface = "bold") +
  coord_cartesian(xlim = c(2.9, max(xs)), ylim = c(0, 1)) +
  scale_x_continuous(breaks = breaks_dias) +
  labs(title = "Duración", x = "Duración (días)", y = "Probabilidad acumulada", subtitle = sprintf("Gris: datos. Azul: marginal ajustada (%s)", toupper(C$x$mejor_ajuste_y_dist))) + tema
pm2 <- ggplot() + stat_ecdf(data = C$e, aes(intensidad), geom = "step", colour = "grey55", linewidth = .8) +
  geom_line(data = data.frame(x = xi, F = C$FI(xi)), aes(x, F), colour = col_calor, linewidth = 1.1) +
  annotate("segment", x = I0, xend = I0, y = 0, yend = pc0["v0"], linetype = "dashed") + annotate("segment", x = min(xi), xend = I0, y = pc0["v0"], yend = pc0["v0"], linetype = "dashed") +
  annotate("point", x = I0, y = pc0["v0"], shape = 8, size = 4.5, stroke = 1.3) +
  annotate("label", x = 3.55, y = .3, label = sprintf("F(intensidad 2) = %s", coma(pc0["v0"], 2)), size = 4.6, fontface = "bold") +
  coord_cartesian(xlim = range(xi), ylim = c(0, 1)) +
  labs(title = "Intensidad", x = "Intensidad", y = NULL, subtitle = sprintf("Gris: datos. Naranja: marginal ajustada (%s)", toupper(C$x$mejor_ajuste_x_dist))) + tema
ggsave(file.path(out, "caso_marginales_uv.png"), pm1 + pm2, width = 12, height = 5, dpi = 130)

pts <- data.frame(u = C$FD(C$e$duracion), v = C$FI(C$e$intensidad))
gr <- expand.grid(u = seq(.005, .995, length.out = 200), v = seq(.005, .995, length.out = 200))
gr$C <- pCopula(cbind(gr$v, gr$u), C$mv@copula); gr$ind <- gr$u * gr$v
br <- c(.05, .1, .2, .3, .5, .7)
pc <- ggplot() + geom_contour(data = gr, aes(u, v, z = ind), breaks = br, colour = "grey55", linetype = "dashed", linewidth = .6) +
  geom_contour(data = gr, aes(u, v, z = C), breaks = br, colour = azul, linewidth = 1) +
  geom_point(data = pts, aes(u, v), size = 2, alpha = .45, colour = col_calor) +
  annotate("point", x = pc0["u0"], y = pc0["v0"], shape = 8, size = 5, stroke = 1.4) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "F de la duración (u)", y = "F de la intensidad (v)", title = sprintf("Cópula ajustada: Normal, ρ = %s", coma(C$rho, 2)),
       subtitle = "Continuas: cópula ajustada. Punteadas: independencia") + tema + theme(plot.subtitle = element_text(size = 11))
ggsave(file.path(out, "caso_copula_ajustada.png"), pc, width = 8.1, height = 7.4, dpi = 130)

# Etiquetas de las curvas de nivel: punto medio de cada linea (igual criterio que el grafico del pipeline)
etiquetas_iso <- function(g, x, y, niveles) {
  l <- grDevices::contourLines(sort(unique(g[[x]])), sort(unique(g[[y]])), matrix(g$T[order(g[[y]], g[[x]])], nrow = length(unique(g[[x]]))), levels = niveles)
  bind_rows(lapply(l, function(z) { m <- ceiling(length(z$x) / 2); data.frame(d = z$x[m], i = z$y[m], nivel = z$level) }))
}
# ---- 6. isolineas (calor sin marcar, calor con el evento de diciembre de 2013, frio)
iso <- function(K, ttl, marca = NULL, etiqueta = NULL, dmax = NULL, imax = NULL) {
  dmax <- dmax %||% (max(K$e$duracion) + .5); imax <- imax %||% (max(K$e$intensidad) + .4)
  g <- expand.grid(d = seq(3, dmax, length.out = 220), i = seq(min(K$e$intensidad) - .1, imax, length.out = 220))
  g$T <- K$Tc(g$d, g$i)
  p <- ggplot(g, aes(d, i)) + geom_raster(aes(fill = T)) +
    scale_fill_viridis_c(trans = "log10", limits = c(2, 500), oob = scales::squish, breaks = c(2, 5, 10, 20, 50, 100, 200, 500), name = "Período de\nretorno (años)") +
    geom_contour(aes(z = T), breaks = c(2, 5, 10, 20, 50, 100, 200, 500), colour = "white", linewidth = .5) +
    geom_label(data = etiquetas_iso(g, "d", "i", c(2, 5, 10, 20, 50, 100, 200, 500)), aes(d, i, label = nivel), size = 4.2, linewidth = 0, label.padding = unit(.14, "lines"), fill = rgb(1, 1, 1, .8), colour = "grey20") +
    geom_point(data = K$e, aes(duracion, intensidad), shape = 21, fill = "white", colour = "grey20", size = 2.2, stroke = .4) +
    labs(x = "Duración (días)", y = "Intensidad", title = "Período de retorno combinado", subtitle = ttl) +
    theme_minimal(base_size = 16) + theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 17), plot.subtitle = element_text(size = 13, colour = "#56606B"), plot.background = element_rect(fill = "white", colour = NA)) +
    coord_cartesian(expand = FALSE) + scale_x_continuous(breaks = breaks_dias)
  if (!is.null(marca)) p <- p + geom_point(data = marca, aes(duracion, intensidad), shape = 8, size = 6, stroke = 1.5, colour = "#E0501F") +
    annotate("label", x = marca$duracion - 2.1, y = marca$intensidad + .5, label = etiqueta, size = 4.6, fontface = "bold", colour = "#B4491F")
  p
}
`%||%` <- function(a, b) if (is.null(a)) b else a
dic13 <- C$e %>% filter(fecha_inicio == as.Date("2013-12-22"))
Tdic <- C$Tc(dic13$duracion, dic13$intensidad)
u13 <- C$FD(dic13$duracion - .5); v13 <- C$FI(dic13$intensidad)
cat(sprintf("\nDic 2013 (TX): duracion %d, intensidad %.2f, T combinado %.0f anios; duracion sola: T=%.0f; intensidad sola: T=%.0f; indep: %.0f\n",
            dic13$duracion, dic13$intensidad, Tdic, 1 / (C$lam * (1 - u13)), 1 / (C$lam * (1 - v13)), 1 / (C$lam * (1 - u13) * (1 - v13))))
ggsave(file.path(out, "caso_isolineas_calor.png"), iso(C, "Períodos cálidos de TX, Buenos Aires (87585)"), width = 8.1, height = 5.94, dpi = 140)
ggsave(file.path(out, "caso_isolineas_evento.png"), iso(C, "Períodos cálidos de TX, Buenos Aires (87585)", dic13,
       sprintf("Dic. 2013: %d días, %s\nT combinado ≈ %.0f años", dic13$duracion, coma(dic13$intensidad, 2), Tdic)), width = 8.1, height = 5.94, dpi = 140)
ggsave(file.path(out, "caso_isolineas_frio.png"), iso(F, "Períodos fríos de TX, Buenos Aires (87585)"), width = 8.1, height = 5.94, dpi = 140)

# ---- 7. ENSO
cl <- ClasificarMesesEnso(LeerEventosEnso("data/input/enso_roni.txt"))
pl <- bind_rows(C$e %>% mutate(tipo = "Períodos cálidos de TX"), F$e %>% mutate(tipo = "Períodos fríos de TX")) %>%
  mutate(evento_enso = factor(AsignarEventoEnso(fecha_inicio, cl), levels = enso$niveles))
cat("\nENSO:\n"); print(table(pl$tipo, pl$evento_enso, useNA = "ifany"))
set.seed(3)
pe <- ggplot(pl, aes(duracion, intensidad, fill = evento_enso, shape = evento_enso)) +
  geom_point(size = 3.2, colour = "black", stroke = .4, alpha = .85) + scale_x_continuous(breaks = breaks_dias) +
  scale_fill_manual(values = enso$colores, drop = FALSE, name = "Fase ENSO") + scale_shape_manual(values = enso$formas, drop = FALSE, name = "Fase ENSO") +
  facet_wrap(~tipo, scales = "free") + labs(x = "Duración (días)", y = "Intensidad", subtitle = "Buenos Aires (87585), fase ENSO al inicio de cada período") +
  tema + theme(strip.text = element_text(face = "bold", size = 16), legend.position = "right")
ggsave(file.path(out, "caso_scatter_enso.png"), pe, width = 11.5, height = 5.4, dpi = 130)

# ---- 8. intensidad-magnitud (rapido)
mag <- lapply(list(list(C, "periodocalidotx"), list(F, "periodofriotx")), function(z) {
  K <- z[[1]]; m <- readRDS(sprintf("data/output/mejores_ajustes_multivariados_%s.rds", z[[2]])) %>% filter(station_id == "87585", variable_x == "intensidad", variable_y == "magnitud")
  tm <- cor(K$e$intensidad, K$e$magnitud, method = "kendall"); cat(sprintf("%s intensidad-magnitud: tau %.2f, familia %s\n", K$tp, tm, m$familia))
  ggplot(K$e, aes(intensidad, magnitud)) + geom_point(size = 2.2, alpha = .55, colour = K$col) +
    labs(title = sub("Períodos", "Períodos", ifelse(K$tp == "PeriodoCalidoTX", "Períodos cálidos de TX", "Períodos fríos de TX")), x = "Intensidad", y = "Magnitud",
         subtitle = sprintf("tau de Kendall = %s · cópula ganadora: %s", coma(tm, 2), m$familia)) + tema
})
ggsave(file.path(out, "caso_magnitud.png"), mag[[1]] + mag[[2]], width = 12, height = 4.8, dpi = 130)

# ---- 9. parsimonia (intensidad de los periodos calidos: ajusta con 20, evalua con el resto)
x <- C$e$intensidad
set.seed(42); idx <- sample(seq_along(x), 20); tr <- x[idx]; te <- x[-idx]
sh <- mean(tr)^2 / var(tr); rt <- mean(tr) / var(tr)
dsimple <- function(z) dgamma(z, shape = sh, rate = rt); bw <- 0.04; dover <- function(z) sapply(z, function(q) mean(dnorm(q, tr, bw)))
ll <- function(f, z) mean(log(pmax(f(z), 1e-12)))
cat(sprintf("\nparsimonia: held-out logL media simple %.2f, sobre %.2f (n_te=%d)\n", ll(dsimple, te), ll(dover, te), length(te)))
g <- data.frame(z = seq(0.2, max(x) + .5, length.out = 1500)); g$s <- dsimple(g$z); g$o <- dover(g$z)
base <- function() ggplot() + geom_histogram(data = data.frame(z = tr), aes(z, after_stat(density)), binwidth = 0.5, fill = "grey85", colour = "white") +
  geom_point(data = data.frame(z = tr), aes(z, y = 0), shape = 16, size = 2.6, colour = "grey30") +
  geom_point(data = data.frame(z = te), aes(z, y = 0), shape = 1, size = 2.8, colour = col_calor, stroke = 1) +
  coord_cartesian(ylim = c(0, 1.0), xlim = c(0.2, max(x) + .5)) + labs(x = "Intensidad", y = "Densidad") + tema
pA <- base() + geom_line(data = g, aes(z, s), colour = azul, linewidth = 1.4) +
  labs(title = "Parsimonioso: 2 parámetros", subtitle = sprintf("Puntaje con datos no usados (más alto, mejor): %s", coma(ll(dsimple, te), 2)))
pB <- base() + geom_line(data = g, aes(z, o), colour = col_calor, linewidth = 1.0) +
  labs(title = "Sobreajustado: 1 por dato", subtitle = sprintf("Puntaje con datos no usados (más alto, mejor): %s", coma(ll(dover, te), 2))) + theme(axis.title.y = element_blank())
ggsave(file.path(out, "caso_parsimonia.png"), pA + pB, width = 11, height = 4.6, dpi = 130)
sink()
