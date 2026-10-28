# Figuras: parsimonioso vs sobreajustado, AIC por familia con n chico, dependencia debil y ajustes sin ganadora (requiere aic_por_familia.R)
# Se corre desde la raiz del repo: Rscript scripts/presentacion/figuras_parsimonia_limites.R
# Entradas: data/output/*.rds (corrida Olas_20), data/input/eventos_identificados_frio_calor.csv
# Salidas: DIR_PRES (por defecto data/output/presentacion); las figuras van en DIR_PRES/figuras
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
dir.create(file.path(DIR_PRES, "figuras"), recursive = TRUE, showWarnings = FALSE)
REGIONES_CSV <- Sys.getenv("REGIONES_CSV", "data/input/regiones.csv")
suppressMessages({library(dplyr); library(ggplot2); library(patchwork); library(copula)})
source("lib/funciones_ajuste_distribuciones.R", encoding="UTF-8")
out <- file.path(DIR_PRES, "figuras")
azul <- "#1F5F8B"; calor <- "#B4491F"; frio <- "#1F6F9F"; gris <- "#56606B"
tema <- theme_minimal(base_size = 15) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15),
  plot.subtitle = element_text(size = 12, colour = gris), plot.title.position = "plot", plot.background = element_rect(fill = "white", colour = NA))
ev <- readr::read_csv("data/input/eventos_identificados_frio_calor.csv", show_col_types = FALSE) %>% filter(duracion >= 3)

# ---- A: parsimonioso vs sobreajustado (intensidad, olas de calor 87585)
x <- ev %>% filter(station_id == "87585", tipo_evento == "OlaCalorTXyTN") %>% pull(intensidad)
set.seed(42); idx <- sample(seq_along(x), 20); tr <- x[idx]; te <- x[-idx]
sh <- mean(tr)^2 / var(tr); rt <- mean(tr) / var(tr)            # gamma por momentos: 2 parametros
dsimple <- function(z) dgamma(z, shape = sh, rate = rt)
bw <- 0.04; dover <- function(z) sapply(z, function(q) mean(dnorm(q, tr, bw)))
ll <- function(f, z) mean(log(pmax(f(z), 1e-12)))
cat(sprintf("held-out logL media: simple %.2f, sobre %.2f (n_te=%d)\n", ll(dsimple, te), ll(dover, te), length(te)))
g <- data.frame(z = seq(0.3, 6, length.out = 1500)); g$s <- dsimple(g$z); g$o <- dover(g$z)
base <- function() ggplot() + geom_histogram(data = data.frame(z = tr), aes(z, after_stat(density)), binwidth = 0.5, fill = "grey85", colour = "white") +
  geom_point(data = data.frame(z = tr), aes(z, y = 0), shape = 16, size = 2.6, colour = "grey30") +
  geom_point(data = data.frame(z = te), aes(z, y = 0), shape = 1, size = 2.8, colour = calor, stroke = 1) +
  coord_cartesian(ylim = c(0, 1.0), xlim = c(0.3, 6)) + labs(x = "Intensidad", y = "Densidad") + tema
pA <- base() + geom_line(data = g, aes(z, s), colour = azul, linewidth = 1.4) +
  labs(title = "Parsimonioso: 2 parámetros", subtitle = sprintf("Puntaje con datos no usados (más alto, mejor): %s", format(round(ll(dsimple, te), 2), decimal.mark = ",")))
pB <- base() + geom_line(data = g, aes(z, o), colour = calor, linewidth = 1.0) +
  labs(title = "Sobreajustado: 1 por dato", subtitle = sprintf("Puntaje con datos no usados (más alto, mejor): %s", format(round(ll(dover, te), 2), decimal.mark = ","))) + theme(axis.title.y = element_blank())
ggsave(file.path(out, "parsimonia.png"), pA + pB, width = 11, height = 4.6, dpi = 130)

# ---- B1: AIC por familia, estación 87679, olas de frío
r <- readRDS(file.path(DIR_PRES, "aic_olafrio.rds")) %>% filter(station_id == "87679") %>% mutate(familia = reorder(familia, aic))
r$ganadora <- r$aic == min(r$aic)
n87679 <- nrow(ev %>% filter(station_id == "87679", tipo_evento == "OlaFrioTXyTN"))
pB1 <- ggplot(r, aes(familia, aic, fill = ganadora)) + geom_col(width = .65) +
  geom_text(aes(label = format(round(aic, 2), decimal.mark = ",")), vjust = -.4, size = 4.2) +
  scale_fill_manual(values = c(`TRUE` = frio, `FALSE` = "grey75"), guide = "none") +
  labs(title = "AIC por familia (menor es mejor)", subtitle = sprintf("Olas de frío, estación 87679, %d olas", n87679), x = NULL, y = "AIC") +
  coord_cartesian(ylim = c(0, 3.4)) + tema
ggsave(file.path(out, "pocos_aic.png"), pB1, width = 5.4, height = 4.0, dpi = 130)

# ---- B2: período de retorno con cópula vs independencia (frío BA)
d2 <- data.frame(m = factor(c("Con la cópula", "Suponiendo independencia"), levels = c("Con la cópula", "Suponiendo independencia")), T = c(19.3, 24.8))
pB2 <- ggplot(d2, aes(m, T, fill = m)) + geom_col(width = .55) +
  geom_text(aes(label = paste(format(T, decimal.mark = ","), "años")), vjust = -.4, size = 4.6, fontface = "bold") +
  scale_fill_manual(values = c(frio, "grey70"), guide = "none") +
  labs(title = "Período de retorno combinado", subtitle = "Olas de frío, Buenos Aires: ≥ 4 días e intensidad ≥ 3 (tau 0,05)", x = NULL, y = "Años") +
  coord_cartesian(ylim = c(0, 30)) + tema
ggsave(file.path(out, "dependencia_debil.png"), pB2, width = 5.4, height = 4.0, dpi = 130)

# ---- B3: tau de Kendall de los 6 ajustes sin ganadora
d3 <- data.frame(est = c("Iguazú, frío", "Roque Sáenz Peña, frío", "B. de Irigoyen, frío", "Monte Caseros, frío", "Concordia, frío", "Río Gallegos, cálido TN"),
                 tau = c(-0.195, -0.264, -0.132, -0.088, -0.118, -0.078))
d3$est <- factor(d3$est, levels = rev(d3$est[order(d3$tau)]))
pB3 <- ggplot(d3, aes(tau, est)) + geom_col(fill = "#8A5A44", width = .6) +
  geom_text(aes(label = format(round(tau, 2), decimal.mark = ",")), hjust = 1.1, size = 4.2) +
  scale_x_continuous(limits = c(-0.34, 0.05)) +
  labs(title = "Ajustes sin cópula ganadora", subtitle = "Los 6 tienen tau de Kendall negativo", x = "Tau de Kendall", y = NULL) + tema +
  theme(axis.text.y = element_text(size = 11))
ggsave(file.path(out, "sin_ganadora.png"), pB3, width = 5.4, height = 4.0, dpi = 130)
