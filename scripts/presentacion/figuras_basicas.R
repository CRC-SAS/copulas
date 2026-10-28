# Figuras: scatter 87585, de los datos a la copula, familias y colas simuladas, copula ajustada y scatter por fase ENSO
# Se corre desde la raiz del repo: Rscript scripts/presentacion/figuras_basicas.R
# Entradas: data/output/*.rds (corrida Olas_20), data/input/eventos_identificados_frio_calor.csv
# Salidas: DIR_PRES (por defecto data/output/presentacion); las figuras van en DIR_PRES/figuras
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
dir.create(file.path(DIR_PRES, "figuras"), recursive = TRUE, showWarnings = FALSE)
REGIONES_CSV <- Sys.getenv("REGIONES_CSV", "data/input/regiones.csv")
suppressMessages({library(dplyr); library(ggplot2); library(patchwork); library(copula); library(lubridate)})
for (f in c("lib/funciones_enso.R")) source(f, encoding = "UTF-8")
txt <- readLines("lib/funciones_periodo_retorno.R", encoding = "UTF-8")
i <- grep("^ConstantesEnso <- function", txt); j <- i + which(txt[(i+1):length(txt)] == "}")[1]
eval(parse(text = txt[i:j])); enso <- ConstantesEnso()
out <- file.path(DIR_PRES, "figuras")
tema <- theme_minimal(base_size = 16) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 17),
  plot.background = element_rect(fill = "white", colour = NA))
col_calor <- "#B4491F"; col_frio <- "#1F6F9F"; azul <- "#1F5F8B"

ev <- readr::read_csv("data/input/eventos_identificados_frio_calor.csv", show_col_types = FALSE) %>%
  filter(station_id == "87585", duracion >= 3) %>% mutate(fecha_inicio = as.Date(fecha_inicio))
calor <- ev %>% filter(tipo_evento == "OlaCalorTXyTN"); frio <- ev %>% filter(tipo_evento == "OlaFrioTXyTN")
cat("n calor", nrow(calor), "n frio", nrow(frio), "\n")
tau_c <- cor(calor$intensidad, calor$duracion, method = "kendall")
cat("tau", tau_c, "\n")

# 5: scatter
set.seed(1)
p5 <- ggplot(calor, aes(duracion, intensidad)) +
  geom_jitter(width = 0.12, height = 0, size = 3.2, alpha = .75, colour = col_calor) +
  labs(title = "Olas de calor, Buenos Aires (87585)", x = "Duración (días)", y = "Intensidad",
       subtitle = sprintf("n = %d olas · tau de Kendall = %s", nrow(calor), format(round(tau_c, 2), decimal.mark = ","))) + tema
ggsave(file.path(out, "scatter_calor.png"), p5, width = 8.1, height = 5.94, dpi = 140)

# 8: de los datos a la copula
u <- pobs(cbind(calor$duracion, calor$intensidad)); colnames(u) <- c("u", "v")
fit <- fitCopula(gumbelCopula(2), u, method = "mpl")
th <- 1.165608; cat("theta gumbel (pipeline)", th, "\n")
a <- ggplot(calor, aes(duracion, intensidad)) + geom_jitter(width = .12, size = 2.6, alpha = .75, colour = col_calor) +
  labs(title = "1. Los datos", x = "Duración", y = "Intensidad") + tema
b <- ggplot(as.data.frame(u), aes(u, v)) + geom_point(size = 2.6, alpha = .75, colour = azul) +
  coord_fixed() + labs(title = "2. Sus rangos", x = "Rango de duración", y = "Rango de intensidad") + tema
set.seed(11); sm <- as.data.frame(rCopula(500, gumbelCopula(th))); names(sm) <- c("u", "v")
c3 <- ggplot(sm, aes(u, v)) + geom_point(size = 1.4, alpha = .45, colour = azul) +
  coord_fixed() + labs(title = "3. La cópula ajustada", x = "Rango de duración", y = "Rango de intensidad") + tema
ggsave(file.path(out, "de_datos_a_copula.png"), a + b + c3 + plot_layout(ncol = 3), width = 13, height = 4.6, dpi = 130)

# 9: familias y colas
set.seed(7); n <- 600
sim <- function(cop, nom) { x <- rCopula(n, cop); data.frame(u = x[,1], v = x[,2], fam = nom) }
d <- bind_rows(sim(normalCopula(iTau(normalCopula(), .5)), "Normal: sin colas"),
               sim(gumbelCopula(iTau(gumbelCopula(), .5)), "Gumbel: cola superior"),
               sim(claytonCopula(iTau(claytonCopula(), .5)), "Clayton: cola inferior"),
               sim(tCopula(iTau(tCopula(df = 3, df.fixed = TRUE), .5), df = 3, df.fixed = TRUE), "t: ambas colas"))
d$fam <- factor(d$fam, levels = unique(d$fam))
p9 <- ggplot(d, aes(u, v)) + geom_point(size = .9, alpha = .5, colour = azul) + facet_wrap(~fam, nrow = 2) + coord_fixed() +
  labs(x = "Rango de la primera variable", y = "Rango de la segunda variable", subtitle = "Mismo tau de Kendall (0,5) en los cuatro paneles") +
  tema + theme(strip.text = element_text(face = "bold", size = 16))
ggsave(file.path(out, "familias_colas.png"), p9, width = 8.1, height = 7.4, dpi = 130)

# 20: scatter coloreado por ENSO
cl <- ClasificarMesesEnso(LeerEventosEnso("data/input/enso_roni.txt"))
pl <- bind_rows(calor %>% mutate(tipo = "Olas de calor"), frio %>% mutate(tipo = "Olas de frío")) %>%
  mutate(evento_enso = AsignarEventoEnso(fecha_inicio, cl)) %>%
  mutate(evento_enso = factor(evento_enso, levels = enso$niveles))
print(table(pl$tipo, pl$evento_enso, useNA = "ifany"))
set.seed(3)
p20 <- ggplot(pl, aes(duracion, intensidad, fill = evento_enso, shape = evento_enso)) +
  geom_jitter(width = .12, height = 0, size = 3.4, colour = "black", stroke = .4, alpha = .95) +
  scale_fill_manual(values = enso$colores, drop = FALSE, name = "Fase ENSO") +
  scale_shape_manual(values = enso$formas, drop = FALSE, name = "Fase ENSO") +
  facet_wrap(~tipo, scales = "free_y") + labs(x = "Duración (días)", y = "Intensidad", subtitle = "Buenos Aires (87585), fase ENSO al inicio de cada ola") +
  tema + theme(strip.text = element_text(face = "bold", size = 16), legend.position = "right")
ggsave(file.path(out, "scatter_enso.png"), p20, width = 11.5, height = 5.4, dpi = 130)

# Cópula ajustada (parámetros del pipeline) con ejemplo
source("lib/funciones_ajuste_distribuciones.R", encoding = "UTF-8")
x_ <- readRDS("data/output/copulas_olacalor.rds") %>% filter(station_id == "87585", variable_x == "intensidad", variable_y == "duracion", tipo_serie == "observada")
mv <- x_$ajuste_multivariado[[1]]
FI <- function(q) do.call("pglo", c(list(q = q), mv@paramMargins[[1]]))
FD <- function(q) do.call("ppe3", c(list(q = q), mv@paramMargins[[2]]))
pts <- data.frame(u = FI(calor$intensidad), v = FD(calor$duracion))
gr <- expand.grid(u = seq(.005, .995, length.out = 200), v = seq(.005, .995, length.out = 200))
gr$C <- pCopula(cbind(gr$u, gr$v), mv@copula); gr$ind <- gr$u * gr$v
u0 <- FI(3); v0 <- FD(5); C0 <- pCopula(c(u0, v0), mv@copula)
cat(sprintf("u0=%.3f v0=%.3f C0=%.3f ind=%.3f\n", u0, v0, C0, u0*v0))
br <- c(.05, .1, .2, .3, .5, .7)
pc <- ggplot() +
  geom_contour(data = gr, aes(u, v, z = ind), breaks = br, colour = "grey55", linetype = "dashed", linewidth = .6) +
  geom_contour(data = gr, aes(u, v, z = C), breaks = br, colour = azul, linewidth = 1) +
  geom_point(data = pts, aes(u, v), size = 2.2, alpha = .55, colour = col_calor) +
  annotate("point", x = u0, y = v0, shape = 8, size = 5, stroke = 1.4) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "F de la intensidad", y = "F de la duración",
       title = "Cópula ajustada: Gumbel, θ = 1,17",
       subtitle = "Continuas: cópula ajustada. Punteadas: independencia") + tema +
  theme(plot.subtitle = element_text(size = 11))
ggsave(file.path(out, "copula_ajustada.png"), pc, width = 8.1, height = 7.4, dpi = 130)
