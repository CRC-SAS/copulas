# Figura: cuatro tipos de copulas por cola, con muestras simuladas
# Se corre desde la raiz del repo: Rscript scripts/presentacion/figura_tipos_copulas.R
# Entradas: data/output/*.rds (corrida Olas_20), data/input/eventos_identificados_frio_calor.csv
# Salidas: DIR_PRES (por defecto data/output/presentacion); las figuras van en DIR_PRES/figuras
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
dir.create(file.path(DIR_PRES, "figuras"), recursive = TRUE, showWarnings = FALSE)
REGIONES_CSV <- Sys.getenv("REGIONES_CSV", "data/input/regiones.csv")
suppressMessages({library(dplyr); library(ggplot2); library(copula)})
out <- file.path(DIR_PRES, "figuras")
set.seed(21); n <- 450; tau <- 0.5
sim <- function(cop, g) { x <- rCopula(n, cop); data.frame(u = x[,1], v = x[,2], g = g) }
d <- bind_rows(sim(normalCopula(iTau(normalCopula(), tau)), "Sin cola"),
               sim(gumbelCopula(iTau(gumbelCopula(), tau)), "Cola superior"),
               sim(claytonCopula(iTau(claytonCopula(), tau)), "Cola inferior"),
               sim(tCopula(iTau(tCopula(df = 3, df.fixed = TRUE), tau), df = 3, df.fixed = TRUE), "Ambas colas"))
niv <- c("Sin cola", "Cola superior", "Cola inferior", "Ambas colas"); d$g <- factor(d$g, levels = niv)
cols <- c("Sin cola" = "#7C8A97", "Cola superior" = "#B4491F", "Cola inferior" = "#1F5F8B", "Ambas colas" = "#6B4E9B")
sub <- c("Sin cola" = "Normal, Frank", "Cola superior" = "Gumbel, Joe", "Cola inferior" = "Clayton", "Ambas colas" = "t de Student")
lab <- setNames(paste0(niv, "\n", sub[niv]), niv)
rect <- bind_rows(data.frame(g = "Cola superior", x0 = .8, x1 = 1, y0 = .8, y1 = 1), data.frame(g = "Cola inferior", x0 = 0, x1 = .2, y0 = 0, y1 = .2),
                  data.frame(g = "Ambas colas", x0 = .8, x1 = 1, y0 = .8, y1 = 1), data.frame(g = "Ambas colas", x0 = 0, x1 = .2, y0 = 0, y1 = .2))
rect$g <- factor(rect$g, levels = niv)
p <- ggplot(d, aes(u, v, colour = g)) +
  geom_rect(data = rect, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = g), inherit.aes = FALSE, alpha = .18) +
  geom_point(size = 1.1, alpha = .55) + facet_wrap(~g, nrow = 1, labeller = as_labeller(lab)) +
  scale_colour_manual(values = cols, guide = "none") + scale_fill_manual(values = cols, guide = "none") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) + labs(x = "Valores bajos  →  altos de la primera variable", y = "Segunda variable") +
  theme_minimal(base_size = 15) + theme(panel.grid.minor = element_blank(), axis.text = element_blank(), panel.grid.major = element_line(colour = "grey92"),
    strip.text = element_text(face = "bold", size = 16, lineheight = 1.1), plot.background = element_rect(fill = "white", colour = NA), axis.title = element_text(size = 12, colour = "#56606B"))
ggsave(file.path(out, "tipos_copulas.png"), p, width = 13.6, height = 4.1, dpi = 130)
