# Figuras de los filtros del pipeline: criterio de estacionariedad (Bonferroni) y embudo de pares que llegan a tener copula
# Se corre desde la raiz del repo, despues de resumen_ajustes.R y resumen_ajustes_filtrados.R
# Entradas: DIR_PRES/est_todos_tests.rds (p-valores de los 5 tests de estacionariedad en las series observadas, calculados fuera del pipeline
#   porque el pipeline corta al primer test que rechaza), data/output/ajustes_multivariados_*.rds, DIR_PRES/id_filtrado.csv
# Salidas: DIR_PRES/figuras (png) y DIR_PRES/filtros.txt
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
out <- file.path(DIR_PRES, "figuras")
suppressMessages({library(dplyr); library(ggplot2); library(patchwork)})
azul <- "#1F5F8B"; calor <- "#B4491F"; gris <- "#56606B"
tema <- theme_minimal(base_size = 15) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15),
  plot.subtitle = element_text(size = 12, colour = gris), plot.title.position = "plot", plot.background = element_rect(fill = "white", colour = NA))
sink(file.path(DIR_PRES, "filtros.txt"), split = TRUE)

# ---- Bonferroni: 5 tests de estacionariedad sobre las series observadas
e <- readRDS(file.path(DIR_PRES, "est_todos_tests.rds")) %>% filter(tipo_serie == "observada", as.character(station_id) != "87420")
nombres <- c(BoxPierceLjungBox = "Box-Ljung", PuntoCambioUni = "Cambio\nde nivel", PuntoCambioAutocopula = "Cambio de\nautocópula", MannKendall = "Mann-Kendall", EstacionaridadEmpirico = "Empírico")
rech <- bind_rows(lapply(c(.05, .01), function(a) e %>% group_by(test) %>% summarise(p_rech = 100 * mean(p <= a), .groups = "drop") %>% mutate(alfa = a)))
rech$test <- factor(nombres[rech$test], levels = nombres); rech$alfa <- factor(rech$alfa, levels = c(.05, .01), labels = c("α = 0,05 en cada test", "α = 0,01 en cada test (Bonferroni)"))
glob <- bind_rows(lapply(c(.05, .01), function(a) e %>% group_by(station_id, variable, tipo, n_serie) %>% summarise(est = all(p > a), .groups = "drop") %>% summarise(est = 100 * mean(est)) %>% mutate(alfa = a)))
cat("estacionarias (todos los tests pasan) en series observadas:\n"); print(as.data.frame(glob)); print(as.data.frame(rech %>% tidyr::pivot_wider(names_from = alfa, values_from = p_rech)))
cat(sprintf("probabilidad de rechazar alguna vez por azar con 5 tests independientes: %.1f %% (alfa .05), %.1f %% (alfa .01)\n", 100 * (1 - .95^5), 100 * (1 - .99^5)))
cat(sprintf("series observadas: %d\n", n_distinct(paste(e$station_id, e$variable, e$tipo))))
p1 <- ggplot(rech, aes(test, p_rech, fill = alfa)) + geom_col(position = position_dodge(.75), width = .7) +
  geom_text(aes(label = paste0(round(p_rech), " %")), position = position_dodge(.75), vjust = -.4, size = 4) +
  scale_fill_manual(values = c("grey70", azul), name = NULL) + coord_cartesian(ylim = c(0, 20)) +
  labs(x = NULL, y = "% de series observadas que rechaza", title = "Cada test rechaza menos con el umbral corregido", subtitle = "Series observadas de los 6 tipos de evento (H0: estacionaria)") + tema +
  theme(legend.position = "top", axis.text.x = element_text(size = 10.5))
d2 <- data.frame(umbral = factor(c("α = 0,05\nen cada test", "α = 0,01 en cada test\n(Bonferroni)"), levels = c("α = 0,05\nen cada test", "α = 0,01 en cada test\n(Bonferroni)")),
                 estac = glob$est, azar = c(100 * (1 - .95^5), 100 * (1 - .99^5)))
p2 <- ggplot(d2, aes(umbral, estac, fill = umbral)) + geom_col(width = .55) +
  geom_text(aes(label = paste0(round(estac), " %")), vjust = -.4, size = 5.4, fontface = "bold") +
  scale_fill_manual(values = c("grey70", azul), guide = "none") + coord_cartesian(ylim = c(0, 105)) +
  labs(x = NULL, y = "% de series que pasan los 5 tests", title = "Series consideradas estacionarias",
       subtitle = sprintf("Falso rechazo por azar con 5 tests:\n%d %% con 0,05 y %d %% con 0,01", round(d2$azar[1]), round(d2$azar[2]))) + tema
ggsave(file.path(out, "bonferroni.png"), p1 + p2 + plot_layout(widths = c(1.5, 1)), width = 13, height = 5.2, dpi = 130)

# ---- Embudo: de los pares estacion x tipo a los que tienen copula (par intensidad-duracion, serie observada)
r <- list(); for (s in c("olacalor", "olafrio", "periodocalidotn", "periodocalidotx", "periodofriotn", "periodofriotx")) {
  a <- readRDS(sprintf("data/output/ajustes_multivariados_%s.rds", s))
  r[[s]] <- a %>% filter(variable_x == "intensidad", variable_y == "duracion", tipo_serie == "observada") %>%
    filter(as.character(station_id) != "87420") %>% distinct(station_id, x_es_estacionaria, y_es_estacionaria, son_dependientes) %>% mutate(tipo = s) }
d <- bind_rows(r); id <- read.csv(file.path(DIR_PRES, "id_filtrado.csv")); stopifnot(nrow(id) == nrow(d))
dep <- sum(d$son_dependientes); est <- sum(d$x_es_estacionaria & d$y_es_estacionaria); amb <- sum(d$son_dependientes & d$x_es_estacionaria & d$y_es_estacionaria); cop <- sum(!is.na(id$familia))
cat(sprintf("\nembudo: %d pares; pasan dependencia %d; pasan estacionariedad %d; ambas %d; con copula ganadora %d\n", nrow(d), dep, est, amb, cop))
olas <- d %>% filter(grepl("^ola", tipo)); per <- d %>% filter(!grepl("^ola", tipo))
f <- function(x) c(pares = nrow(x), dep = sum(x$son_dependientes), est = sum(x$x_es_estacionaria & x$y_es_estacionaria), ambas = sum(x$son_dependientes & x$x_es_estacionaria & x$y_es_estacionaria))
cat("olas:", f(olas), " periodos:", f(per), "\n")
fu <- data.frame(paso = factor(c("Pares estación × tipo", "Pasa dependencia\n(3 tests)", "Pasa estacionariedad\n(5 tests, Bonferroni)", "Pasa las dos", "Con cópula ganadora"),
                               levels = rev(c("Pares estación × tipo", "Pasa dependencia\n(3 tests)", "Pasa estacionariedad\n(5 tests, Bonferroni)", "Pasa las dos", "Con cópula ganadora"))),
                 n = c(nrow(d), dep, est, amb, cop))
p3 <- ggplot(fu, aes(n, paso)) + geom_col(fill = c(azul, "#D98A66", "#D98A66", calor, calor)[5:1], width = .65) +
  geom_text(aes(label = n), hjust = -.2, size = 5.4, fontface = "bold") + scale_x_continuous(limits = c(0, 480), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Cuántos pares llegan a tener cópula", subtitle = "Par intensidad-duración; filtros sobre la serie observada") + tema + theme(axis.text.y = element_text(size = 13))
ggsave(file.path(out, "embudo.png"), p3, width = 8.4, height = 4.8, dpi = 130)
sink()

# ---- Dependencia: cuantos pares estacion x tipo pasan los tres tests (Kendall, Spearman, copula empirica; en secuencia)
# Salida: DIR_PRES/figuras/dependencia_tests.png (cifras de dependencia_tests.R, que lee data/partial/dependencia_*.rds)
dp <- readr::read_csv(file.path(DIR_PRES, "dependencia_tests.csv"), show_col_types = FALSE) %>% mutate(grupo = ifelse(grepl("^ola", tipo), "Olas", "Períodos"))
dg <- bind_rows(dp %>% group_by(grupo) %>% summarise(across(c(pares, pasan_los_tres, rechaza_kendall, rechaza_spearman, rechaza_empirica), sum), .groups = "drop"),
                dp %>% summarise(grupo = "Todos", across(c(pares, pasan_los_tres, rechaza_kendall, rechaza_spearman, rechaza_empirica), sum)))
dl <- dg %>% transmute(grupo = sprintf("%s\n(%d pares)", grupo, pares), `Pasan los tres tests` = pasan_los_tres / pares * 100,
                       `Rechaza Kendall` = rechaza_kendall / pares * 100, `Rechaza Spearman` = rechaza_spearman / pares * 100, `Rechaza cópula empírica` = rechaza_empirica / pares * 100) %>%
  tidyr::pivot_longer(-grupo, names_to = "resultado", values_to = "pct") %>% filter(resultado != "Rechaza Spearman")
dl$resultado <- factor(dl$resultado, levels = c("Rechaza cópula empírica", "Rechaza Kendall", "Pasan los tres tests"))
dl$grupo <- factor(dl$grupo, levels = rev(unique(dl$grupo[order(match(sub("\n.*", "", dl$grupo), c("Todos", "Olas", "Períodos")))])))
p4 <- ggplot(dl, aes(pct, grupo, fill = resultado)) + geom_col(width = .6) +
  geom_text(aes(label = ifelse(pct >= 4, paste0(round(pct), " %"), "")), position = position_stack(vjust = .5), colour = "white", fontface = "bold", size = 5) +
  scale_fill_manual(values = c("Pasan los tres tests" = azul, "Rechaza Kendall" = "#B4491F", "Rechaza cópula empírica" = "#E0A07F"), breaks = c("Pasan los tres tests", "Rechaza Kendall", "Rechaza cópula empírica"), name = NULL) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 100), labels = function(x) paste0(x, " %")) +
  labs(x = NULL, y = NULL, title = "Dependencia: la mayoría de las olas no pasa los tests", subtitle = "Par intensidad–duración, serie observada; Spearman no rechaza ningún par que pasó Kendall") + tema +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 13))
ggsave(file.path(out, "dependencia_tests.png"), p4, width = 9.2, height = 3.9, dpi = 130)
