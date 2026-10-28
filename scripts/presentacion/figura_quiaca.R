# Limitacion del metodo: tendencia en la frecuencia de eventos (ejemplo La Quiaca, olas de calor) y cuantos pares con copula la tienen
# Se corre desde la raiz del repo, despues de resumen_ajustes_filtrados.R
# Entradas: data/input/eventos_identificados_frio_calor.csv, DIR_PRES/id_filtrado.csv
# Salidas: DIR_PRES/figuras/quiaca_frecuencia.png y DIR_PRES/quiaca.txt
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
out <- file.path(DIR_PRES, "figuras")
suppressMessages({library(dplyr); library(ggplot2); library(Kendall)})
azul <- "#1F5F8B"; calor <- "#B4491F"; gris <- "#56606B"
tema <- theme_minimal(base_size = 16) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 16),
  plot.subtitle = element_text(size = 12.5, colour = gris), plot.title.position = "plot", plot.background = element_rect(fill = "white", colour = NA))
sink(file.path(DIR_PRES, "quiaca.txt"), split = TRUE)
ev <- readr::read_csv("data/input/eventos_identificados_frio_calor.csv", show_col_types = FALSE) %>% filter(duracion >= 3) %>%
  mutate(station_id = as.character(station_id), anio = as.integer(format(as.Date(fecha_inicio), "%Y")))
anual <- function(e) data.frame(anio = 1961:2025) %>% left_join(e %>% count(anio), by = "anio") %>% mutate(n = ifelse(is.na(n), 0L, n))

# ---- La Quiaca (87007), olas de calor
q <- anual(ev %>% filter(station_id == "87007", tipo_evento == "OlaCalorTXyTN"))
mk <- MannKendall(q$n); pm <- glm(n ~ anio, poisson, q)
tasa_total <- sum(q$n) / nrow(q); tasa_10 <- sum(q$n[q$anio >= 2016]) / 10
tasa_91 <- sum(q$n[q$anio >= 1991 & q$anio <= 2020]) / 30   # normal climatologica 1991-2020
cat(sprintf("olas: %d; por anio: %.2f; ultimos 10 anios (2016-2025): %.2f; antes de 1998: %d olas en 37 anios\n", sum(q$n), tasa_total, tasa_10, sum(q$n[q$anio < 1998])))
cat(sprintf("Mann-Kendall sobre conteos anuales: tau = %.2f, p = %.1e; Poisson: pendiente %.3f por anio (%.1f %% anual)\n", mk$tau, mk$sl, coef(pm)[2], 100 * (exp(coef(pm)[2]) - 1)))
P <- 0.2
cat(sprintf("Tasa 1991-2020: %.2f por anio (%d olas en 30 anios); con P = %.1f: T = %.1f anios; desde 2010: %d olas, 2010-2019: %d\n", tasa_91, sum(q$n[q$anio >= 1991 & q$anio <= 2020]), P, 1 / (tasa_91 * P), sum(q$n[q$anio >= 2010]), sum(q$n[q$anio >= 2010 & q$anio <= 2019])))
cat(sprintf("Ejemplo con P = %.1f: T = %.1f anios con la tasa de todo el registro; %.1f con la de los ultimos 10 anios\n", P, 1 / (tasa_total * P), 1 / (tasa_10 * P)))
q$ajuste <- predict(pm, type = "response")
p <- ggplot(q, aes(anio, n)) + geom_col(fill = calor, width = .8) +
  geom_hline(yintercept = tasa_total, colour = azul, linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = ajuste), colour = "black", linewidth = 1) +
  annotate("segment", x = 1991 - .4, xend = 2020 + .4, y = tasa_91, yend = tasa_91, colour = "#2E8B57", linewidth = 1.6) +
  annotate("text", x = 1961, y = 3.3, hjust = 0, colour = "#2E8B57", size = 5, fontface = "bold", label = paste0("Tasa 1991-2020 (línea verde): ", sub(".", ",", sprintf("%.2f", tasa_91), fixed = TRUE), " olas por año")) +
  annotate("text", x = 1961, y = tasa_total + .75, hjust = 0, colour = azul, size = 5, label = paste0("Tasa supuesta: ", sub(".", ",", sprintf("%.2f", tasa_total), fixed = TRUE), " olas por año")) +
  annotate("text", x = 1961, y = 2.2, hjust = 0, colour = "black", size = 5, label = "Tasa ajustada (regresión de Poisson)") +
  scale_x_continuous(breaks = seq(1965, 2025, 10)) + scale_y_continuous(breaks = 0:6, limits = c(0, 6)) +
  labs(x = NULL, y = "Olas de calor por año", title = "La Quiaca: la frecuencia de olas de calor crece",
       subtitle = sprintf("25 olas en 65 años; Mann-Kendall sobre los conteos anuales: p < 0,001 (tau = %s)", sub(".", ",", sprintf("%.2f", mk$tau), fixed = TRUE))) + tema
ggsave(file.path(out, "quiaca_frecuencia.png"), p, width = 9.2, height = 5.4, dpi = 130)

# ---- cuantos de los 196 pares con copula tienen tendencia en la frecuencia (Mann-Kendall sobre conteos anuales)
id <- read.csv(file.path(DIR_PRES, "id_filtrado.csv")) %>% filter(!is.na(familia)) %>% mutate(station_id = as.character(station_id))
r <- bind_rows(lapply(seq_len(nrow(id)), function(i) { z <- anual(ev %>% filter(station_id == id$station_id[i], tipo_evento == id$tipo[i]))
  m <- MannKendall(z$n); data.frame(station_id = id$station_id[i], tipo = id$tipo[i], clase = id$clase[i], n = sum(z$n), tau = m$tau, p = m$sl) }))
cat(sprintf("\npares con copula: %d\n", nrow(r)))
for (a in c(.05, .01)) cat(sprintf("tendencia en frecuencia (MK sobre conteos anuales) p < %.2f: %d de %d (%.0f %%)\n", a, sum(r$p < a), nrow(r), 100 * mean(r$p < a)))
cat("p < 0,01 por clase y signo:\n"); print(as.data.frame(r %>% filter(p < .01) %>% count(clase, signo = ifelse(tau > 0, "aumenta", "disminuye"))))
cat("p < 0,01 por tipo:\n"); print(as.data.frame(r %>% group_by(tipo) %>% summarise(pares = n(), con_tendencia = sum(p < .01))))
write.csv(r, file.path(DIR_PRES, "tendencia_frecuencia.csv"), row.names = FALSE)
sink()
