# Figuras de "Donde hay que ser cuidadosos": AIC por familia con n chico y periodo de retorno con copula contra independencia (87585)
# Se corre desde la raiz del repo, despues de resumen_ajustes.R, resumen_ajustes_filtrados.R y figuras_caso_87585.R
# Entradas: data/output/ajustes_multivariados_*.rds, DIR_PRES/id_filtrado.csv, data/output/copulas_periodocalidotx.rds
# Salidas: DIR_PRES/figuras (png) y DIR_PRES/limites.txt
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
out <- file.path(DIR_PRES, "figuras")
suppressMessages({library(dplyr); library(ggplot2)})
azul <- "#1F5F8B"; calor <- "#B4491F"; frio <- "#1F6F9F"; gris <- "#56606B"
tema <- theme_minimal(base_size = 15) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15),
  plot.subtitle = element_text(size = 12, colour = gris), plot.title.position = "plot", plot.background = element_rect(fill = "white", colour = NA))
coma <- function(x, d = 1) format(round(x, d), decimal.mark = ",", nsmall = d)
sink(file.path(DIR_PRES, "limites.txt"), split = TRUE)
id <- read.csv(file.path(DIR_PRES, "id_filtrado.csv")) %>% mutate(station_id = as.character(station_id))
chicos <- id %>% filter(!is.na(familia), n < 30)
cat("ajustes con copula y n < 30:\n"); print(chicos %>% select(station_id, tipo, familia, n, tau_id))
# AIC por familia (mediana sobre las series perturbadas que se ajustaron) para esos ajustes
res <- list()
for (i in seq_len(nrow(chicos))) {
  s <- c(OlaCalorTXyTN = "olacalor", OlaFrioTXyTN = "olafrio")[[chicos$tipo[i]]]
  a <- readRDS(sprintf("data/output/ajustes_multivariados_%s.rds", s)) %>% filter(station_id == chicos$station_id[i], variable_x == "intensidad", variable_y == "duracion")
  res[[i]] <- bind_rows(lapply(seq_len(nrow(a)), function(k) { b <- a$bondad.ajuste[[k]]; if (!is.list(b)) return(NULL); e <- b$estadisticos
    e <- e[e$test == "AIC", , drop = FALSE]; if (!nrow(e)) return(NULL); data.frame(familia = a$familia[k], aic = e$valor) })) %>% group_by(familia) %>% summarise(aic = median(aic), .groups = "drop") %>%
    arrange(aic) %>% mutate(station_id = chicos$station_id[i], tipo = chicos$tipo[i], n = chicos$n[i], dif = aic[2] - aic[1])
}
r <- bind_rows(res); print(as.data.frame(r %>% group_by(station_id, tipo) %>% summarise(n = first(n), g = first(familia), dif = first(dif), .groups = "drop")))
SEL <- Sys.getenv("EST_AIC", "")  # estacion elegida, se fija despues de ver la tabla
if (nzchar(SEL)) {
  z <- r %>% filter(station_id == SEL) %>% mutate(ganadora = aic == min(aic), familia = reorder(familia, aic))
  nm <- z$n[1]; tp <- z$tipo[1]
  p <- ggplot(z, aes(familia, aic, fill = ganadora)) + geom_col(width = .65) +
    geom_text(aes(label = coma(aic, 2)), vjust = ifelse(z$aic >= 0, -.4, 1.3), size = 4.2) +
    scale_fill_manual(values = c(`TRUE` = ifelse(grepl("Calor", tp), calor, frio), `FALSE` = "grey75"), guide = "none") +
    labs(title = "AIC por familia (menor es mejor)", subtitle = sprintf("%s, estación %s, %d eventos", ifelse(grepl("Calor", tp), "Olas de calor", "Olas de frío"), SEL, nm), x = NULL, y = "AIC") + tema
  ggsave(file.path(out, "pocos_aic.png"), p, width = 5.4, height = 4.0, dpi = 130)
}
# Periodo de retorno con copula e independencia: 87585, periodos calidos de TX, al menos 5 dias (evaluada en 4,5) e intensidad mayor a 2 (valores de caso_87585.txt)
d2 <- data.frame(m = factor(c("Con la cópula", "Suponiendo independencia"), levels = c("Con la cópula", "Suponiendo independencia")), T = c(3.0, 3.7))
pB2 <- ggplot(d2, aes(m, T, fill = m)) + geom_col(width = .55) +
  geom_text(aes(label = paste(coma(T, 1), "años")), vjust = -.4, size = 4.6, fontface = "bold") +
  scale_fill_manual(values = c(calor, "grey70"), guide = "none") +
  labs(title = "Período de retorno combinado", subtitle = "Períodos cálidos de TX, Buenos Aires: al menos 5 días e intensidad mayor a 2 (tau 0,16)", x = NULL, y = "Años") +
  coord_cartesian(ylim = c(0, 5)) + tema + theme(plot.subtitle = element_text(size = 10.5))
ggsave(file.path(out, "dependencia_debil.png"), pB2, width = 5.4, height = 4.0, dpi = 130)
sink()
