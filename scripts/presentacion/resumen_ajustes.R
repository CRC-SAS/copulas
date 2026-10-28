# Resume los mejores ajustes (par intensidad-duracion, 6 tipos de evento) por estacion, region, cola de la copula, n de eventos y tau de Kendall. Genera DIR_PRES/resumen_ajustes.csv
# Se corre desde la raiz del repo: Rscript scripts/presentacion/resumen_ajustes.R
# Entradas: data/output/*.rds (corrida Olas_20), data/input/eventos_identificados_frio_calor.csv
# Salidas: DIR_PRES (por defecto data/output/presentacion); las figuras van en DIR_PRES/figuras
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
dir.create(file.path(DIR_PRES, "figuras"), recursive = TRUE, showWarnings = FALSE)
REGIONES_CSV <- Sys.getenv("REGIONES_CSV", "data/input/regiones.csv")
suppressMessages(library(dplyr)); suppressMessages(library(readr))
reg <- read.csv2(REGIONES_CSV, fileEncoding="latin1")
reg <- data.frame(station_id=as.character(reg[[1]]), region=toupper(reg[[ncol(reg)]]))
reg <- bind_rows(reg, data.frame(station_id=c("87345","87420","87569","87934"), region=c("CENTRO","CUYO","BUENOS AIRES","PATAGONIA")))
ev <- read_csv("data/input/eventos_identificados_frio_calor.csv", show_col_types=FALSE) %>% mutate(station_id=as.character(station_id))
tipos <- c(olacalor="OlaCalorTXyTN", olafrio="OlaFrioTXyTN", periodocalidotn="PeriodoCalidoTN",
           periodocalidotx="PeriodoCalidoTX", periodofriotn="PeriodoFrioTN", periodofriotx="PeriodoFrioTX")
res <- bind_rows(lapply(names(tipos), function(s) {
  m <- readRDS(sprintf("data/output/mejores_ajustes_multivariados_%s.rds", s)) %>% mutate(station_id=as.character(station_id), tipo=tipos[[s]], slug=s)
  m
}))
tau <- ev %>% group_by(station_id, tipo=tipo_evento) %>%
  summarise(n=n(), tau_id=suppressWarnings(cor(intensidad, duracion, method="kendall")),
            tau_im=suppressWarnings(cor(intensidad, magnitud, method="kendall")), .groups="drop")
# 87420 (Mendoza Observatorio) se excluye de todos los resultados (decision de Natalia, 2026-10-23)
EXCLUIDA <- "87420"
res <- res %>% filter(station_id != EXCLUIDA) %>% left_join(tau, by=c("station_id","tipo")) %>% left_join(reg, by="station_id") %>%
  mutate(cola=case_when(familia %in% c("normal","frank")~"Sin cola", familia %in% c("gumbel","joe")~"Superior",
                        familia=="clayton"~"Inferior", familia=="t"~"Ambas", TRUE~familia),
         clase=ifelse(grepl("^Ola",tipo),"Ola","Periodo"), estacion_tipo=ifelse(grepl("Calor|Calido",tipo),"Calor","Frio"))
write_csv(res, file.path(DIR_PRES, "resumen_ajustes.csv"))
id <- res %>% filter(variable_x=="intensidad", variable_y=="duracion")
cat("sin region:", sum(is.na(id$region)), "\n")
pct <- function(d, ...) d %>% group_by(...) %>% summarise(est=n(), sin=round(100*mean(cola=="Sin cola")), sup=round(100*mean(cola=="Superior")), inf=round(100*mean(cola=="Inferior")), amb=round(100*mean(cola=="Ambas")), n_med=median(n), tau_med=round(median(tau_id,na.rm=TRUE),2), .groups="drop")
print(as.data.frame(pct(id, region)))
print(as.data.frame(pct(id, estacion_tipo)))
print(as.data.frame(pct(id, estacion_tipo, clase)))
print(as.data.frame(pct(id, region, estacion_tipo)))
cat("\n--- 87585 ---\n"); print(as.data.frame(res %>% filter(station_id=="87585") %>% select(tipo,variable_x,variable_y,familia,AIC,n,tau_id,tau_im)))
cat("\n--- n eventos por tipo ---\n"); print(summary(id$n)); print(id %>% group_by(tipo) %>% summarise(min=min(n), med=median(n), max=max(n)))
print(table(id$familia, id$n<30))
