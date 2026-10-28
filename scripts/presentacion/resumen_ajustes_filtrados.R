# Tablas del par intensidad-duracion con los filtros corregidos (dependencia + estacionariedad Bonferroni).
# Cuenta pares con cópula ganadora y calcula % de cola solo sobre esos. Requiere resumen_ajustes.R (DIR_PRES/resumen_ajustes.csv)
# Se corre desde la raiz: Rscript scripts/presentacion/resumen_ajustes_filtrados.R
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
suppressMessages(library(dplyr)); suppressMessages(library(readr))
res <- read_csv(file.path(DIR_PRES, "resumen_ajustes.csv"), show_col_types=FALSE)
id <- res %>% filter(variable_x=="intensidad", variable_y=="duracion") %>% mutate(gana=!is.na(familia))
pct <- function(d, ...) d %>% group_by(...) %>% summarise(pares=n(), con_cop=sum(gana), pct_gana=round(100*mean(gana)),
  sin=round(100*mean(cola[gana]=="Sin cola")), sup=round(100*mean(cola[gana]=="Superior")), inf=round(100*mean(cola[gana]=="Inferior")), amb=round(100*mean(cola[gana]=="Ambas")),
  con_cola=round(100*mean(cola[gana]!="Sin cola")), n_med=median(n), tau_med=round(median(tau_id,na.rm=TRUE),2), .groups="drop")
cat("TOTAL\n"); print(as.data.frame(pct(id, todos="todos")))
for (g in list("tipo","clase","estacion_tipo",c("estacion_tipo","clase"),"region",c("region","estacion_tipo"))) { cat("\n", paste(g,collapse="x"), "\n"); print(as.data.frame(do.call(pct, c(list(id), lapply(g, as.name))))) }
cat("\nfamilias\n"); print(table(id$familia, id$clase, useNA="ifany"))
cat("\nn<30:", sum(id$n<30), " con cop:", sum(id$gana[id$n<30]), "\n")
write_csv(id, file.path(DIR_PRES, "id_filtrado.csv"))
