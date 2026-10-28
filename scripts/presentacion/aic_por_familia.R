# AIC por familia de copula (mediana sobre series perturbadas) para las olas de frio. Genera DIR_PRES/aic_olafrio.rds
# Se corre desde la raiz del repo: Rscript scripts/presentacion/aic_por_familia.R
# Entradas: data/output/*.rds (corrida Olas_20), data/input/eventos_identificados_frio_calor.csv
# Salidas: DIR_PRES (por defecto data/output/presentacion); las figuras van en DIR_PRES/figuras
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
dir.create(file.path(DIR_PRES, "figuras"), recursive = TRUE, showWarnings = FALSE)
REGIONES_CSV <- Sys.getenv("REGIONES_CSV", "data/input/regiones.csv")
suppressMessages(library(dplyr))
a <- readRDS("data/output/ajustes_multivariados_olafrio.rds") %>% filter(variable_x=="intensidad", variable_y=="duracion")
m <- readRDS("data/output/mejores_ajustes_multivariados_olafrio.rds") %>% filter(variable_x=="intensidad", variable_y=="duracion")
ev <- readr::read_csv("data/input/eventos_identificados_frio_calor.csv", show_col_types=FALSE) %>% filter(duracion>=3, tipo_evento=="OlaFrioTXyTN") %>% mutate(station_id=as.character(station_id)) %>% count(station_id, name="n")
get <- function(i) { e <- a$bondad.ajuste[[i]]$estadisticos; if (is.null(e) || !"test" %in% names(e)) return(NULL)
  e <- e[e$test=="AIC", , drop=FALSE]; if (!nrow(e)) return(NULL); data.frame(station_id=as.character(a$station_id[i]), familia=a$familia[i], aic=e$valor) }
r <- bind_rows(lapply(seq_len(nrow(a)), get)) %>% group_by(station_id, familia) %>% summarise(aic=median(aic), .groups="drop")
w <- r %>% group_by(station_id) %>% arrange(aic, .by_group=TRUE) %>% mutate(rk=row_number()) %>% summarise(n_fam=n(), a1=aic[1], f1=familia[1], a2=aic[2], f2=familia[2], dif=a2-a1, .groups="drop") %>% left_join(ev, by="station_id") %>% left_join(m %>% mutate(station_id=as.character(station_id)) %>% select(station_id, ganadora=familia), by="station_id")
print(as.data.frame(w %>% filter(n<=26) %>% arrange(n)))
saveRDS(r, file.path(DIR_PRES, "aic_olafrio.rds"))
