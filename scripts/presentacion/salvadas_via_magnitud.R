# Estaciones "salvadas" por la via duracion-magnitud: tienen copula duracion-magnitud pero no intensidad-duracion (olas de calor y de frio)
# Se corre desde la raiz del repo, despues de resumen_ajustes_filtrados.R y periodo_retorno_via_magnitud.R
# Entradas: DIR_PRES/id_filtrado.csv (pares estacion x tipo y su copula intensidad-duracion), DIR_PRES/via_magnitud_validacion.csv (pares con copula duracion-magnitud)
# Salidas: DIR_PRES/salvadas_via_magnitud.csv y DIR_PRES/salvadas_via_magnitud.txt
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
suppressMessages(library(dplyr))
sink(file.path(DIR_PRES, "salvadas_via_magnitud.txt"), split = TRUE)
id <- read.csv(file.path(DIR_PRES, "id_filtrado.csv")) %>% mutate(station_id = as.character(station_id))
dm <- readr::read_csv(file.path(DIR_PRES, "via_magnitud_validacion.csv"), show_col_types = FALSE) %>% distinct(slug, station_id) %>% mutate(station_id = as.character(station_id), dm = TRUE)
d <- id %>% left_join(dm, by = c("slug", "station_id")) %>%
  mutate(di = !is.na(familia), dm = !is.na(dm), solo_dm = dm & !di, alguna = di | dm)
stopifnot(nrow(d) == nrow(id), sum(dm$dm) == sum(d$dm))
res <- function(x) x %>% summarise(estaciones = n(), con_ID = sum(di), solo_DM = sum(solo_dm), solo_ID = sum(di & !dm), con_ID_o_DM = sum(alguna), .groups = "drop")
cat("== Totales por tipo\n"); print(as.data.frame(d %>% group_by(tipo = slug) %>% res()))
cat("\n== Olas por region\n"); o <- d %>% filter(grepl("^ola", slug)) %>% mutate(tipo = ifelse(slug == "olacalor", "Calor", "Frio"))
print(as.data.frame(o %>% group_by(tipo, region) %>% res() %>% arrange(tipo, desc(con_ID_o_DM / estaciones))))
cat("\n== Periodos (todos los tipos) por region\n"); print(as.data.frame(d %>% filter(!grepl("^ola", slug)) %>% group_by(region) %>% res()))
write.csv(d %>% select(station_id, slug, region, di, dm, solo_dm, alguna), file.path(DIR_PRES, "salvadas_via_magnitud.csv"), row.names = FALSE)
sink()
