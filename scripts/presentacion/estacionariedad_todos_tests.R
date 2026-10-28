# P-valores de los 5 tests de estacionariedad en todas las series (observadas y perturbadas) de los 6 tipos de evento
# Se necesitan porque el pipeline aplica los tests en secuencia y corta al primero que rechaza (EsEstacionaria): en los
# resultados guardados (estacionariedad_<tipo>.rds) faltan los p-valores de los tests posteriores al primer rechazo.
# Aca se calculan los cinco en cada serie, con las mismas funciones y la misma semilla por tarea que el pipeline
# (FijarSemillaTarea), asi que los p-valores coinciden con los del pipeline en los tests que este llego a ejecutar.
# Se corre desde la raiz del repo: Rscript scripts/presentacion/estacionariedad_todos_tests.R
# Entradas: data/partial/eventos_completos_<tipo>.rds (series observadas y perturbadas de cada corrida)
# Variables de entorno opcionales: DIR_PRES (salida), SEMILLA (s.series.perturbadas, 123), NCORES (nucleos, 4),
#   ESTACIONES (ids separados por coma, para probar con pocas estaciones)
# Salidas: DIR_PRES/est_todos_tests.rds (columnas: station_id, variable, tipo_serie, n_serie, n, test, p, tipo), que usa figuras_filtros.R
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
semilla  <- as.integer(Sys.getenv("SEMILLA", "123"))
ncores   <- as.integer(Sys.getenv("NCORES", "4"))
estaciones <- Sys.getenv("ESTACIONES", "")
dir.create(DIR_PRES, showWarnings = FALSE, recursive = TRUE)
suppressMessages({library(dplyr); library(purrr); library(tibble); library(xts); library(parallel)})
for (f in c("lib/funciones_auxiliares.R", "lib/funciones_bondad_ajuste.R", "lib/funciones_test_estacionaridad.R"))
  source(f, encoding = "latin1")   # los archivos de lib mezclan Latin-1 y UTF-8 en los comentarios
tests <- c("BoxPierceLjungBox", "PuntoCambioUni", "PuntoCambioAutocopula", "MannKendall", "EstacionaridadEmpirico")
tipos <- c("olacalor", "olafrio", "periodocalidotn", "periodocalidotx", "periodofriotn", "periodofriotx")

p_valores <- function(x, fechas, clave) {
  FijarSemillaTarea(semilla, "estacionariedad", clave$station_id, clave$variable, clave$tipo_serie, clave$n_serie, clave$realizacion)
  purrr::map_dfr(tests, function(t) {
    p <- tryCatch(unname(do.call(paste0("Test", t), list(x = x, fechas = fechas))$p.value), error = function(e) NA_real_)
    tibble(test = t, p = p)
  })
}

res <- list()
for (tp in tipos) {
  ev <- readRDS(sprintf("data/partial/eventos_completos_%s.rds", tp))
  if (nzchar(estaciones)) ev <- ev %>% filter(station_id %in% strsplit(estaciones, ",")[[1]])
  grupos <- ev %>% group_by(station_id, variable, tipo_serie, n_serie, realizacion) %>% group_split()
  cat(sprintf("%s: %d series\n", tp, length(grupos)))
  res[[tp]] <- parallel::mclapply(grupos, function(g) {
    clave <- g[1, c("station_id", "variable", "tipo_serie", "n_serie", "realizacion")]
    p_valores(g$valor, g$fecha_inicio, clave) %>%
      mutate(station_id = clave$station_id, variable = clave$variable, tipo_serie = clave$tipo_serie, n_serie = clave$n_serie, n = nrow(g), tipo = tp)
  }, mc.cores = ncores) %>% bind_rows()
}
out <- bind_rows(res) %>% select(station_id, variable, tipo_serie, n_serie, n, test, p, tipo) %>%
  arrange(tipo, station_id, variable, tipo_serie, n_serie, match(test, tests))
saveRDS(out, file.path(DIR_PRES, "est_todos_tests.rds"))
cat(sprintf("%d filas; rechazos con p <= 0,01 por test (series observadas):\n", nrow(out)))
print(as.data.frame(out %>% filter(tipo_serie == "observada") %>% group_by(test) %>% summarise(rechaza = mean(p <= 0.01, na.rm = TRUE))))
