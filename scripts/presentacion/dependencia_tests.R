# Cuantos pares estacion x tipo pasan cada test de dependencia (Kendall, Spearman, copula empirica), par intensidad-duracion, serie observada
# Los tests se aplican en secuencia: el primero que rechaza corta (no se calculan los siguientes)
# Se corre desde la raiz del repo. Entradas: data/partial/dependencia_<tipo>.rds. Salida: DIR_PRES/dependencia_tests.csv y .txt
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
suppressMessages({library(dplyr); library(purrr)})
sink(file.path(DIR_PRES, "dependencia_tests.txt"), split = TRUE)
tipos <- c("olacalor", "olafrio", "periodocalidotx", "periodocalidotn", "periodofriotx", "periodofriotn")
r <- bind_rows(lapply(tipos, function(t) readRDS(sprintf("data/partial/dependencia_%s.rds", t)) %>%
  filter(tipo_serie == "observada", variable_x == "intensidad", variable_y == "duracion", as.character(station_id) != "87420") %>%
  transmute(tipo = t, station_id, grupo = ifelse(grepl("^ola", t), "Olas", "Períodos"), det = detalles_dependencia)))
res <- function(l) { p <- l$pasa.test; n <- l$test.name
  c(kendall = isTRUE(p[n == "Kendall"]), spearman = isTRUE(p[n == "Spearman"]), empirica = isTRUE(p[n == "CopulaEmpirica"]),
    llega_spearman = "Spearman" %in% n, llega_empirica = "CopulaEmpirica" %in% n, error = any(l$exec.error)) }
m <- as.data.frame(t(sapply(r$det, res))); r <- bind_cols(r %>% select(-det), m)
r$todos <- r$kendall & r$spearman & r$empirica
r$primer_rechazo <- ifelse(r$todos, "ninguno", ifelse(!r$kendall, "Kendall", ifelse(!r$spearman, "Spearman", "Copula empirica")))
f <- function(x) x %>% summarise(pares = n(), pasa_kendall = sum(kendall), pasa_spearman = sum(spearman), pasa_empirica = sum(empirica), pasan_los_tres = sum(todos),
  rechaza_kendall = sum(primer_rechazo == "Kendall"), rechaza_spearman = sum(primer_rechazo == "Spearman"), rechaza_empirica = sum(primer_rechazo == "Copula empirica"), errores = sum(error), .groups = "drop")
cat("== Total\n"); print(as.data.frame(f(r %>% mutate(g = "todos") %>% group_by(g))))
cat("\n== Por grupo\n"); print(as.data.frame(f(r %>% group_by(grupo))))
cat("\n== Por tipo\n"); print(as.data.frame(f(r %>% group_by(tipo))))
write.csv(f(r %>% group_by(tipo)), file.path(DIR_PRES, "dependencia_tests.csv"), row.names = FALSE)
sink()
