# Tests de los filtros de entrada de AjustarCopulas (lib/funciones_worker.R)
suppressMessages({library(dplyr); library(purrr); library(tibble)})
raiz <- if (file.exists("lib")) "." else "../.."
src <- readLines(file.path(raiz, "lib/funciones_worker.R"), encoding = "latin1")
ini <- grep("^AjustarCopulas <- function", src); fin <- grep("^[A-Za-z.]+ <- function", src)
fin <- min(fin[fin > ini]) - 1
IdentificarIdColumn <- function(x) "station_id"
eval(parse(text = src[ini:fin]))
script <- list(info = function(...) invisible(NULL))
fila <- function(...) tibble(station_id = "1", nombre = "E", variable_x = "intensidad", variable_y = "duracion",
                             tipo_serie = "perturbada", n_serie = 1, familia = "frank", funcion_ajuste = "AjustarCopulaFrank",
                             x_es_estacionaria = TRUE, y_es_estacionaria = TRUE, son_dependientes = TRUE,
                             observada_valida = TRUE, ...)
sale_sin_ajuste <- function(f, ...) {
  r <- AjustarCopulas(f, script, NULL, 0.05, ...)
  is.na(r$copula[[1]]) && is.na(r$bondad.ajuste[[1]])
}
test_that("no se ajusta cuando la serie observada del par no es valida", {
  expect_true(sale_sin_ajuste(fila() %>% mutate(observada_valida = FALSE)))
})
test_that("no se ajusta si no son dependientes o no son estacionarias", {
  expect_true(sale_sin_ajuste(fila() %>% mutate(son_dependientes = FALSE)))
  expect_true(sale_sin_ajuste(fila() %>% mutate(x_es_estacionaria = FALSE)))
  expect_true(sale_sin_ajuste(fila() %>% mutate(y_es_estacionaria = FALSE)))
})
test_that("con exigir.serie.observada = FALSE el filtro de la observada no corta (falla mas adelante, no en el control)", {
  expect_error(AjustarCopulas(fila() %>% mutate(observada_valida = FALSE), script, NULL, 0.05, exigir.serie.observada = FALSE))
})
