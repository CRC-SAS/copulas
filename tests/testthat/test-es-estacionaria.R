# Tests de EsEstacionaria (lib/funciones_test_estacionaridad.R)
suppressMessages({library(dplyr); library(purrr); library(tibble)})
raiz <- if (file.exists("lib")) "." else "../.."
source(file.path(raiz, "lib/funciones_bondad_ajuste.R"), encoding = "latin1")  # ParametrosADataFrame
source(file.path(raiz, "lib/funciones_test_estacionaridad.R"), encoding = "latin1")

set.seed(2)
n <- 120
fechas <- seq(as.Date("1961-01-01"), by = "180 days", length.out = n)
x.est   <- rnorm(n)                          # ruido blanco: estacionaria
x.tend  <- 0.05 * seq_len(n) + rnorm(n)      # tendencia fuerte: no estacionaria
x.salto <- c(rnorm(n/2), rnorm(n/2, mean = 4)) # cambio de nivel: no estacionaria

test.si <- c("BoxPierceLjungBox", "PuntoCambioUni", "MannKendall")

test_that("una serie estacionaria pasa todos los tests", {
  r <- EsEstacionaria(x.est, fechas, test.si, 0.05)
  expect_equal(r$test.name, test.si)
  expect_true(all(r$pasa.test))
})

test_that("una serie con tendencia no es estacionaria", {
  r <- EsEstacionaria(x.tend, fechas, test.si, 0.05)
  expect_false(all(r$pasa.test))
})

test_that("una serie con cambio de nivel no es estacionaria", {
  r <- EsEstacionaria(x.salto, fechas, test.si, 0.05)
  expect_false(all(r$pasa.test))
})

test_that("se corta en el primer test que falla", {
  r <- EsEstacionaria(x.tend, fechas, c("MannKendall", "PuntoCambioUni"), 0.05)
  expect_equal(r$test.name, "MannKendall")
  expect_false(r$pasa.test)
})

test_that("umbral: con umbral 0 una serie estacionaria pasa, con umbral 1 no pasa", {
  expect_true(all(EsEstacionaria(x.est, fechas, "MannKendall", 0)$pasa.test))
  expect_false(all(EsEstacionaria(x.est, fechas, "MannKendall", 1)$pasa.test))
})

test_that("un test que falla por error cuenta como no estacionaria", {
  r <- EsEstacionaria(x.est, fechas, "TestQueNoExiste", 0.05)
  expect_false(all(r$pasa.test))
  expect_true(r$exec.error)
})

test_that("UmbralEstacionariedad aplica Bonferroni por defecto y permite desactivarlo", {
  expect_equal(UmbralEstacionariedad(0.05, 5), 0.01)
  expect_equal(UmbralEstacionariedad(0.05, 5, "bonferroni"), 0.01)
  expect_equal(UmbralEstacionariedad(0.05, 5, NULL), 0.01)
  expect_equal(UmbralEstacionariedad(0.05, 5, "ninguna"), 0.05)
  expect_equal(UmbralEstacionariedad(0.05, 1), 0.05)
  expect_error(UmbralEstacionariedad(0.05, 5, "holm"))
})

test_that("con umbral Bonferroni una serie con deriva leve deja de rechazarse", {
  # buscar (de forma reproducible) una serie con p-valor de Mann-Kendall entre 0,01 y 0,05
  pv <- function(y) ParametrosADataFrame(TestMannKendall(y, fechas)) %>% dplyr::filter(parametro == "p.value") %>% dplyr::pull(valor)
  y <- NULL
  for (semilla in 1:200) {
    set.seed(semilla); cand <- 0.004 * seq_len(n) + rnorm(n)
    if (pv(cand) > 0.01 && pv(cand) <= 0.05) { y <- cand; break }
  }
  expect_false(is.null(y))
  expect_false(all(EsEstacionaria(y, fechas, "MannKendall", 0.05)$pasa.test))
  expect_true(all(EsEstacionaria(y, fechas, "MannKendall", UmbralEstacionariedad(0.05, 5))$pasa.test))
})

test_that("ConsolidarEstacionariedad no da por estacionaria una serie observada que rechaza el test", {
  expect_false(ConsolidarEstacionariedad(FALSE))
  expect_true(ConsolidarEstacionariedad(TRUE))
  expect_true(ConsolidarEstacionariedad(c(TRUE, FALSE)))        # empate
  expect_false(ConsolidarEstacionariedad(c(TRUE, FALSE, FALSE)))
  expect_false(ConsolidarEstacionariedad(c(TRUE, TRUE, FALSE, FALSE, FALSE)))
  expect_true(ConsolidarEstacionariedad(c(TRUE, TRUE, TRUE, FALSE, FALSE)))
})
