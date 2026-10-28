# DeterminarMejorCopula desempata por BIC, MRE, RMSE y validacion cruzada (lib/funciones_mejor_ajuste_copula.R)
suppressMessages(library(dplyr))
raiz <- if (file.exists("lib")) "." else "../.."
source(file.path(raiz, "lib/funciones_mejor_ajuste_copula.R"), encoding = "UTF-8")
ajuste <- function(familia, aic, bic, mre = 1, rmse = 1, xv = 1)
  list(copula = list(familia = familia, copula = familia, parametro = 1, dependencia = 0.5),
       bondad.ajuste = list(pasa.tests = TRUE, estadisticos = data.frame(
         test = c("AIC", "BIC", "MRE", "RMSE", "Validacion cruzada"), parametro = "x", valor = c(aic, bic, mre, rmse, xv))))
test_that("con AIC igual gana el menor BIC", {
  r <- DeterminarMejorCopula(list(ajuste("a", -10, -5), ajuste("b", -10, -8)))
  expect_equal(r$familia, "b")
  r <- DeterminarMejorCopula(list(ajuste("b", -10, -8), ajuste("a", -10, -5)))
  expect_equal(r$familia, "b")
})
test_that("con AIC y BIC iguales gana el menor MRE", {
  r <- DeterminarMejorCopula(list(ajuste("a", -10, -8, mre = 3), ajuste("b", -10, -8, mre = 2)))
  expect_equal(r$familia, "b")
})
test_that("gana el menor AIC", {
  expect_equal(DeterminarMejorCopula(list(ajuste("a", -9, -9), ajuste("b", -10, -1)))$familia, "b")
})
