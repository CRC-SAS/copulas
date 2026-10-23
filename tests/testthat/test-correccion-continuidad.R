# Test de la correccion de continuidad en CalcularGrillaPeriodoRetorno (lib/funciones_periodo_retorno.R)
suppressMessages({library(copula); library(dplyr); library(tidyr)})
raiz <- if (file.exists("lib")) "." else "../.."
src <- readLines(file.path(raiz, "lib/funciones_periodo_retorno.R"), encoding = "latin1")
ini <- grep("^CalcularGrillaPeriodoRetorno <- function", src); fin <- grep("^GraficarPeriodoRetorno <- function", src)
eval(parse(text = src[ini:(fin - 1)]))
ini <- grep("^CalcularGrillaPeriodoRetornoUV <- function", src); fin <- grep("^GraficarPeriodoRetornoUV <- function", src)
eval(parse(text = src[ini:(fin - 1)]))
mv <- mvdc(normalCopula(0.3), c("gamma", "gamma"), list(list(shape = 2, rate = 0.5), list(shape = 4, rate = 0.25)))
test_that("desplazar_x = 0 no cambia el resultado y desplazar_x = 0.5 evalua en x - 0.5", {
  gx <- c(3, 4, 5); gy <- c(2, 4)
  g0 <- CalcularGrillaPeriodoRetorno(mv, 60, 100, gx, gy)
  g1 <- CalcularGrillaPeriodoRetorno(mv, 60, 100, gx, gy, desplazar_x = 0.5)
  gref <- CalcularGrillaPeriodoRetorno(mv, 60, 100, gx - 0.5, gy)
  expect_equal(g1$x, g0$x)                       # los valores mostrados no cambian
  expect_equal(g1$T, gref$T)                     # se evalua en x - 0.5
  expect_true(all(g1$T <= g0$T))                 # P(D > k - 0.5) >= P(D > k): al evaluar mas a la izquierda el periodo de retorno no sube
})
test_that("el caso univariado aplica el mismo desplazamiento", {
  g0 <- CalcularGrillaPeriodoRetornoUV("gamma", list(shape = 2, rate = 0.5), 60, 100, c(3, 4, 5))
  g1 <- CalcularGrillaPeriodoRetornoUV("gamma", list(shape = 2, rate = 0.5), 60, 100, c(3, 4, 5), desplazar_x = 0.5)
  expect_equal(g1$T, CalcularGrillaPeriodoRetornoUV("gamma", list(shape = 2, rate = 0.5), 60, 100, c(2.5, 3.5, 4.5))$T)
  expect_true(all(g1$T <= g0$T))
})
