# El BIC penaliza con ln(n), n = cantidad de eventos (lib/funciones_bondad_ajuste_copulas.R)
raiz <- if (file.exists("lib")) "." else "../.."
src <- readLines(file.path(raiz, "lib/funciones_bondad_ajuste_copulas.R"), encoding = "latin1")
ini <- grep("^TestBIC <- function", src); fin <- grep("^[A-Za-z.]+ <- function|^# --- Consolidacion", src)
fin <- min(fin[fin > ini]) - 1
eval(parse(text = src[ini:fin]))

test_that("TestBIC usa ln(n eventos) y no ln(cantidad de variables)", {
  x <- c(1.2, 3.4, 2.2, 5.1, 4.4, 3.9, 2.8, 6.0, 1.9, 4.1)
  y <- c(2.0, 3.1, 2.5, 5.5, 4.0, 4.4, 2.2, 5.8, 1.5, 4.9)
  cop <- list(parametro = list(theta = 2), bondad.ajuste = list(loglike = 3.5))
  esperado <- -2 * 3.5 + log(length(x)) * 1
  expect_equal(TestBIC(copula = cop, x = x, y = y), esperado)
  # con mas eventos la penalizacion crece
  expect_gt(TestBIC(copula = cop, x = rep(x, 3), y = rep(y, 3)), TestBIC(copula = cop, x = x, y = y))
})
