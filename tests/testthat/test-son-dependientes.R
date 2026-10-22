# Tests de SonDependientes (lib/funciones_test_independencia.R)
# Correr desde la raiz del repo: Rscript -e 'testthat::test_dir("tests/testthat")'

suppressMessages({library(dplyr); library(purrr); library(tibble)})
raiz <- if (file.exists("lib")) "." else "../.."
source(file.path(raiz, "lib/funciones_bondad_ajuste.R"), encoding = "latin1")  # define ParametrosADataFrame
source(file.path(raiz, "lib/funciones_test_independencia.R"), encoding = "latin1")

tests.dep <- c("Kendall", "Spearman", "CopulaEmpirica")

set.seed(1)
n <- 80
x.indep <- runif(n); y.indep <- runif(n)
x.dep   <- runif(n); y.dep   <- x.dep + rnorm(n, sd = 0.1)

test_that("devuelve un data.frame y no una lista de resultados parciales", {
  r <- SonDependientes(x.dep, y.dep, tests.dep, 0.05)
  expect_s3_class(r, "data.frame")
  expect_named(r, c("pasa.test", "test.name", "p.value", "exec.error"))
})

test_that("series dependientes pasan los tres tests", {
  r <- SonDependientes(x.dep, y.dep, tests.dep, 0.05)
  expect_equal(r$test.name, tests.dep)
  expect_true(all(r$pasa.test))
})

test_that("series independientes no pasan y el resultado es FALSE", {
  r <- SonDependientes(x.indep, y.indep, tests.dep, 0.05)
  expect_false(all(r$pasa.test))
  expect_false(dplyr::last(r$pasa.test))
})

test_that("se corta en el primer test que falla (no corre los siguientes)", {
  r <- SonDependientes(x.indep, y.indep, tests.dep, 0.05)
  expect_equal(r$test.name, "Kendall")
  expect_equal(nrow(r), 1)
})

test_that("umbral de p-valor: con umbral 1 todo pasa, con umbral 0 nada pasa", {
  expect_false(all(SonDependientes(x.dep, y.dep, c("Kendall"), 0)$pasa.test))
  expect_true(all(SonDependientes(x.indep, y.indep, c("Kendall"), 1)$pasa.test))
})

test_that("un test que falla por error cuenta como no dependiente", {
  r <- SonDependientes(x.dep, y.dep, c("Kendall", "TestQueNoExiste"), 0.05)
  expect_false(all(r$pasa.test))
  expect_true(dplyr::last(r$exec.error))
})

test_that("la expresion del worker (all(...$pasa.test)) refleja el resultado", {
  expect_true(all(SonDependientes(x.dep, y.dep, tests.dep, 0.05)$pasa.test))
  expect_false(all(SonDependientes(x.indep, y.indep, tests.dep, 0.05)$pasa.test))
})

test_that("TestIndCopulaEmpirica es reproducible y no altera la semilla global", {
  set.seed(3); x <- rnorm(60); y <- 0.2 * x + rnorm(60)
  set.seed(10); a <- TestIndCopulaEmpirica(x, y); s1 <- runif(1)
  set.seed(10); b <- TestIndCopulaEmpirica(x, y); s2 <- runif(1)
  set.seed(99); d <- TestIndCopulaEmpirica(x, y)
  expect_equal(a$p.value, b$p.value)
  expect_equal(a$p.value, d$p.value)
  set.seed(10); invisible(TestIndCopulaEmpirica(x, y)); expect_equal(runif(1), s1)
})
