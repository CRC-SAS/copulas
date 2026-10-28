# FijarSemillaTarea: semilla determinista por tarea (lib/funciones_auxiliares.R)
raiz <- if (file.exists("lib")) "." else "../.."
src <- readLines(file.path(raiz, "lib/funciones_auxiliares.R"), encoding = "latin1")
ini <- grep("^FijarSemillaTarea <- function", src); fin <- grep("^[A-Za-z.]+ <- function", src)
fin <- min(c(fin[fin > ini], length(src) + 1)) - 1
eval(parse(text = src[ini:fin]))

test_that("la misma clave da los mismos numeros aleatorios", {
  FijarSemillaTarea(123, "copula", "87585", "intensidad", "duracion", "observada", 0, "gumbel")
  a <- runif(5)
  FijarSemillaTarea(123, "copula", "87585", "intensidad", "duracion", "observada", 0, "gumbel")
  expect_equal(runif(5), a)
})
test_that("claves o semillas distintas dan numeros distintos", {
  FijarSemillaTarea(123, "copula", "87585", "intensidad", "duracion", "observada", 0, "gumbel"); a <- runif(5)
  FijarSemillaTarea(123, "copula", "87585", "intensidad", "duracion", "observada", 0, "joe"); expect_false(isTRUE(all.equal(runif(5), a)))
  FijarSemillaTarea(124, "copula", "87585", "intensidad", "duracion", "observada", 0, "gumbel"); expect_false(isTRUE(all.equal(runif(5), a)))
})
test_that("con semilla NULL no se toca el generador", {
  set.seed(7); a <- runif(3)
  set.seed(7); FijarSemillaTarea(NULL, "x"); expect_equal(runif(3), a)
})
