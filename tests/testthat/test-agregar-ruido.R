# AgregarRuido suma ruido uniforme en (0, valor minimo de deteccion) (lib/funciones_auxiliares.R)
raiz <- if (file.exists("lib")) "." else "../.."
src <- readLines(file.path(raiz, "lib/funciones_auxiliares.R"), encoding = "latin1")
ini <- grep("^AgregarRuido <- function", src); fin <- grep("^[A-Za-z.]+ <- function", src)
fin <- min(c(fin[fin > ini], length(src) + 1)) - 1
eval(parse(text = src[ini:fin]))

test_that("el ruido esta en (0, v) para cualquier valor minimo de deteccion", {
  set.seed(1)
  for (v in c(1, 0.1)) {
    r <- AgregarRuido(rep(10, 5000), v) - 10
    expect_true(all(r >= 0 & r <= v))
    expect_gt(max(r), 0.95 * v)       # antes con v = 0.1 el maximo era 0.01
    expect_equal(mean(r), v / 2, tolerance = 0.05 * v)
  }
})
test_that("rompe los empates", {
  expect_equal(length(unique(AgregarRuido(rep(3, 100), 1))), 100)
})
