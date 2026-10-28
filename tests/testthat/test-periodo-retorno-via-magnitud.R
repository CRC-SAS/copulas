# Test de ProbabilidadConjuntaViaMagnitud (lib/funciones_periodo_retorno.R): la simulacion debe coincidir con la integral exacta
suppressMessages({library(copula)})
raiz <- if (file.exists("lib")) "." else "../.."
src <- readLines(file.path(raiz, "lib/funciones_periodo_retorno.R"), encoding = "latin1")
ini <- grep("^ProbabilidadConjuntaViaMagnitud <- function", src); fin <- grep("^PeriodoRetornoViaMagnitud <- function", src)
eval(parse(text = src[ini:(fin - 1)])); eval(parse(text = src[fin:length(src)]))
pgamma_ <- stats::pgamma; qgamma_ <- stats::qgamma
mv <- mvdc(normalCopula(0.6), c("gamma", "gamma"), list(list(shape = 2, rate = 0.5), list(shape = 4, rate = 0.25)))
# P(D > d, M > i * D) = int_{d}^{inf} P(M > i t | D = t) f_D(t) dt; con copula normal P(V <= v | U = u) = pnorm((qnorm(v) - rho qnorm(u)) / sqrt(1 - rho^2))
exacta <- function(d, i, rho = 0.6) {
  f <- function(t) { u <- pmin(pgamma(t, 2, 0.5), 1 - 1e-12); v <- pmin(pgamma(i * t, 4, 0.25), 1 - 1e-12)
    (1 - pnorm((qnorm(v) - rho * qnorm(u)) / sqrt(1 - rho^2))) * dgamma(t, 2, 0.5) }
  stats::integrate(f, d, Inf)$value
}
test_that("la simulacion coincide con la integral exacta", {
  for (pt in list(c(2, 2), c(4, 3), c(6, 5))) {
    p <- ProbabilidadConjuntaViaMagnitud(mv, pt[1], pt[2], nsim = 600000)
    expect_equal(p, exacta(pt[1], pt[2]), tolerance = 0.03)
  }
})
test_that("es reproducible y deja intacta la semilla global", {
  set.seed(1); a <- runif(1); set.seed(1)
  p1 <- ProbabilidadConjuntaViaMagnitud(mv, 3, 3, nsim = 20000); p2 <- ProbabilidadConjuntaViaMagnitud(mv, 3, 3, nsim = 20000)
  expect_equal(p1, p2); expect_equal(runif(1), a)
})
test_that("el periodo de retorno es N / (n * P)", {
  p <- ProbabilidadConjuntaViaMagnitud(mv, 3, 3, nsim = 20000)
  expect_equal(PeriodoRetornoViaMagnitud(mv, 3, 3, N = 60, n = 100, nsim = 20000), 60 / (100 * p))
})
