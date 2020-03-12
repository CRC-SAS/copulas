# ---------------------------------------------------------------------------- #
# AEP 4
raep4 <- function(n, xi, alpha, kappa, h) {
  obj.paraep4 <- list(
    type = "aep4",
    para = c("xi" = xi, "alpha" = alpha, "kappa" = kappa, "h" = h),
    source = "paraep4"
  )
  return (lmomco::rlmomco(n, para = obj.paraep4))
}
daep4 <- function(x, xi, alpha, kappa, h) {
  obj.paraep4 <- list(
    type = "aep4",
    para = c("xi" = xi, "alpha" = alpha, "kappa" = kappa, "h" = h),
    source = "paraep4"
  )
  return (lmomco::dlmomco(x, para = obj.paraep4))
}
paep4 <- function(q, xi, alpha, kappa, h) {
  obj.paraep4 <- list(
    type = "aep4",
    para = c("xi" = xi, "alpha" = alpha, "kappa" = kappa, "h" = h),
    source = "paraep4"
  )
  return (lmomco::plmomco(q, para = obj.paraep4))
}
qaep4 <- function(f, xi, alpha, kappa, h) {
  obj.paraep4 <- list(
    type = "aep4",
    para = c("xi" = xi, "alpha" = alpha, "kappa" = kappa, "h" = h),
    source = "paraep4"
  )
  return (lmomco::qlmomco(f, para = obj.paraep4))
}
# ---------------------------------------------------------------------------- #
# Cauchy
rcau <- function(n, xi, alpha) {
  obj.parcau <- list(
    type = "cau",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parcau"
  )
  return (lmomco::rlmomco(n, para = obj.parcau))
}
dcau <- function(x, xi, alpha) {
  obj.parcau <- list(
    type = "cau",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parcau"
  )
  return (lmomco::dlmomco(x, para = obj.parcau))
}
pcau <- function(q, xi, alpha) {
  obj.parcau <- list(
    type = "cau",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parcau"
  )
  return (lmomco::plmomco(q, para = obj.parcau))
}
qcau <- function(f, xi, alpha) {
  obj.parcau <- list(
    type = "cau",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parcau"
  )
  return (lmomco::qlmomco(f, para = obj.parcau))
}
# ---------------------------------------------------------------------------- #
# Exponential
rexl <- function(n, xi, alpha) {
  obj.parexp <- list(
    type = "exp",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parexp"
  ) 
  return (lmomco::rlmomco(n, para = obj.parexp))
}
dexl <- function(x, xi, alpha) {
  obj.parexp <- list(
    type = "exp",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parexp"
  )
  return (lmomco::dlmomco(x, para = obj.parexp))
}
pexl <- function(q, xi, alpha) {
  obj.parexp <- list(
    type = "exp",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parexp"
  )
  return (lmomco::plmomco(q, para = obj.parexp))
}
qexl <- function(f, xi, alpha) {
  obj.parexp <- list(
    type = "exp",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parexp"
  )
  return (lmomco::qlmomco(f, para = obj.parexp))
}
# ---------------------------------------------------------------------------- #
# Gamma
rgam <- function(n, alpha, beta) {
  obj.pargam <- list(
    type = "gam",
    para = c("alpha" = alpha, 'beta' = beta),
    source = "pargam"
  )
  return (lmomco::rlmomco(n, para = obj.pargam))
}
dgam <- function(x, alpha, beta) {
  obj.pargam <- list(
    type = "gam",
    para = c("alpha" = alpha, 'beta' = beta),
    source = "pargam"
  )
  return (lmomco::dlmomco(x, para = obj.pargam))
}
pgam <- function(q, alpha, beta) {
  obj.pargam <- list(
    type = "gam",
    para = c("alpha" = alpha, 'beta' = beta),
    source = "pargam"
  )
  return (lmomco::plmomco(q, para = obj.pargam))
}
qgam <- function(f, alpha, beta) {
  obj.pargam <- list(
    type = "gam",
    para = c("alpha" = alpha, 'beta' = beta),
    source = "pargam"
  )
  return (lmomco::qlmomco(f, para = obj.pargam))
}
# ---------------------------------------------------------------------------- #
# Normal
rnor <- function(n, mu, sigma) {
  obj.parnor <- list(
    type = "nor",
    para = c('mu' = mu, "sigma"= sigma),
    source = "parnor"
  )
  return (lmomco::rlmomco(n, para = obj.parnor))
}
dnor <- function(x, mu, sigma) {
  obj.parnor <- list(
    type = "nor",
    para = c('mu' = mu, "sigma"= sigma),
    source = "parnor"
  )
  return (lmomco::dlmomco(x, para = obj.parnor))
}
pnor <- function(q, mu, sigma) {
  obj.parnor <- list(
    type = "nor",
    para = c('mu' = mu, "sigma"= sigma),
    source = "parnor"
  )
  return (lmomco::plmomco(q, para = obj.parnor))
}
qnor <- function(f, mu, sigma) {
  obj.parnor <- list(
    type = "nor",
    para = c('mu' = mu, "sigma"= sigma),
    source = "parnor"
  )
  return (lmomco::qlmomco(f, para = obj.parnor))
}
# ---------------------------------------------------------------------------- #
# Generalized Extreme Value
rgev <- function(n, xi, alpha, kappa) {
  obj.pargev <- list(
    type = "gev",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "pargev"
  )
  return (lmomco::rlmomco(n, para = obj.pargev))
}
dgev <- function(x, xi, alpha, kappa) {
  obj.pargev <- list(
    type = "gev",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "pargev"
  )
  return (lmomco::dlmomco(x, para = obj.pargev))
}
pgev <- function(q, xi, alpha, kappa) {
  obj.pargev <- list(
    type = "gev",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "pargev"
  )
  return (lmomco::plmomco(q, para = obj.pargev))
}
qgev <- function(f, xi, alpha, kappa) {
  obj.pargev <- list(
    type = "gev",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "pargev"
  )
  return (lmomco::qlmomco(f, para = obj.pargev))
}
# ---------------------------------------------------------------------------- #
# Generalized Lambda
rgld <- function(n, xi, alpha, kappa, h) {
  obj.pargld <- list(
    type = "gld",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa, "h" = h),
    source = "pargld"
  )
  return (lmomco::rlmomco(n, para = obj.pargld))
}
dgld <- function(x, xi, alpha, kappa, h) {
  obj.pargld <- list(
    type = "gld",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa, "h" = h),
    source = "pargld"
  )
  return (lmomco::dlmomco(x, para = obj.pargld))
}
pgld <- function(q, xi, alpha, kappa, h) {
  obj.pargld <- list(
    type = "gld",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa, "h" = h),
    source = "pargld"
  )
  return (lmomco::plmomco(q, para = obj.pargld))
}
qgld <- function(f, xi, alpha, kappa, h) {
  obj.pargld <- list(
    type = "gld",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa, "h" = h),
    source = "pargld"
  )
  return (lmomco::qlmomco(f, para = obj.pargld))
}
# ---------------------------------------------------------------------------- #
# Generalized Logistic
rglo <- function(n, xi, alpha, kappa) {
  obj.parglo <- list(
    type = "glo",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "parglo"
  )
  return (lmomco::rlmomco(n, para = obj.parglo))
}
dglo <- function(x, xi, alpha, kappa) {
  obj.parglo <- list(
    type = "glo",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "parglo"
  )
  return (lmomco::dlmomco(x, para = obj.parglo))
}
pglo <- function(q, xi, alpha, kappa) {
  obj.parglo <- list(
    type = "glo",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "parglo"
  )
  return (lmomco::plmomco(q, para = obj.parglo))
}
qglo <- function(f, xi, alpha, kappa) {
  obj.parglo <- list(
    type = "glo",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "parglo"
  )
  return (lmomco::qlmomco(f, para = obj.parglo))
}
# ---------------------------------------------------------------------------- #
# Generalized Normal
rgno <- function(n, xi, alpha, kappa) {
  obj.pargno <- list(
    type = "gno",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "pargno"
  )
  return (lmomco::rlmomco(n, para = obj.pargno))
}
dgno <- function(x, xi, alpha, kappa) {
  obj.pargno <- list(
    type = "gno",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "pargno"
  )
  return (lmomco::dlmomco(x, para = obj.pargno))
}
pgno <- function(q, xi, alpha, kappa) {
  obj.pargno <- list(
    type = "gno",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "pargno"
  )
  return (lmomco::plmomco(q, para = obj.pargno))
}
qgno <- function(f, xi, alpha, kappa) {
  obj.pargno <- list(
    type = "gno",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "pargno"
  )
  return (lmomco::qlmomco(f, para = obj.pargno))
}
# ---------------------------------------------------------------------------- #
# Govindarajulu
rgov <- function(n, xi, alpha, beta) {
  obj.pargov <- list(
    type = "gov",
    para = c('xi' = xi, "alpha" = alpha, "beta"= beta),
    source = "pargov"
  )
  return (lmomco::rlmomco(n, para = obj.pargov))
}
dgov <- function(x, xi, alpha, beta) {
  obj.pargov <- list(
    type = "gov",
    para = c('xi' = xi, "alpha" = alpha, "beta"= beta),
    source = "pargov"
  )
  return (lmomco::dlmomco(x, para = obj.pargov))
}
pgov <- function(q, xi, alpha, beta) {
  obj.pargov <- list(
    type = "gov",
    para = c('xi' = xi, "alpha" = alpha, "beta"= beta),
    source = "pargov"
  )
  return (lmomco::plmomco(q, para = obj.pargov))
}
qgov <- function(f, xi, alpha, beta) {
  obj.pargov <- list(
    type = "gov",
    para = c('xi' = xi, "alpha" = alpha, "beta"= beta),
    source = "pargov"
  )
  return (lmomco::qlmomco(f, para = obj.pargov))
}
# ---------------------------------------------------------------------------- #
# Generalized Pareto
rgpa <- function(n, xi, alpha, kappa) {
  obj.pargpa <- list(
    type = "gpa",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "pargpa"
  )
  return (lmomco::rlmomco(n, para = obj.pargpa))
}
dgpa <- function(x, xi, alpha, kappa) {
  obj.pargpa <- list(
    type = "gpa",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "pargpa"
  )
  return (lmomco::dlmomco(x, para = obj.pargpa))
}
pgpa <- function(q, xi, alpha, kappa) {
  obj.pargpa <- list(
    type = "gpa",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "pargpa"
  )
  return (lmomco::plmomco(q, para = obj.pargpa))
}
qgpa <- function(f, xi, alpha, kappa) {
  obj.pargpa <- list(
    type = "gpa",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa),
    source = "pargpa"
  )
  return (lmomco::qlmomco(f, para = obj.pargpa))
}
# ---------------------------------------------------------------------------- #
# Gumbel
rgum <- function(n, xi, alpha) {
  obj.pargum <- list(
    type = "gum",
    para = c('xi' = xi, "alpha" = alpha),
    source = "pargum"
  )
  return (lmomco::rlmomco(n, para = obj.pargum))
}
dgum <- function(x, xi, alpha) {
  obj.pargum <- list(
    type = "gum",
    para = c('xi' = xi, "alpha" = alpha),
    source = "pargum"
  )
  return (lmomco::dlmomco(x, para = obj.pargum))
}
pgum <- function(q, xi, alpha) {
  obj.pargum <- list(
    type = "gum",
    para = c('xi' = xi, "alpha" = alpha),
    source = "pargum"
  )
  return (lmomco::plmomco(q, para = obj.pargum))
}
qgum <- function(f, xi, alpha) {
  obj.pargum <- list(
    type = "gum",
    para = c('xi' = xi, "alpha" = alpha),
    source = "pargum"
  )
  return (lmomco::qlmomco(f, para = obj.pargum))
}
# ---------------------------------------------------------------------------- #
# Kappa
rkap <- function(n, xi, alpha, kappa, h) {
  obj.parkap <- list(
    type = "kap",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa, "h" = h),
    source = "parkap"
  )
  return (lmomco::rlmomco(n, para = obj.parkap))
}
dkap <- function(x, xi, alpha, kappa, h) {
  obj.parkap <- list(
    type = "kap",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa, "h" = h),
    source = "parkap"
  )
  return (lmomco::dlmomco(x, para = obj.parkap))
}
pkap <- function(q, xi, alpha, kappa, h) {
  obj.parkap <- list(
    type = "kap",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa, "h" = h),
    source = "parkap"
  )
  return (lmomco::plmomco(q, para = obj.parkap))
}
qkap <- function(f, xi, alpha, kappa, h) {
  obj.parkap <- list(
    type = "kap",
    para = c('xi' = xi, "alpha" = alpha, "kappa"= kappa, "h" = h),
    source = "parkap"
  )
  return (lmomco::qlmomco(f, para = obj.parkap))
}
# ---------------------------------------------------------------------------- #
# Kumaraswamy
rkur <- function(n, alpha, beta) {
  obj.parkur <- list(
    type = "kur",
    para = c("alpha" = alpha,"beta" = beta),
    source = "parkur"
  )
  return (lmomco::rlmomco(n, para = obj.parkur))
}
dkur <- function(x, alpha, beta) {
  obj.parkur <- list(
    type = "kur",
    para = c("alpha" = alpha,"beta" = beta),
    source = "parkur"
  )
  return (lmomco::dlmomco(x, para = obj.parkur))
}
pkur <- function(q, alpha, beta) {
  obj.parkur <- list(
    type = "kur",
    para = c("alpha" = alpha,"beta" = beta),
    source = "parkur"
  )
  return (lmomco::plmomco(q, para = obj.parkur))
}
qkur <- function(f, alpha, beta) {
  obj.parkur <- list(
    type = "kur",
    para = c("alpha" = alpha,"beta" = beta),
    source = "parkur"
  )
  return (lmomco::qlmomco(f, para = obj.parkur))
}
# ---------------------------------------------------------------------------- #
# Laplace
rlap <- function(n, xi, alpha) {
  obj.parlap <- list(
    type = "lap",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parlap"
  )
  return (lmomco::rlmomco(n, para = obj.parlap))
}
dlap <- function(x, xi, alpha) {
  obj.parlap <- list(
    type = "lap",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parlap"
  )
  return (lmomco::dlmomco(x, para = obj.parlap))
}
plap <- function(q, xi, alpha) {
  obj.parlap <- list(
    type = "lap",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parlap"
  )
  return (lmomco::plmomco(q, para = obj.parlap))
}
qlap <- function(f, xi, alpha) {
  obj.parlap <- list(
    type = "lap",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parlap"
  )
  return (lmomco::qlmomco(f, para = obj.parlap))
}
# ---------------------------------------------------------------------------- #
# Linear Mean Residual Quantile 
rlmrq <- function(n, mu, alpha) {
  obj.parlmrq <- list(
    type = "lmrq",
    para = c("mu" = mu, "alpha" = alpha),
    source = "parlmrq"
  )
  return (lmomco::rlmomco(n, para = obj.parlmrq))
}
dlmrq <- function(x, mu, alpha) {
  obj.parlmrq <- list(
    type = "lmrq",
    para = c("mu" = mu, "alpha" = alpha),
    source = "parlmqr"
  )
  return (lmomco::dlmomco(x, para = obj.parlmrq))
}
plmrq <- function(q, mu, alpha) {
  obj.parlmrq <- list(
    type = "lmrq",
    para = c("mu" = mu, "alpha" = alpha),
    source = "parlmqr"
  )
  return (lmomco::plmomco(q, para = obj.parlmrq))
}
qlmrq <- function(f, mu, alpha) {
  obj.parlmrq <- list(
    type = "lmrq",
    para = c("mu" = mu, "alpha" = alpha),
    source = "parlmqr"
  )
  return (lmomco::qlmomco(f, para = obj.parlmrq))
}
# ---------------------------------------------------------------------------- #
# Generalized Exponential Poisson 
rln3 <- function(n, zeta = 0, mulog, sigmalog) {
  obj.parln3 <- list(
    type = "ln3",
    para = c('zeta' = zeta, "mulog"= mulog, "sigmalog" = sigmalog),
    source = "parln3"
  )
  return (lmomco::rlmomco(n, para = obj.parln3))
}
dln3 <- function(x, zeta = 0, mulog, sigmalog) {
  obj.parln3 <- list(
    type = "ln3",
    para = c("zeta" = zeta, "mulog" = mulog, "sigmalog" = sigmalog),
    source = "parln3"
  )
  return (lmomco::dlmomco(x, para = obj.parln3))
}
pln3 <- function(q, zeta = 0, mulog, sigmalog) {
  obj.parln3 <- list(
    type = "ln3",
    para = c('zeta' = zeta, "mulog"= mulog, "sigmalog" = sigmalog),
    source = "parln3"
  )
  return (lmomco::plmomco(q, para = obj.parln3))
}
qln3 <- function(f, zeta = 0, mulog, sigmalog) {
  obj.parln3 <- list(
    type = "ln3",
    para = c('zeta' = zeta, "mulog"= mulog, "sigmalog" = sigmalog),
    source = "parln3"
  )
  return (lmomco::qlmomco(f, para = obj.parln3))
}
# ---------------------------------------------------------------------------- #
# Pearson III
rpe3 <- function(n, mu, sigma, gamma) {
  obj.parpe3 <- list(
    type = "pe3",
    para = c('mu' = mu, "sigma" = sigma, "gamma"= gamma),
    source = "parpe3"
  )
  return (lmomco::rlmomco(n, para = obj.parpe3))
}
dpe3 <- function(x, mu, sigma, gamma) {
  obj.parpe3 <- list(
    type = "pe3",
    para = c('mu' = mu, "sigma" = sigma, "gamma"= gamma),
    source = "parpe3"
  )
  return (lmomco::dlmomco(x, para = obj.parpe3))
}
ppe3 <- function(q, mu, sigma, gamma) {
  obj.parpe3 <- list(
    type = "pe3",
    para = c('mu' = mu, "sigma" = sigma, "gamma"= gamma),
    source = "parpe3"
  )
  return (lmomco::plmomco(q, para = obj.parpe3))
}
qpe3 <- function(f, mu, sigma, gamma) {
  obj.parpe3 <- list(
    type = "pe3",
    para = c('mu' = mu, "sigma" = sigma, "gamma"= gamma),
    source = "parpe3"
  )
  return (lmomco::qlmomco(f, para = obj.parpe3))
}
# ---------------------------------------------------------------------------- #
# Rayleigh
rray <- function(n, xi, alpha) {
  obj.parray <- list(
    type = "ray",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parray"
  )
  return (lmomco::rlmomco(n, para = obj.parray))
}
dray <- function(x, xi, alpha) {
  obj.parray <- list(
    type = "ray",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parray"
  )
  return (lmomco::dlmomco(x, para = obj.parray))
}
pray <- function(q, xi, alpha) {
  obj.parray <- list(
    type = "ray",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parray"
  )
  return (lmomco::plmomco(q, para = obj.parray))
}
qray <- function(f, xi, alpha) {
  obj.parray <- list(
    type = "ray",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parray"
  )
  return (lmomco::qlmomco(f, para = obj.parray))
}
# ---------------------------------------------------------------------------- #
# Rice
rrice <- function(n, nu, alpha) {
  obj.parrice <- list(
    type = "rice",
    para = c("nu" = nu, "alpha" = alpha),
    source = "parrice"
  )
  return (lmomco::rlmomco(n, para = obj.parrice))
}
drice <- function(x, nu, alpha) {
  obj.parrice <- list(
    type = "rice",
    para = c("nu" = nu, "alpha" = alpha),
    source = "parrice"
  )
  return (lmomco::dlmomco(x, para = obj.parrice))
}
price <- function(q, nu, alpha) {
  obj.parrice <- list(
    type = "rice",
    para = c("nu" = nu, "alpha" = alpha),
    source = "parrice"
  )
  return (lmomco::plmomco(q, para = obj.parrice))
}
qrice <- function(f, nu, alpha) {
  obj.parrice <- list(
    type = "rice",
    para = c("nu" = nu, "alpha" = alpha),
    source = "parrice"
  )
  return (lmomco::qlmomco(f, para = obj.parrice))
}
# ---------------------------------------------------------------------------- #
# Slash
rsla <- function(n, xi, alpha) {
  obj.parsla <- list(
    type = "sla",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parsla"
  )
  return (lmomco::rlmomco(n, para = obj.parsla))
}
dsla <- function(x, xi, alpha) {
  obj.parsla <- list(
    type = "sla",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parsla"
  )
  return (lmomco::dlmomco(x, para = obj.parsla))
}
psla <- function(q, xi, alpha) {
  obj.parsla <- list(
    type = "sla",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parsla"
  )
  return (lmomco::plmomco(q, para = obj.parsla))
}
qsla <- function(f, xi, alpha) {
  obj.parsla <- list(
    type = "sla",
    para = c("xi" = xi, "alpha" = alpha),
    source = "parsla"
  )
  return (lmomco::qlmomco(f, para = obj.parsla))
}
# ---------------------------------------------------------------------------- #
# 3-Parameter Student t
rst3 <- function(n, xi, alpha, nu) {
  obj.parst3 <- list(
    type = "st3",
    para = c('xi' = xi, "alpha" = alpha, "nu"= nu),
    source = "parst3"
  )
  return (lmomco::rlmomco(n, para = obj.parst3))
}
dst3 <- function(x, xi, alpha, nu) {
  obj.parst3 <- list(
    type = "st3",
    para = c('xi' = xi, "alpha" = alpha, "nu"= nu),
    source = "parst3"
  )
  return (lmomco::dlmomco(x, para = obj.parst3))
}
pst3 <- function(q, xi, alpha, nu) {
  obj.parst3 <- list(
    type = "st3",
    para = c('xi' = xi, "alpha" = alpha, "nu"= nu),
    source = "parst3"
  )
  return (lmomco::plmomco(q, para = obj.parst3))
}
qst3 <- function(f, xi, alpha, nu) {
  obj.parst3 <- list(
    type = "st3",
    para = c('xi' = xi, "alpha" = alpha, "nu"= nu),
    source = "parst3"
  )
  return (lmomco::qlmomco(f, para = obj.parst3))
}
# ---------------------------------------------------------------------------- #
# 3-Parameter Student t
rtri <- function(n, min, mode, max) {
  obj.partri <- list(
    type = "tri",
    para = c('min' = min, "mode" = mode, "max"= max),
    source = "partri"
  )
  return (lmomco::rlmomco(n, para = obj.parst3))
}
dtri <- function(x, min, mode, max) {
  obj.partri <- list(
    type = "tri",
    para = c('min' = min, "mode" = mode, "max"= max),
    source = "partri"
  )
  return (lmomco::dlmomco(x, para = obj.partri))
}
ptri <- function(q, min, mode, max) {
  obj.partri <- list(
    type = "tri",
    para = c('min' = min, "mode" = mode, "max"= max),
    source = "partri"
  )
  return (lmomco::plmomco(q, para = obj.partri))
}
qtri <- function(f, min, mode, max) {
  obj.partri <- list(
    type = "tri",
    para = c('min' = min, "mode" = mode, "max"= max),
    source = "partri"
  )
  return (lmomco::qlmomco(f, para = obj.partri))
}
# ---------------------------------------------------------------------------- #
# Wakeby
rwak <- function(n, xi, alpha, beta, gamma, delta) {
  obj.parwak <- list(
    type = "wak",
    para = c('xi' = xi, "alpha" = alpha, "beta"= beta, "gamma" = gamma, "delta" = delta),
    source = "parwak"
  )
  return (lmomco::rlmomco(n, para = obj.parwak))
}
dwak <- function(x, xi, alpha, beta, gamma, delta) {
  obj.parwak <- list(
    type = "wak",
    para = c('xi' = xi, "alpha" = alpha, "beta"= beta, "gamma" = gamma, "delta" = delta),
    source = "parwak"
  )
  return (lmomco::dlmomco(x, para = obj.parwak))
}
pwak <- function(q, xi, alpha, beta, gamma, delta) {
  obj.parwak <- list(
    type = "wak",
    para = c('xi' = xi, "alpha" = alpha, "beta"= beta, "gamma" = gamma, "delta" = delta),
    source = "parwak"
  )
  return (lmomco::plmomco(q, para = obj.parwak))
}
qwak <- function(f, xi, alpha, beta, gamma, delta) {
  obj.parwak <- list(
    type = "wak",
    para = c('xi' = xi, "alpha" = alpha, "beta"= beta, "gamma" = gamma, "delta" = delta),
    source = "parwak"
  )
  return (lmomco::plmomco(f, para = obj.parwak))
}
# ---------------------------------------------------------------------------- #
# Weibull
rwei <- function(n, zeta, beta, delta) {
  obj.parwei <- list(
    type = "wei",
    para = c('zeta' = zeta, "beta" = beta, "delta"= delta),
    source = "parwei"
  )
  return (lmomco::rlmomco(n, para = obj.parwei))
}
dwei <- function(x, zeta, beta, delta) {
  obj.parwei <- list(
    type = "wei",
    para = c('zeta' = zeta, "beta" = beta, "delta"= delta),
    source = "parwei"
  )
  return (lmomco::dlmomco(x, para = obj.parwei))
}
pwei <- function(q, zeta, beta, delta) {
  obj.parwei <- list(
    type = "wei",
    para = c('zeta' = zeta, "beta" = beta, "delta"= delta),
    source = "parwei"
  )
  return (lmomco::plmomco(q, para = obj.parwei))
}
qwei <- function(f, zeta, beta, delta) {
  obj.parwei <- list(
    type = "wei",
    para = c('zeta' = zeta, "beta" = beta, "delta"= delta),
    source = "parwei"
  )
  return (lmomco::qlmomco(f, para = obj.parwei))
}