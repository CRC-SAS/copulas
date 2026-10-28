# Periodo de retorno combinado (duracion, intensidad) a partir de la copula duracion-magnitud, y su validacion.
# Como magnitud = duracion * intensidad (exacto), P(D > d, I > i) = P(D > d, M > i * D). Se calcula simulando de la copula + marginales
# ajustadas de (duracion, magnitud) y contando la fraccion de pares que cumple la condicion.
# Validacion: (1) donde existen ambas copulas (intensidad-duracion y duracion-magnitud) se comparan los dos calculos y las frecuencias observadas;
#             (2) en todos los pares con copula duracion-magnitud se compara el numero esperado de eventos (n * P) con el observado,
#                 contra el supuesto de independencia con marginales empiricas.
# Se corre desde la raiz del repo. Entradas: data/output/copulas_*.rds, data/input/eventos_identificados_frio_calor.csv
# Salidas: DIR_PRES/via_magnitud_validacion.csv y DIR_PRES/via_magnitud.txt
DIR_PRES <- Sys.getenv("DIR_PRES", "data/output/presentacion")
NSIM <- as.integer(Sys.getenv("NSIM", "400000"))
suppressMessages({library(dplyr); library(copula)})
source("lib/funciones_ajuste_distribuciones.R", encoding = "UTF-8")
EXCLUIDA <- "87420"

# Probabilidad conjunta P(D > d, I > i) para (d, i) vectores, simulando del Mvdc ajustado (x = duracion, y = magnitud)
ProbViaMagnitud <- function(mv, d, i, nsim = NSIM, semilla = 20261023) {
  set.seed(semilla)
  u <- copula::rCopula(nsim, mv@copula)
  qx <- get(paste0("q", mv@margins[1])); qy <- get(paste0("q", mv@margins[2]))
  D <- do.call(qx, c(list(u[, 1]), mv@paramMargins[[1]])); M <- do.call(qy, c(list(u[, 2]), mv@paramMargins[[2]]))
  ok <- is.finite(D) & is.finite(M); D <- D[ok]; M <- M[ok]; I <- M / D
  mapply(function(dd, ii) mean(D > dd & I > ii), d, i)
}
# Probabilidad conjunta con la copula intensidad-duracion (la que usa el pipeline): 1 - u - v + C(u, v)
ProbDirecta <- function(mv, d, i) {
  Fi <- do.call(get(paste0("p", mv@margins[1])), c(list(q = i), mv@paramMargins[[1]]))
  Fd <- do.call(get(paste0("p", mv@margins[2])), c(list(q = d), mv@paramMargins[[2]]))
  pmax(1 - Fi - Fd + copula::pCopula(cbind(Fi, Fd), mv@copula), 0)
}

if (sys.nframe() == 0) {
  sink(file.path(DIR_PRES, "via_magnitud.txt"), split = TRUE)
  slugs <- c(olacalor = "OlaCalorTXyTN", olafrio = "OlaFrioTXyTN", periodocalidotn = "PeriodoCalidoTN", periodocalidotx = "PeriodoCalidoTX",
             periodofriotn = "PeriodoFrioTN", periodofriotx = "PeriodoFrioTX")
  ev <- readr::read_csv("data/input/eventos_identificados_frio_calor.csv", show_col_types = FALSE) %>% filter(duracion >= 3) %>% mutate(station_id = as.character(station_id))
  res <- list()
  for (s in names(slugs)) {
    x <- readRDS(sprintf("data/output/copulas_%s.rds", s)) %>% filter(tipo_serie == "observada", as.character(station_id) != EXCLUIDA) %>% mutate(station_id = as.character(station_id))
    x$ok <- vapply(x$mejor_copula, isS4, logical(1)) & vapply(x$ajuste_multivariado, isS4, logical(1))  # un caso (87467 olas de calor) tiene copula pero no el Mvdc guardado
    dm <- x %>% filter(ok, variable_x == "duracion", variable_y == "magnitud"); di <- x %>% filter(ok, variable_x == "intensidad", variable_y == "duracion")
    for (k in seq_len(nrow(dm))) {
      st <- dm$station_id[k]; e <- ev %>% filter(station_id == st, tipo_evento == slugs[[s]]); n <- nrow(e)
      if (n < 20) next
      # umbrales: duraciones enteras 3..8 (las que tienen eventos) x cuantiles de intensidad
      # d entero (D > d, el caso de uso del pipeline, que trata la duracion entera como continua) y d + 0,5 (con correccion de continuidad:
      # D > k + 0,5 equivale a D >= k + 1 en datos enteros)
      ds <- sort(unique(pmin(e$duracion, 8))); ds <- ds[ds <= quantile(e$duracion, .95)]
      qi <- quantile(e$intensidad, c(.25, .5, .75, .9))
      g <- bind_rows(expand.grid(d = ds, i = as.numeric(qi)) %>% mutate(corr = "d entero"), expand.grid(d = ds + .5, i = as.numeric(qi)) %>% mutate(corr = "d + 0,5"))
      g$obs <- mapply(function(dd, ii) sum(e$duracion > dd & e$intensidad > ii), g$d, g$i)
      g$p_dm <- ProbViaMagnitud(dm$ajuste_multivariado[[k]], g$d, g$i)
      g$p_ind_emp <- mapply(function(dd, ii) mean(e$duracion > dd) * mean(e$intensidad > ii), g$d, g$i)
      g$p_emp <- g$obs / n
      j <- which(di$station_id == st)
      g$p_di <- if (length(j) == 1) ProbDirecta(di$ajuste_multivariado[[j]], g$d, g$i) else NA_real_
      g$slug <- s; g$station_id <- st; g$n <- n; g$ambas <- length(j) == 1
      res[[length(res) + 1]] <- g
    }
    cat(s, "listo:", nrow(dm), "pares con copula duracion-magnitud\n")
  }
  r <- bind_rows(res)
  readr::write_csv(r, file.path(DIR_PRES, "via_magnitud_validacion.csv"))
  sink()
}
