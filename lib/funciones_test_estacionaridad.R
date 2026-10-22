DetectarEmpates <- function(x) {
  # Contar la cantidad de e,pates
  cantidad.empates <- length(x) - length(unique(x))
  
  # Devolver restulado
  return(cantidad.empates)
}


# Definici?n de funci?n para el test de estacionaridad de Box-Pierce y Ljung-Box (Univariado)
TestBoxPierceLjungBox <- function(x, fechas) {
  
  # Inicializar objeto a devolver
  resultado <- list(p.value = NA, statistic = NA)
  # Correr test de indepedencia (No correlacion)
  resultado.test <- Box.test(x, lag =  3, type = "Ljung-Box")
  
  # Devolver objeto con los resultados
  resultado$p.value <- unname(resultado.test$p.value)
  resultado$statistic <- unname(resultado.test$statistic)
  
  return(resultado)
}

# Definici?n de funci?n para el test de estacionaridad emp?rico
TestEstacionaridadEmpirico <- function(x, fechas) {
  
  # Inicializar objeto a devolver
  resultado <- list(p.value = NA, statistic = NA)
  # Comprobar que el vector sea del tipo xts
  if(!xts::is.xts(x)) {
    x.xts <- xts::xts(x = x, order.by = fechas)
  } else {
    x.xts <- x
  }
  
  # Crear distribuci;on del estad?stico
  sI.distribucion <- copula::serialIndepTestSim(nrow(x.xts), lag.max = 10, N = 100) 
  
  # Correr test
  resultado.test <- copula::serialIndepTest(x.xts, d = sI.distribucion)
  
  # Devolver objeto con los resultados
  resultado$p.value <- unname(resultado.test$global.statistic.pvalue)
  resultado$statistic <- unname(resultado.test$global.statistic)
  
  # Devolver resultados
  return(resultado)
}

# Definicion de funcion para el test de punto de cambio
TestPuntoCambioUni <- function(x, fechas) {
  
  # Inicializar objeto a devolver
  resultado <- list(p.value = NA, statistic = NA)
  # Comprobar que el vector sea del tipo xts
  if(!xts::is.xts(x)) {
    x.xts <- xts::xts(x = x, order.by = fechas)
  } else {
    x.xts <- x
  }
  
  # Correr test
  resultado.test <- npcp::cpDist(x.xts, b = NULL, method = 'seq')
  
  # Devolver objeto con los resultados
  resultado$p.value <- unname(resultado.test$p.value)
  resultado$statistic <- unname(resultado.test$statistic)
  
  # Devolver resultados
  return(resultado)
  
}

# Definicion de funcion para el test de punto de cambio en la autocopula
TestPuntoCambioAutocopula <- function(x, fechas) {
  
  # Inicializar objeto a devolver
  resultado <- list(p.value = NA, statistic = NA)
  # Comprobar que el vector sea del tipo xts
  if(!xts::is.xts(x)) {
    x.xts <- xts::xts(x = x, order.by = fechas)
  } else {
    x.xts <- x
  }
  
  # Correr test
  resultado.test <- npcp::cpAutocop(x.xts, lag = 1)
  
  # Devolver objeto con los resultados
  resultado$p.value <- unname(resultado.test$p.value)
  resultado$statistic <- unname(resultado.test$statistic)
  
  # Devolver resultados
  return(resultado)
  
}

# Definicion de funcion para el test de Mann-Kendall
TestMannKendall <- function(x, fechas) {
  
  # Inicializar objeto a devolver
  resultado <- list(p.value = NA, statistic = NA)
  # Comprobar que el vector sea del tipo xts
  if(!xts::is.xts(x)) {
    x.xts <- xts::xts(x = x, order.by = fechas)
  } else {
    x.xts <- x
  }
  
  # Correr test
  resultado.test <- Kendall::MannKendall(x.xts)
  
  # Devolver objeto con los resultados
  resultado$p.value <- unname(resultado.test$sl)
  resultado$statistic <- unname(resultado.test$tau)
  
  # Devolver resultados
  return(resultado)
  
}

# Definicion de funcion para consolidar tests de estacionariedad 
EsEstacionaria <- function(x, fechas, tests.estacionariedad, umbral.p.valor) {
  
  resultados.tests <- purrr::reduce(
    .x = tests.estacionariedad,
    .f = function(prev.result, test.name) {
      if (!purrr::is_empty(prev.result$pasa.test))
        if (!dplyr::last(prev.result$pasa.test)) 
          return (prev.result)
      func.name <- paste0("Test", test.name)
      tryCatch({
        p.value <- ParametrosADataFrame(
          do.call(what = func.name, args = list(x = x, fechas = fechas))) %>%
          dplyr::filter(parametro == "p.value") %>%
          dplyr::pull(valor)
        # En todos los tests de estacionariedad la hipotesis nula es que la serie es estacionaria
        # (sin autocorrelacion, sin punto de cambio, sin tendencia): el test se pasa si NO se rechaza H0
        if (is.na(p.value) || p.value <= umbral.p.valor) 
          return (prev.result %>% dplyr::bind_rows(
            tibble::tibble(pasa.test = FALSE, test.name = test.name, p.value = p.value, exec.error = FALSE)))
        else
          return (prev.result %>% dplyr::bind_rows(
            tibble::tibble(pasa.test = TRUE, test.name = test.name, p.value = p.value, exec.error = FALSE)))
      }, error = function(e) {
        cat(e$message, "\n")
        return (prev.result %>% dplyr::bind_rows(
          tibble::tibble(pasa.test = FALSE, test.name = test.name, p.value = NA, exec.error = TRUE)))
      })
    },
    .init = tibble::tibble(pasa.test = logical(), test.name = character(), p.value = double(), exec.error = logical())
  )
  
  # Devolver resultados
  return (resultados.tests)
}

# Definicion de funcion para el test de punto de cambio multivariado
TestPuntoCambioMulti <- function(x, variable = c('duracion.prima', 'intensidad')) {
  
  # Inicializar objeto a devolver
  resultado <- list(p.value = NA, statistic = NA)
  # Comprobar que el vector sea del tipo xts
  if(!xts::is.xts(x)) {
    x.xts <- xts::xts(x = x[,variable], order.by = x$fecha_inicio)
  } else {
    x.xts <- x
  }
  
  # Correr test
  resultado.test <- npcp::cpCopula(x.xts, b = NULL, method = "seq")
  
  # Devolver objeto con los resultados
  resultado$p.value <- unname(resultado.test$p.value)
  resultado$statistic <- unname(resultado.test$statistic)
  
  # Devolver resultados
  return(resultado)
  
}




# Umbral de p-valor a aplicar a CADA test de estacionariedad. Como la serie se considera estacionaria
# solo si pasa los n.tests tests, con umbral individual alfa el error global (rechazar una serie
# estacionaria) crece hasta 1-(1-alfa)^n.tests. La correccion de Bonferroni usa alfa/n.tests por test
# para mantener el error global en alfa.
#   correccion = "bonferroni": umbral / n.tests (valor por defecto si no se especifica)
#   correccion = "ninguna":    umbral sin corregir
UmbralEstacionariedad <- function(umbral, n.tests, correccion = "bonferroni") {
  if (is.null(correccion)) correccion <- "bonferroni"
  if (!correccion %in% c("bonferroni", "ninguna"))
    stop("correccion.estacionariedad debe ser 'bonferroni' o 'ninguna'")
  if (correccion == "bonferroni") return (umbral / max(n.tests, 1))
  return (umbral)
}

# Consolida las realizaciones de una misma serie: estacionaria si lo es la mayoria (empate = estacionaria).
# Atencion: con round(n/2) una sola realizacion (serie observada) siempre resultaba estacionaria, porque round(0.5) = 0
ConsolidarEstacionariedad <- function(es.estacionaria) {
  mean(es.estacionaria) >= 0.5
}
