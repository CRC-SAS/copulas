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
  sI.distribucion <- copula::serialIndepTestSim(nrow(x.xts), lag.max = 10) 
  
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

# Definicion de funcion para consolidar tests de estacionaridad 
EsEstacionaria <- function(x, fechas, umbral.p.valor) {
  # Inicializar objeto a devolver
  resultados.tests <- list(pasa.tests = TRUE)
  
  estadisticos <- purrr::map_dfr(
    .x = c("BoxPierceLjungBox", "EstacionaridadEmpirico", "PuntoCambioUni", 'PuntoCambioAutocopula', 'MannKendall'),
    .f = function(test.name) {
      func.name <- paste0("Test", test.name)
      tryCatch({
        estadisticos.test <- ParametrosADataFrame(do.call(what = func.name, args = list(x = x, fechas = fechas))) %>%
          dplyr::mutate(test = test.name) %>%
          dplyr::select(test, parametro, valor)
      }, error = function(e) {
        cat(e$message, "\n")
        return (NULL)
      })
    }
  )
  
  # Determinar si pasan los tests o no
  p.values <- estadisticos %>%
    dplyr::filter(parametro == "p.value") %>%
    dplyr::pull(valor)
  if (any(is.na(p.values)) || any(p.values < umbral.p.valor)) {
    resultados.tests$pasa.tests <- FALSE  
  } else {
    resultados.tests$pasa.tests <- TRUE 
  }
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



