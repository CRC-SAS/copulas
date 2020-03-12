# -----------------------------------------------------------------------------#
# ---- Funciones de bondad de ajuste ----

# Definición de función para el test Sn
TestSn <- function(copula = NA, x, y) {
  
  # Test de Cramer-von Mises S[n] definido en la ecuación (2) en Genest,
  # Remillard and Beaudoin (2009).
  
  # Inicializar objeto resultado
  resultado <- list(p.value = NA, estimate = NA)
  
  # Creacion de matriz de eventos en valores absolutos para evitar dependencias negativas
  matriz.eventos <- as.matrix(cbind(abs(x), abs(y)))
  # Convertir la matriz de eventos en pseudoobservaciones
  matriz.eventos.pobs <- copula::pobs(matriz.eventos, ties.method = 'random')
  
  # La matriz debe tener al menos dos variables
  base::stopifnot(ncol(matriz.eventos.pobs) >= 2)
  
  # Prueba de bondad de ajuste
  resultados.test <-   copula::gofCopula(copula = copula$copula, x = matriz.eventos.pobs,
    N = 1000, simulation = "pb", estim.method = 'mpl', optim.method = 'BFGS',
    method = 'Sn')
  
  # Devolver objeto resultado
  resultado$p.value <- unname(resultados.test$p.value)
  resultado$estimate <- unname(resultados.test$statistic)
  
  return(resultado)
}

# Definición de función para el test SnB
TestSnB <- function(copula = NA, x, y) {
  # Inicializar objeto resultado
  resultado <- list(p.value = NA, estimate = NA)
  
  # Creacion de matriz de eventos en valores absolutos para evitar dependencias negativas
  matriz.eventos <- as.matrix(cbind(abs(x), abs(y)))
  # Convertir la matriz de eventos en pseudoobservaciones
  matriz.eventos.pobs <- copula::pobs(matriz.eventos, ties.method = 'random')
  
  # La matriz debe tener al menos dos variables
  base::stopifnot(ncol(matriz.eventos.pobs) >= 2)
  
  # Prueba de bondad de ajuste
  resultados.test <-   copula::gofCopula(copula = copula$copula, x = matriz.eventos.pobs,
    N = 1000, simulation = "pb", estim.method = 'mpl', optim.method = 'BFGS',
    method = 'SnB')
  
  # Devolver objeto resultado
  resultado$p.value <- unname(resultados.test$p.value)
  resultado$estimate <- unname(resultados.test$statistic)
  
  return(resultado)
}

# Definición de función para el test Anderson-Darling
TestSnC <- function(copula = NA, x, y) {
  # Inicializar objeto resultado
  resultado <- list(p.value = NA, estimate = NA)
  
  # Creacion de matriz de eventos en valores absolutos para evitar dependencias negativas
  matriz.eventos <- as.matrix(cbind(abs(x), abs(y)))
  # Convertir la matriz de eventos en pseudoobservaciones
  matriz.eventos.pobs <- copula::pobs(matriz.eventos, ties.method = 'random')
  
  # La matriz debe tener al menos dos variables
  base::stopifnot(ncol(matriz.eventos.pobs) >= 2)
  
  # Prueba de bondad de ajuste
  resultados.test <-   copula::gofCopula(copula = copula$copula, x = matriz.eventos.pobs,
    N = 1000, simulation = "pb", estim.method = 'mpl', optim.method = 'BFGS',
    method = 'SnC')
  
  # Devolver objeto resultado
  resultado$p.value <- unname(resultados.test$p.value)
  resultado$estimate <- unname(resultados.test$statistic)
  
  return(resultado)
  
}

# Definicion de funcion para el AIC
TestAIC <- function(copula = NA) {
  
  # Inicializar objeto resultado
  resultado <- list(AIC = NA)
  
  # Calcular el crterio de Akaike
  AIC <- -2 * copula$bondad.ajuste$loglike + 2 * length(copula$parametro)
  
  # Devolver objeto resultado
  resultado$AIC <- AIC
  
  return(AIC)
  
}

# Definicion de funcion para el BIC
TestBIC <- function(copula = NA, x, y) {
  
  # Inicializar objeto resultado
  resultado <- list(BIC = NA)
  
  # Creacion de matriz de eventos en valores absolutos para evitar dependencias negativas
  matriz.eventos <- as.matrix(cbind(abs(x), abs(y)))
  
  # La matriz debe tener al menos dos variables
  base::stopifnot(ncol(matriz.eventos) >= 2)
  
  # Calcular el crterio de Akaike
  BIC <- -2 * copula$bondad.ajuste$loglike + log(ncol(matriz.eventos)) * length(copula$parametro)
  
  # Devolver objeto resultado
  resultado$BIC <- BIC
  
  return(BIC)
}

# Definicion de funcion para el RMSE
TestError  <- function(copula = NA, n = 1000) {
  
  # Inicializar objeto resultado
  resultado <- list(mre = NA, rmse = NA)
  
  # Estimación de la copula empírica
  d <- dim(copula$copula)
  U <- copula::rCopula(n, copula = copula$copula) # Muestra de puntos de la copula teórica
  v <- matrix(runif(n * d), nrow = n, ncol = d) # Matriz de puntos a evaluar
  copula.empirica    <- copula::C.n(v, X = U) # Estimación de la copula empitica
  
  # Estimaci/on de la copula teorica
  copula.teorica <- copula::pCopula(u = v, copula = copula$copula) 
  
  # Calculo del error relativo medio
  MRE <-  mean(abs(copula.teorica - copula.empirica) / copula.teorica) * 100
  
  # Calculo del RMSE
  RMSE <- hydroGOF::rmse(sim = copula.teorica, obs = copula.empirica) * 100
  
  # Devolver resultados
  resultado$mre <- MRE
  resultado$rmse <- RMSE
  
  return(resultado)
  
}

# Definicion de funcion para la validacion cruzada
TestValidacionCruzada <- function(copula, x, y) {
  
  # Inicializar objetos resultados
  resultado <- list(CIC = NA)
  
  # Creacion de matriz de eventos en valores absolutos para evitar dependencias negativas
  matriz.eventos <- as.matrix(cbind(abs(x), abs(y)))
  
  # La matriz debe tener al menos dos variables
  base::stopifnot(ncol(matriz.eventos) >= 2)
  
  # Creacion de objeto con los parametros del test
  parametros.test <- list(x = matriz.eventos, copula = copula$copula)
  
  # Realizar validacion curzada
  resultado.test <- do.call(what = copula::xvCopula, args = parametros.test)
  
  # Devolver resultados
  resultado$CIC <- resultado.test
  
  return(resultado)
  
}


# --- Consolidacion de tests
TestearBondadAjusteCopulas <- function(x, y,  umbral.p.valor, copula = NULL) {
  
  
  # 2. Inicializar objeto a devolver
  resultados.tests <- list(pasa.tests = TRUE, estadisticos = NULL)
  
  # 3. Aplicar tests exclusivos para caso parametrico
  #    Si alguno de los tests devuelve NA o un valor de p-value menor al umbral,
  #    interpretar el resultado del test como un fallo. Luego, si alguno de los tests falla, entonces
  #    interpretar como malo el ajuste y devolver todos los parametros en NA.
  if (! is.null(copula$copula)) {
    falla.ajuste <- any(is.na(copula$parametro))
    if (! falla.ajuste) {
      estadisticos <- purrr::map_dfr(
        .x = c("Sn"),
        .f = function(test.name) {
          func.name <- paste0("Test", test.name)
          tryCatch({
            estadisticos.test <- ParametrosADataFrame(do.call(what = func.name, args = list(x = x, y = y,  copula = copula))) %>%
              dplyr::mutate(test = test.name) %>%
              dplyr::select(test, parametro, valor)
          }, error = function(e) {
            cat(e$message, "\n")
            return (NULL)
          })
        }
      )
      resultados.tests$estadisticos <- estadisticos
      
      # Determinar si pasan los tests o no
      p.values <- estadisticos %>%
        dplyr::filter(parametro == "p.value") %>%
        dplyr::pull(valor)
      if (any(is.na(p.values)) || any(p.values < umbral.p.valor)) {
        resultados.tests$pasa.tests <- FALSE  
      } else {
        resultados.tests$pasa.tests <- TRUE 
      }
    }
  }
  
  # 4. Aplicar tests de AIC, BIC, RME, RMSE y Validacion cruzada 
  #    Guardar los valores resultantes, pero no dictaminar en base a esos tests.
  if (! falla.ajuste) {
    # i. AIC
    aic <- TestAIC(copula = copula)
    resultados.tests$estadisticos <- rbind(resultados.tests$estadisticos, data.frame(test = 'AIC', parametro = 'aic', valor = aic))
    
    # ii. BIC
    bic <- TestBIC(x = x, y = y, copula = copula)
    resultados.tests$estadisticos <- rbind(resultados.tests$estadisticos, data.frame(test = 'BIC', parametro = 'bic', valor = bic))
    
    # iii. Error
    error <- TestError(copula = copula, n = 1000)
    resultados.tests$estadisticos <- rbind(resultados.tests$estadisticos, 
      data.frame(test = c('MRE', 'RMSE'), 
        parametro = c('mre', 'rmse'), 
        valor = c(error[[1]], error[[2]])))
    
    # iv. Validacion cruzada
    validacion.cruzada <- TestValidacionCruzada(copula = copula, x = x, y = y)
    resultados.tests$estadisticos <- rbind(resultados.tests$estadisticos, data.frame(test = 'Validacion cruzada', parametro = 'cic', valor = unlist(unname(validacion.cruzada))))
  }
  
  # 5. Devolver resultados
  return (resultados.tests)
}

