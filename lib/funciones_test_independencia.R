# --- Funciones para testear la independencia de variables ----

# Definici?n de funci?n para el test de independencia de Kendall (Bivariado)
TestIndKendall <- function(x, y) {
  
  # Inicializar objeto a devolver
  resultado <- list(p.value = NA, statistic = NA)
  # Correr test de indepedencia (No correlacion)
  resultado.test <- cor.test(x, y, method = 'kendall')
  
  # Devolver objeto con los resultados
  resultado$p.value <- unname(resultado.test$p.value)
  resultado$statistic <- unname(resultado.test$estimate)
  
  return(resultado)
}

# Definici?n de funci?n para el test de independencia de Spearman (Bivariado)
TestIndSpearman <- function(x, y) {
  
  # Inicializar objeto a devolver
  resultado <- list(p.value = NA, statistic = NA)
  # Correr test de indepedencia (No correlacion)
  resultado.test <- cor.test(x, y, method = 'spearman')
  
  # Devolver objeto con los resultados
  resultado$p.value <- unname(resultado.test$p.value)
  resultado$statistic <- unname(resultado.test$estimate)
  
  return(resultado)
}

# Definici?n de funci?n para el test de independencia basado en la copula empirica (Bivariado)
TestIndCopulaEmpirica <- function(x, y) {
  # Inicializar objeto a devolver
  resultado <- list(p.value = NA, statistic = NA)
  
  # Correr test de independencia
  # Creaci?n de la distribuci?n del estad?stico
  distribucion.estadistico <- copula::indepTestSim(n = length(x), 
                                           p = 2, verbose = FALSE)
  
  resultado.test <- copula::indepTest(matrix(c(x, y), ncol = 2, byrow = F), 
                              d = distribucion.estadistico)
  
  # Devolver objeto con los resultados
  resultado$p.value <- unname(resultado.test$pvalues)
  resultado$statistic <- unname(resultado.test$statistics)
  
  return(resultado)
}

# Definicion de funcion para consolidar tests de independencia
SonDependientes <- function(x, y, umbral.p.valor) {
  # Inicializar objeto a devolver
  resultados.tests <- list(pasa.tests = TRUE)
  
  estadisticos <- purrr::map_dfr(
    .x = c("Kendall", "Spearman", "CopulaEmpirica"),
    .f = function(test.name) {
      func.name <- paste0("TestInd", test.name)
      tryCatch({
        estadisticos.test <- ParametrosADataFrame(do.call(what = func.name, args = list(x = x, y = y))) %>%
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
  if (any(is.na(p.values)) || any(p.values > umbral.p.valor)) {
    resultados.tests$pasa.tests <- FALSE  
  } else {
    resultados.tests$pasa.tests <- TRUE 
  }
  # Devolver resultados
  return (resultados.tests)
}

