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
  distribucion.estadistico <- copula::indepTestSim(n = length(x), p = 2, verbose = FALSE, N = 100)
  
  resultado.test <- copula::indepTest(matrix(c(x, y), ncol = 2, byrow = F), d = distribucion.estadistico)
  
  # Devolver objeto con los resultados
  resultado$p.value <- unname(resultado.test$pvalues)
  resultado$statistic <- unname(resultado.test$statistics)
  
  return(resultado)
}

# Definicion de funcion para consolidar tests de independencia
SonDependientes <- function(x, y, tests.dependencia, umbral.p.valor) {
  
  resultados.tests <- purrr::accumulate(
    .x = tests.dependencia,
    .f = function(prev.result, test.name) {
      if (!purrr::is_empty(prev.result$pasa.test))
        if (!dplyr::last(prev.result$pasa.test))
          return (prev.result)
      func.name <- paste0("TestInd", test.name)
      tryCatch({
        p.value <- ParametrosADataFrame(
          do.call(what = func.name, args = list(x = x, y = y))) %>%
          dplyr::filter(parametro == "p.value") %>%
          dplyr::pull(valor)
        if (is.na(p.value) || p.value > umbral.p.valor) 
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

