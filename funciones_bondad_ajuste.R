# --- Funciones para estimar la bondad del ajuste

# Definición de función para el test de Kolmogorov-Smirnov
TestKS <- function(x, ajuste) {
  parametros.test <- list(x = x, y = paste0("p", ajuste$distribucion))
  for (nombre.parametro in names(ajuste$parametros)) {
    parametros.test[[nombre.parametro]] <- ajuste$parametros[[nombre.parametro]]
  }
  resultado.test <- do.call(what = stats::ks.test, args = parametros.test)
  resultados     <- list(p.value = unname(resultado.test$p.value),
    statistic = unname(resultado.test$statistic))
  
  return(resultados)
}

# Definición de función para el test de Anderson-Darling
TestAD <- function(x, ajuste) {
  parametros.test <- list(x = x, null = paste0("p", ajuste$distribucion))
  for (nombre.parametro in names(ajuste$parametros)) {
    parametros.test[[nombre.parametro]] <- ajuste$parametros[[nombre.parametro]]
  }
  resultado.test <- do.call(what = goftest::ad.test, args = parametros.test)
  resultados     <- list(p.value = unname(resultado.test$p.value),
    statistic = unname(resultado.test$statistic))
  
  return(resultados)
}

# Definición de función para el test de Cramer von Mises
TestCvM <- function(x, ajuste) {
  parametros.test <- list(x = x, null = paste0("p", ajuste$distribucion))
  for (nombre.parametro in names(ajuste$parametros)) {
    parametros.test[[nombre.parametro]] <- ajuste$parametros[[nombre.parametro]]
  }
  resultado.test <- do.call(what = goftest::cvm.test, args = parametros.test)
  resultados     <- list(p.value = unname(resultado.test$p.value),
    statistic = unname(resultado.test$statistic))
  
  return(resultados)
}

# Funciones para estimar bondad de ajuste en base a metricas de cuantiles
CalcularCuantilesAjustados <- function(p = NULL, ajuste = NULL) {
  parametros.cuantiles <- list(f = p)
  for (nombre.parametro in names(ajuste$parametros)) {
    parametros.cuantiles[[nombre.parametro]] <- ajuste$parametros[[nombre.parametro]]
  }
  cuantiles <- do.call(what = paste0("q", ajuste$distribucion), 
    args = parametros.cuantiles)
  
  return(cuantiles)
}

# Definición de función para el cálculo del RMSE estandarizado
TestRMSEIQR <- function(x, probs = seq(from = 0.01, to = 0.99, by = 0.1),
  ajuste = NULL) {
  # Calculo del RMSE entre los cuantiles observados y los estimados a partir de una distribucion
  # Validar parametros
  base::stopifnot(is.numeric(x))
  base::stopifnot(all(probs >= 0.0 && probs <= 1.0))
  base::stopifnot(
    (! is.null(ajuste) && (class(ajuste) == "list")))
  
  # Devolver NA si no hubo ajuste
  if (! is.null(ajuste) && any(is.na(ajuste))) {
    return (NA)
  }
  
  # Calcular los cuantiles empiricos
  cuant.observados <- quantile(x, probs = probs)
  
  # Calcular los cuantiles de la PDF ajustada
  cuant.ajustados <- CalcularCuantilesAjustados(probs, ajuste)
  
  # Calcular RMSE entre cuantiles empiricos y teoricos
  rmse <- caret::RMSE(pred = cuant.ajustados, obs = cuant.observados)
  
  # Calcular IQR de la muestra
  iqr <- stats::IQR(x)
  
  return (rmse/iqr)
}

# Definición de función para el test de concordancia
TestCCC <- function(x, probs = seq(from = 0.01, to = 0.99, by = 0.01), 
  ajuste = NULL) {
  # Coeficiente de correlacion de concordancia es una medida de correlacion, consistencia y precisión basada en Lin, L. (1989). 
  # A concordance correlation coefficient to evaluate reproducibility. Biometrics, 45 (1), 255???268.
  
  # Validar parametros
  base::stopifnot(is.numeric(x))
  base::stopifnot(all(probs >= 0.0 && probs <= 1.0))
  base::stopifnot(
    (! is.null(ajuste) && (class(ajuste) == "list"))
  ) 
  
  # Devolver NA si no hubo ajuste
  if (! is.null(ajuste) && any(is.na(ajuste))) {
    return (NA)
  }
  
  # Calcular los cuantiles empiricos
  cuant.observados <- quantile(x, probs = probs)
  
  # Calcular los cuantiles de la PDF ajustada
  cuant.ajustados <- CalcularCuantilesAjustados(probs, ajuste)
  
  # Crear objeto para el calculo del RMSE
  cuantiles <- data.frame(observados = cuant.observados, ajustados = cuant.ajustados)
  
  # Calcular RMSE entre cuantiles empiricos y teoricos
  ccc <- yardstick::ccc(cuantiles, truth = observados, estimate = ajustados)
  if (is.tbl(ccc)) {
    return (dplyr::pull(ccc, .estimate))
  } else {
    return (ccc)
  }
}

# Definición de función del test de comparación de cuantiles
TestQCOMHD <- function(x, probs = seq(from = 0.1, to = 1, by = 0.1), 
  numero.muestras, ajuste = NULL) {
  # Prueba de comparación de cuartiles basada en Rand R. Wilcox, David M. Erceg-Hurn, Florence Clark & Michael Carlson (2014) Comparing 
  # two independent groups via the lower and upper quantiles, Journal of Statistical Computation and Simulation, 84:7, 1543-1551, 
  # DOI: 10.1080/00949655.2012.754026
  
  # Validar parametros
  base::stopifnot(is.numeric(x))
  base::stopifnot(all(probs >= 0.0 && probs <= 1.0))
  base::stopifnot(is.integer(numero.muestras))
  base::stopifnot(
    (! is.null(ajuste) && (class(ajuste) == "list"))
  )  
  
  # Devolver NA si no hubo ajuste
  if (! is.null(ajuste) && any(is.na(ajuste))) {
    return (data.frame(q = probs, p.crit = rep(NA, length(probs)), p.value = rep(NA, length(probs))))
  }
  
  parametros.cuantiles <- list(n = length(x))
  for (nombre.parametro in names(ajuste$parametros)) {
    parametros.cuantiles[[nombre.parametro]] <- ajuste$parametros[[nombre.parametro]]
  }
  muestra <- do.call(what = paste0("r", ajuste$distribucion), 
    args = parametros.cuantiles)
  
  
  # Calcular t-student sobre cada percentil definido en probs entre la muestra 
  # original y una serie sintetica a partir del ajuste de una distribucion
  resultado <- WRS2::Dqcomhd(x = x, y = muestra, 
    q = probs, nboot = numero.muestras)
  qcomhd.resultado    <- resultado$partable %>%
    dplyr::select(q, p.crit, p.value)
  return(qcomhd.resultado)
}

ParametrosADataFrame <- function(parametros.ajuste) {
  parametros <- tibble::enframe(parametros.ajuste) %>%
    tidyr::unnest() %>%
    dplyr::rename(parametro = name, valor = value)
  return (parametros)
}

ParametrosALista <- function(parametros.ajuste) {
  parametros.lista        <- as.list(parametros.ajuste$valor)
  names(parametros.lista) <- parametros.ajuste$parametro
  return (parametros.lista)
}

# --- Consolidacion de tests
TestearBondadAjuste <- function(x, umbral.p.valor, ajuste = NULL) {
  
  
  # 2. Inicializar objeto a devolver
  resultados.tests <- list(pasa.tests = TRUE, estadisticos = NULL)
  
  # 3. Aplicar tests exclusivos para caso parametrico
  #    Si alguno de los tests devuelve NA o un valor de p-value menor al umbral,
  #    interpretar el resultado del test como un fallo. Luego, si alguno de los tests falla, entonces
  #    interpretar como malo el ajuste y devolver todos los parametros en NA.
  falla.ajuste.parametrico <- FALSE
  if (! is.null(ajuste)) {
    falla.ajuste.parametrico <- any(is.na(ajuste))
    if (! falla.ajuste.parametrico) {
      estadisticos <- purrr::map_dfr(
        .x = c("KS", "AD", "CvM"),
        .f = function(test.name) {
          func.name <- paste0("Test", test.name)
          tryCatch({
            estadisticos.test <- ParametrosADataFrame(do.call(what = func.name, args = list(x = x, ajuste = ajuste))) %>%
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
    } else {
      resultados.tests$pasa.tests <- FALSE
    }
  }
  
  # 4. Aplicar tests de RMST, CCC y QCOMHD (aplicables a ajustes parametricos como no parametricos)
  #    Guardar los valores resultantes, pero no dictaminar en base a esos tests.
  if (! falla.ajuste.parametrico) {
    # i. RMSE / IQR
    rmseiqr <- TestRMSEIQR(x = x, probs = seq(from=0.01, to=0.99, by=0.1),
      ajuste = ajuste)
    resultados.tests$estadisticos <- rbind(resultados.tests$estadisticos, data.frame(test = 'RMSEIQR', parametro = 'rmseiqr', valor = rmseiqr))
    
    # ii. CCC
    ccc <- TestCCC(x = x, probs = seq(from=0.01, to=0.99, by=0.01),
      ajuste = ajuste)
    resultados.tests$estadisticos <- rbind(resultados.tests$estadisticos, data.frame(test = 'CCC', parametro = 'ccc', valor = ccc))
    
    # iii. QCOMHD
    qcomhd <- TestQCOMHD(x = x, probs = seq(from=0.01, to=0.99, by=0.1),
      numero.muestras = 501L, ajuste = ajuste)
    resultados.tests$estadisticos <- rbind(resultados.tests$estadisticos, qcomhd %>%
        tidyr::gather(key = para, value = valor, -q) %>%
        dplyr::mutate(test = 'QCOMHD', parametro = paste0(para, '-', q)) %>%
        dplyr::select(test, parametro, valor))
  }
  
  # 5. Devolver resultados
  return (resultados.tests)
}
# -----------------------------------------------------------------------------