

AplicarMejorAjusteACopulas <- function(input.value, script, copulas.ajustadas, 
                                       mejor.ajuste.univariado, mejor.ajuste.multivariado,
                                       eventos.completos){
  # Obtener la ubicación, cópula y número de serie perturbasda
  ucs <- input.value
  
  # Identificar la columna con el id de la ubicación (usualmente station_id, o point_id)
  id_column <- IdentificarIdColumn(ucs)
  
  # Informar estado de la ejecución
  script$info(glue::glue("Aplicando el mejor ajuste multivariado a la ",
                         "copula \"{ucs$variable_x}-{ucs$variable_y}\", ",
                         "ubicación = {ucs %>% dplyr::pull(!!id_column)}, ",
                         "tipo_serie = \"{ucs$tipo_serie}\", ",
                         "n_serie = {ucs$n_serie}"))
  
  # Extraer la familia que mejor ajusta esta cópula
  mejor_familia <- mejor.ajuste.multivariado %>% 
    dplyr::filter(!!rlang::sym(id_column) == dplyr::pull(input.value, !!id_column),
                  variable_x == input.value$variable_x, variable_y == input.value$variable_y) %>%
    dplyr::pull(familia)
  
  # Se seleccionan las copulas que fueron ajustas utilizando la que resultó ser la mejor familia
  copulas_ajustadas_familia <- copulas.ajustadas %>%
    dplyr::filter(!!rlang::sym(id_column) == dplyr::pull(input.value, !!id_column),
                  variable_x == input.value$variable_x, variable_y == input.value$variable_y,
                  familia == mejor_familia)
  
  # Se extraen los parámetros obtenidos al ajustar cada una de las series perturbadas
  # utilizando la mejor familia
  parametros_x_serie <- purrr::map2_dfr(
    .x = copulas_ajustadas_familia$n_serie,
    .y = copulas_ajustadas_familia$copula,
    .f = function(n_serie, copula) {
      return(tibble::tibble(n_serie = !!n_serie,
                            parametro = copula$parametro[[1]]))
    })
  
  # Se calcula la mediana de los parametros obtenidos en el paso anterior
  parametro_mejor_copula <- median(parametros_x_serie$parametro)
  
  # A partir de parametro promediado, se crea una nueva mejor cópula
  mejor_copula <- do.call(what = paste0(mejor_familia, "Copula"), 
                          args = list(param = parametro_mejor_copula, dim = 2))
  
  # Aquí se obtiene el mejor ajuste univariado para la variable x
  mejor.ajuste.dsitribucion.x <- mejor.ajuste.univariado %>%
    dplyr::filter(!!rlang::sym(id_column) == dplyr::pull(input.value, !!id_column), 
                  variable == input.value$variable_x) 
  
  # Aquí se obtiene el mejor ajuste univariado para la variable y
  mejor.ajuste.dsitribucion.y <- mejor.ajuste.univariado %>%
    dplyr::filter(!!rlang::sym(id_column) == dplyr::pull(input.value, !!id_column), 
                  variable == input.value$variable_y) 
  
  # Finalmente, realizar el ajuste multivariado
  if (!is.na(mejor.ajuste.dsitribucion.x$distribucion) & !is.na(mejor.ajuste.dsitribucion.y$distribucion)) {
    ajuste_multivariado <- copula::mvdc(copula=mejor_copula,
                                        margins=c(mejor.ajuste.dsitribucion.x$distribucion, 
                                                  mejor.ajuste.dsitribucion.y$distribucion),
                                        paramMargins=list(mejor.ajuste.dsitribucion.x$parametros[[1]], 
                                                          mejor.ajuste.dsitribucion.y$parametros[[1]]))
  } else {
    ajuste_multivariado <- NA
  }

  return(input.value %>% 
           dplyr::mutate(parametro_mejor_copula = !!parametro_mejor_copula,
                         mejor_copula = list(mejor_copula),
                         mejor_ajuste_x_dist = mejor.ajuste.dsitribucion.x$distribucion,
                         mejor_ajuste_x_params = list(mejor.ajuste.dsitribucion.x$parametros[[1]]),
                         mejor_ajuste_y_dist = mejor.ajuste.dsitribucion.y$distribucion,
                         mejor_ajuste_y_params = list(mejor.ajuste.dsitribucion.y$parametros[[1]]),
                         ajuste_multivariado = list(ajuste_multivariado)))
  
}


MejorAjusteMultivariadoUC <- function(input.value, script, copulas.ajustadas, umbral.p.valor){
  # Obtener la ubicación y cópula a analizar
  uc <- input.value
  
  # Identificar la columna con el id de la ubicación (usualmente station_id, o point_id)
  id_column <- IdentificarIdColumn(uc)
  
  # Informar estado de la ejecución
  script$info(glue::glue("Determinando mejor ajuste multivariado para la ",
                         "copula \"{uc$variable_x}-{uc$variable_y}\", ", 
                         "ubicación = {uc %>% dplyr::pull(!!id_column)}"))

  # Se seleccionan las copulas de la familia en input_value
  copulas_ajustadas_ubicacion <- copulas.ajustadas %>%
    dplyr::filter(!!rlang::sym(id_column) == dplyr::pull(input.value, !!id_column),
                  variable_x == input.value$variable_x, variable_y == input.value$variable_y)
  
  # Se determinan las mátricas para la identificación de la cópula que mejor ajusta
  estadisticos <-  copulas_ajustadas_ubicacion %>% purrr::pmap_dfr(
    function(n_serie_perturbada, familia, bondad.ajuste, ...) {
      return(tidyr::crossing(n_serie_perturbada = n_serie_perturbada, familia = familia, 
                             bondad.ajuste$estadisticos))
    })
  metricas.ajuste.copula <- estadisticos %>%
    dplyr::filter(., test != 'estimate') %>%
    dplyr::group_by(familia, test) %>%
    dplyr::summarise(., mediana = median(valor),
                     media = mean(valor),
                     sd = sd(valor)) 
  
  # Inicializar objeto para guardar resultados
  mejor.ajuste <- tibble::tibble(
    familia = NA, AIC = Inf, BIC = Inf, MRE = Inf,
    RMSE = Inf, Sn = umbral.p.valor, Vx = Inf)
  
  # Identificación de la cópula que mejor ajusta
  resultado <- purrr::walk(
    .x = unique(metricas.ajuste.copula$familia),
    .f = function(familia, metricas.ajuste.copula) {
      
      bondad <- metricas.ajuste.copula %>% 
        dplyr::filter(familia == !!familia)
      
      if (!is.null(bondad)) {
        aic.ajuste <- bondad %>%
          dplyr::filter(test == "AIC") %>%
          dplyr::pull(mediana)
        bic.ajuste <- bondad %>%
          dplyr::filter(test == "BIC") %>%
          dplyr::pull(mediana)
        mre.ajuste <- bondad %>%
          dplyr::filter(test == "MRE") %>%
          dplyr::pull(mediana)
        rmse.ajuste <- bondad %>%
          dplyr::filter(test == "RMSE") %>%
          dplyr::pull(mediana)
        Sn.ajuste <- bondad %>%
          dplyr::filter(test == "Sn") %>%
          dplyr::pull(mediana)
        Vx.ajuste <- bondad %>%
          dplyr::filter(test == "Validacion cruzada") %>%
          dplyr::pull(mediana)
        if (!is.na(mejor.ajuste$familia)) {
          es.mejor.ajuste <- (Sn.ajuste > mejor.ajuste$Sn) || 
            ((Sn.ajuste == mejor.ajuste$Sn) && (aic.ajuste <- mejor.ajuste$AIC)) ||
            ((Sn.ajuste == mejor.ajuste$Sn) && (aic.ajuste == mejor.ajuste$AIC) && (bic.ajuste < mejor.ajuste$BIC)) || 
            ((Sn.ajuste == mejor.ajuste$Sn) && (aic.ajuste == mejor.ajuste$AIC) && (bic.ajuste == mejor.ajuste$BIC) && (mre.ajuste < mejor.ajuste$MRE)) ||
            ((Sn.ajuste == mejor.ajuste$Sn) && (aic.ajuste == mejor.ajuste$AIC) && (bic.ajuste == mejor.ajuste$BIC) && (mre.ajuste == mejor.ajuste$MRE) && (rmse.ajuste < mejor.ajuste$RMSE)) || 
            ((Sn.ajuste == mejor.ajuste$Sn) && (aic.ajuste == mejor.ajuste$AIC) && (bic.ajuste == mejor.ajuste$BIC) && (mre.ajuste == mejor.ajuste$MRE) && (rmse.ajuste == mejor.ajuste$RMSE) && (Vx.ajuste <- mejor.ajuste$Vx)) 
        } else {
          es.mejor.ajuste <- TRUE
        }
        
        if (es.mejor.ajuste) {
          mejor.ajuste$familia <<- familia
          mejor.ajuste$AIC <<- aic.ajuste
          mejor.ajuste$BIC <<- bic.ajuste
          mejor.ajuste$MRE <<- mre.ajuste
          mejor.ajuste$RMSE <<- rmse.ajuste
          mejor.ajuste$Sn <<- Sn.ajuste
          mejor.ajuste$Vx <<- Vx.ajuste
          
        }
      }
      
      return(mejor.ajuste)
      
    }, metricas.ajuste.copula
  )
  
  return(input.value %>% dplyr::bind_cols(mejor.ajuste))
    
}


AjustarCopulas <- function(input.value, script, eventos.completos, umbral.p.valor) {
  # Obtener la ubicación y las variables de la copula, para cada serie perturbada 
  # (la distribución es única -la que mejor ajustó- para cada par ubicación, variable)
  vcu <- input.value
  
  # Identificar la columna con el id de la ubicación (usualmente station_id, o point_id)
  id_column <- IdentificarIdColumn(vcu)
  
  # Informar estado de la ejecución
  script$info(glue::glue("Ajustando cópula \"{vcu$variable_x}-{vcu$variable_y}\" ", 
                         "para la ubicación = {vcu %>% dplyr::pull(!!id_column)} ({vcu$nombre}) ",
                         "tipo_serie = {vcu$tipo_serie}, n_serie = {vcu$n_serie}, ",
                         "familia = {vcu$familia}, funcion_ajuste = {vcu$funcion_ajuste} "))
  
  
  # -----------------------------------------------------------------------------#
  # Control INICIAL: salir si las variables en la copula no son dependientes  ----
  # -----------------------------------------------------------------------------#
  
  if(! vcu$son_dependientes) 
    return(input.value %>% dplyr::mutate(copula = list(NA), bondad.ajuste = list(NA)))
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso 0: Determinar valores iniciales ----
  # -----------------------------------------------------------------------------#
  
  eventos.ubic.var <- eventos.completos %>%
    dplyr::filter(!!rlang::sym(id_column) == dplyr::pull(input.value, !!id_column),
                  variable == input.value$variable_x | variable == input.value$variable_y, 
                  tipo_serie == vcu$tipo_serie, n_serie == input.value$n_serie) 
  
  x.prima <- eventos.ubic.var %>% dplyr::filter(variable == input.value$variable_x) %>% dplyr::pull(valor)
  y.prima <- eventos.ubic.var %>% dplyr::filter(variable == input.value$variable_y) %>% dplyr::pull(valor)
  
  # Crear parametros para el ajuste de las copulas
  parametros <- list(x = x.prima, y = y.prima)
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso 1: Realizar ajuste multivariado ----
  # -----------------------------------------------------------------------------#
  # 
  ajuste.copula <- do.call(what = vcu$funcion_ajuste, args = parametros)
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso 2: Realizar test de bondad del ajuste ----
  # -----------------------------------------------------------------------------#
  
  # Calcular bondad de ajuste
  bondad.ajuste.copula <- TestearBondadAjusteCopulas(x = parametros$x,
                                                     y = parametros$y,
                                                     umbral.p.valor = umbral.p.valor,
                                                     copula = ajuste.copula)
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso FINAL: Retornar resultados ----
  # -----------------------------------------------------------------------------#
  
  return (input.value %>% dplyr::mutate(copula = list(ajuste.copula), 
                                        bondad.ajuste = list(bondad.ajuste.copula)))
  
}


DeterminarDependencia <- function(input.value, script, eventos.completos,
                                  tests.dependencia, umbral.p.valor) {
  # Obtener la ubicación y las variables de la copula, para cada serie perturbada 
  # (la distribución es única -la que mejor ajustó- para cada par ubicación, variable)
  vcu <- input.value
  
  # Identificar la columna con el id de la ubicación (usualmente station_id, o point_id)
  id_column <- IdentificarIdColumn(vcu)
  
  # Informar estado de la ejecución
  script$info(glue::glue("Determinando dependencia de las series {vcu$tipo_serie}s que componen la copula ", 
                         "\"{vcu$variable_x}-{vcu$variable_y}\" para la ubicación = {vcu %>% dplyr::pull(!!id_column)} ",
                         "({vcu$nombre}) y serie_perturbada = {vcu$n_serie} "))
  
  
  # -----------------------------------------------------------------------------#
  # Paso 0: Determinar valores iniciales ----
  # -----------------------------------------------------------------------------#
  
  eventos.ubic.var <- eventos.completos %>%
    dplyr::filter(!!rlang::sym(id_column) == dplyr::pull(input.value, !!id_column),
                  variable == input.value$variable_x | variable == input.value$variable_y, 
                  tipo_serie == vcu$tipo_serie, n_serie == input.value$n_serie) 
  
  x.prima <- eventos.ubic.var %>% dplyr::filter(variable == input.value$variable_x) %>% dplyr::pull(valor)
  y.prima <- eventos.ubic.var %>% dplyr::filter(variable == input.value$variable_y) %>% dplyr::pull(valor)
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso 1: Calcular dependencia ----
  # -----------------------------------------------------------------------------#
  
  t0 <- proc.time()
  resultados.dependencia <- SonDependientes(x.prima, y.prima, tests.dependencia, umbral.p.valor)
  t1 <- proc.time() - t0
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso FINAL: Retornar resultados ----
  # -----------------------------------------------------------------------------#
  
  return (input.value %>% dplyr::mutate(son_dependientes = all(resultados.dependencia$pasa.test),
                                        detalles_dependencia = list(resultados.dependencia),
                                        segundos_calc_dependencia = t1[["elapsed"]]))
  
}


DeterminarEstacionariedad <- function(input.value, script, eventos.completos,
                                      tests.estacionariedad, umbral.p.valor) {
  # Identificar la columna con el id de la ubicación (usualmente station_id, o point_id)
  id_column <- IdentificarIdColumn(input.value)
  
  # Informar estado de la ejecución
  script$info(glue::glue("Determinando estacionariedad de la serie {input.value$tipo_serie} número {input.value$n_serie}, ",
                         "para la ubicación = {input.value %>% dplyr::pull(!!id_column)} ({input.value$nombre}), ",
                         "variable = \"{input.value$variable}\", realización = {input.value$realizacion}."))
  
  
  # -----------------------------------------------------------------------------#
  # Paso 0: Determinar valores iniciales ----
  # -----------------------------------------------------------------------------#
  
  eventos.ubic.var <- eventos.completos %>%
    dplyr::filter(!!rlang::sym(id_column) == dplyr::pull(input.value, !!id_column),
                  variable == input.value$variable, tipo_serie == input.value$tipo_serie, 
                  n_serie == input.value$n_serie, realizacion == input.value$realizacion) 
  
  x.prima <- dplyr::pull(eventos.ubic.var, valor)
  fechas <- dplyr::pull(eventos.ubic.var, fecha_inicio)
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso 1: Calcular estacionariedad ----
  # -----------------------------------------------------------------------------#
  
  t0 <- proc.time()
  resultados.estacionariedad.x <- EsEstacionaria(x.prima, fechas, tests.estacionariedad, umbral.p.valor)
  t1 <- proc.time() - t0
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso FINAL: Retornar resultados ----
  # -----------------------------------------------------------------------------#
  
  return (input.value %>% dplyr::mutate(es_estacionaria = all(resultados.estacionariedad.x$pasa.test),
                                        detalles_estacionariedad = list(resultados.estacionariedad.x),
                                        segundos_calc_estacionariedad = t1[["elapsed"]]))
  
}


AplicarMejorAjusteASeriesPerturbadas <- function(input.value, script, series.perturbadas) {
  # Identificar la columna con el id de la ubicación (usualmente station_id, o point_id)
  id_column <- IdentificarIdColumn(input.value)
  
  # Informar estado de la ejecución
  script$info(glue::glue("Realizando mejor ajuste para la ubicación = {input.value %>% dplyr::pull(!!id_column)}, ",
                         "variable = \"{input.value$variable}\", distribucion = \"{input.value$distribucion}\", ",
                         "tipo_serie = \"{input.value$tipo_serie}\", n_serie = {input.value$n_serie} "))
  
  
  # -----------------------------------------------------------------------------#
  # Paso 0: Determinar valores iniciales ----
  # -----------------------------------------------------------------------------#
  
  series.perturbadas.ubic.var <- series.perturbadas %>%
    dplyr::filter(!!rlang::sym(id_column) == dplyr::pull(input.value, !!id_column),
                  variable == input.value$variable, tipo_serie == "perturbada", 
                  n_serie == input.value$n_serie) 
  
  x.prima <- dplyr::pull(series.perturbadas.ubic.var, valor)
  fechas <- dplyr::pull(series.perturbadas.ubic.var, fecha_inicio)
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso 1: Aplicar mejor ajuste univariado ----
  # -----------------------------------------------------------------------------#
  
  # Ajustar distribucion univariada utilizando el mejor ajuste
  if (is.na(input.value$mejor_ajuste)) {
    ajuste <- list(parametros = NA)
  } else if (input.value$mejor_ajuste == 'lmomentos') {
    # Ajuste por L-momentos
    ajuste <- do.call(what = input.value$funcion_mejor_ajuste, 
                      args = list(x.prima, min.cantidad.valores = 50))
    
  } else if (input.value$mejor_ajuste == 'maxima.verosimilitud')   { 
    # Ajustar por Maxima Verosimilitud 
    ajuste <- do.call(what = input.value$funcion_mejor_ajuste, 
                      args = list(x.prima, min.cantidad.valores = 50))
  }
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso FINAL: Retornar resultados ----
  # -----------------------------------------------------------------------------#
  
  resultado <- input.value %>% dplyr::select(-valor_minimo_deteccion) %>%
    dplyr::mutate(parametros_ajuste_serie_perturbada = list(ajuste$parametros))
  
  return (resultado)
  
}


MejorAjusteUnivariadoUV <- function(input.value, script, ajustes.univariados) {
  # Obtener la ubicación y variable a analizar (ojo, no cópula, sino variable)
  uv <- input.value
  
  # Identificar la columna con el id de la ubicación (usualmente station_id, o point_id)
  id_column <- IdentificarIdColumn(uv)
  
  # Informar estado de la ejecución
  script$info(glue::glue("Determinando mejor ajuste para la ubicación {uv %>% dplyr::pull(!!id_column)} ",
                         "y la variable \"{uv$variable}\""))
  
  # -----------------------------------------------------------------------------#
  # Paso 0: Determinar valores iniciales ----
  # -----------------------------------------------------------------------------#
  
  ajustes <- ajustes.univariados %>% 
    dplyr::filter(!!rlang::sym(id_column) == dplyr::pull(uv, !!id_column), 
                  variable == uv$variable) %>%
    dplyr::select(distribucion, lmomentos, maxima.verosimilitud) %>%
    purrr::pmap(.f = function(...) { return (list(...)) })
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso 1: Determinar el mejor ajuste univariado ----
  # -----------------------------------------------------------------------------#
  
  mejor.ajuste   <- DeterminarMejorAjusteUnivariado(ajustes)
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso FINAL: Retornar resultados ----
  # -----------------------------------------------------------------------------#
  return (uv %>% dplyr::mutate(distribucion = mejor.ajuste$distribucion,
                               mejor_ajuste = mejor.ajuste$metodo_ajuste,
                               parametros = list(mejor.ajuste$parametros),
                               rmse = mejor.ajuste$rmse,
                               ccc = mejor.ajuste$ccc,
                               cuantiles = mejor.ajuste$cuantiles))
}


AjusteUnivariadoUVD <- function(input.value, script, serie.observada, umbral.p.valor) {
  # Ubicación, variable y distribución 
  uvd <- input.value
  
  # Identificar la columna con el id de la ubicación (usualmente station_id, o point_id)
  id_column <- IdentificarIdColumn(uvd)
  
  # Informar estado de la ejecución
  script$info(glue::glue("Realizando ajuste para la ubicación = {uvd %>% dplyr::pull(!!id_column)}, ",
                         "variable = \"{uvd$variable}\", distribucion = \"{uvd$distribucion}\""))
  
  
  # -----------------------------------------------------------------------------#
  # Paso 0: Determinar valores iniciales ----
  # -----------------------------------------------------------------------------#
  
  serie.observada.ubic.var <- serie.observada %>%
    dplyr::filter(!!rlang::sym(id_column) == dplyr::pull(uvd, !!id_column),
                  variable == uvd$variable, tipo_serie == "observada") 
  
  x <- dplyr::pull(serie.observada.ubic.var, valor)
  fechas <- dplyr::pull(serie.observada.ubic.var, fecha_inicio)
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso 1: Realizar ajuste univariado ----
  # -----------------------------------------------------------------------------#
  
  configuracion <- uvd %>% 
    dplyr::select(distribucion, funcion_ajuste_lmomentos, funcion_ajuste_maxima_verosimilitud)
  
  parametros.lmomentos            <- list(x = x, min.cantidad.valores = 50)
  parametros.maxima.verosimilitud <- list(x = x, min.cantidad.valores = 50, numero.muestras = NULL)
  
  ajuste.univariado <- AjusteUnivariadoConfig(x = x,
                                              umbral.p.valor = umbral.p.valor,
                                              configuracion, parametros.lmomentos,
                                              parametros.maxima.verosimilitud)
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso FINAL: Retornar resultados ----
  # -----------------------------------------------------------------------------#
  return (uvd %>% dplyr::mutate(lmomentos = list(ajuste.univariado$lmomentos),
                                maxima.verosimilitud = list(ajuste.univariado$maxima.verosimilitud)))
}
