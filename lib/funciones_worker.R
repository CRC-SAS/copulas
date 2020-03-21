

DeterminarDependenciaEntreSeriesPerturbadas <- function(input.value, script, x.prima, y.prima) {
  # Obtener la ubicación y las variables de la copula, para cada serie perturbada 
  # (la distribución en única para cada par ubicación, variable)
  vcu <- input.value
  
  # Identificar la columna con el id de la ubicación (usualmente station_id, o point_id)
  id_column <- IdentificarIdColumn(vcu)
  
  # Informar estado de la ejecución
  script$info(glue::glue("Determinando dependencia de las series perturbadas que componen la copula ", 
                         "{vcu$variable_x}-{vcu$variable_y} para la ubicación = {vcu %>% dplyr::pull(!!id_column)}, ",
                         "distribucion = \"{vcu$distribucion}\", serie_perturbada = {vcu$n_serie_perturbada} "))
  
}


DeterminarEstacionariedadDeSeriesPerturbadas <- function(input.value, script, umbral.p.valor) {
   # Identificar la columna con el id de la ubicación (usualmente station_id, o point_id)
  id_column <- IdentificarIdColumn(input.value)
  
  # Informar estado de la ejecución
  script$info(glue::glue("Determinando estacionariedad de la serie perturbada para la ubicación = {input.value %>% dplyr::pull(!!id_column)}, ",
                         "variable = \"{input.value$variable}\", distribucion = \"{input.value$distribucion}\", ",
                         "serie_perturbada = {input.value$n_serie_perturbada} "))
  
  # Extraer x.prima
  x.prima <- input.value %>% dplyr::pull(serie_perturbada) %>% do.call("c", .)
  fechas  <- input.value %>% dplyr::pull(fechas_serie_perturbada) %>% do.call("c", .)
  
  # Calcular estacionariedad
  t0 <- proc.time()
  resultados.estacionariedad.x <- EsEstacionaria(x = x.prima, fecha = fechas, umbral.p.valor)
  t1 <- proc.time() - t0
  
  return (input.value %>% dplyr::mutate(es_estacionaria = resultados.estacionariedad.x$pasa.tests, 
                                        segundos_calc_estacionariedad = t1[["elapsed"]]))

}


AplicarMejorAjusteASeriesPerturbadas <- function(input.value, script, eventos.completos) {
  # Identificar la columna con el id de la ubicación (usualmente station_id, o point_id)
  id_column <- IdentificarIdColumn(input.value)
  
  # Informar estado de la ejecución
  script$info(glue::glue("Realizando mejor ajuste para la ubicación = {input.value %>% dplyr::pull(!!id_column)}, ",
                         "variable = \"{input.value$variable}\", distribucion = \"{input.value$distribucion}\", ",
                         "serie_perturbada = {input.value$n_serie_perturbada} "))
  
  
  # -----------------------------------------------------------------------------#
  # Paso 0: Determinar valores iniciales ----
  # -----------------------------------------------------------------------------#
  
  eventos.ubicacion <- eventos.completos %>%
    dplyr::filter(!!rlang::sym(id_column) == dplyr::pull(input.value, !!id_column) & tipo_evento == 'seco')
  
  x <- abs(dplyr::pull(eventos.ubicacion, input.value$variable))
  fechas <- dplyr::pull(eventos.ubicacion, fecha_inicio)
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso 1: Aplicar mejor ajuste univariado ----
  # -----------------------------------------------------------------------------#
  
  # Crear variable con ruido
  x.prima <- AgregarRuido(x, input.value$valor_minimo_deteccion)
  
  # Ajustar distribucion univariada utilizando el mejor ajuste
  if (input.value$mejor_ajuste == 'lmomentos') {
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
  
  resultado <- input.value %>% 
    dplyr::select(-valor_minimo_deteccion) %>%
    dplyr::mutate(serie_perturbada = list(x.prima), fechas_serie_perturbada = list(fechas),
                  parametros_ajuste_serie_perturbada = list(ajuste$parametros))
  
  return (resultado)
  
}


MejorAjusteUnivariadoUV <- function(input.value, script, ajustes_univariados) {
  # Obtener la ubicación y variable a ser analizados
  uv <- input.value
  
  # Identificar la columna con el id de la ubicación (usualmente station_id, o point_id)
  id_column <- IdentificarIdColumn(uv)
  
  # Informar estado de la ejecución
  script$info(glue::glue("Determinando mejor ajuste para la ubicación {uv %>% dplyr::pull(!!id_column)} ",
                         "y la variable \"{uv$variable}\""))
  
  
  # -----------------------------------------------------------------------------#
  # CONTROL: salir si uv$variables no existe en ajustes_univariados ----
  # -----------------------------------------------------------------------------#
  
  if(! uv$variable %in% ajustes_univariados$variable) {
    type_of_id_col <- typeof(dplyr::pull(uv,!!id_column))
    return (tibble::tibble(!!id_column := if(type_of_id_col == "integer") integer() else 
                                          if(type_of_id_col == "numeric") double() else 
                                          if(type_of_id_col == "logical") logical() else 
                                          if(type_of_id_col == "character") character() else
                                            character(),
      variable = character(), distribucion = character(),
      mejor_ajuste = character(), parametros = list(),
      rmse = double(), ccc = double(), cuantiles = double()))
  }
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso 0: Determinar valores iniciales ----
  # -----------------------------------------------------------------------------#
  
  ajustes <- ajustes_univariados %>% 
    dplyr::filter(!!rlang::sym(id_column) == dplyr::pull(uv, !!id_column), 
                  variable == uv$variable) %>%
    dplyr::select(-tidyselect::matches(id_column), -variable) %>%
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
  return (tibble::tibble(!!id_column := dplyr::pull(uv, !!id_column),
                         variable = dplyr::pull(uv, variable), 
                         distribucion = mejor.ajuste$distribucion,
                         mejor_ajuste = mejor.ajuste$metodo_ajuste,
                         parametros = list(mejor.ajuste$parametros),
                         rmse = mejor.ajuste$rmse,
                         ccc = mejor.ajuste$ccc,
                         cuantiles = mejor.ajuste$cuantiles))
}


AjusteUnivariadoUVD <- function(input.value, script, config, eventos.completos) {
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
  
  eventos.ubicacion <- eventos.completos %>%
    dplyr::filter(!!rlang::sym(id_column) == dplyr::pull(uvd, !!id_column) & tipo_evento == 'seco')
  
  x <- abs(dplyr::pull(eventos.ubicacion, uvd$variable))
  fechas <- dplyr::pull(eventos.ubicacion, fecha_inicio)
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso 1: Realizar ajuste univariado ----
  # -----------------------------------------------------------------------------#
  
  configuracion <- uvd %>% 
    dplyr::select(distribucion, funcion_ajuste_lmomentos, funcion_ajuste_maxima_verosimilitud)
  
  parametros.lmomentos            <- list(x = x, min.cantidad.valores = 50)
  parametros.maxima.verosimilitud <- list(x = x, min.cantidad.valores = 50, numero.muestras = NULL)
  
  ajuste.univariado <- AjusteUnivariadoConfig(x = x,
                                              umbral.p.valor = config$params$umbral.p.valor,
                                              configuracion, parametros.lmomentos,
                                              parametros.maxima.verosimilitud)
  
  # ------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------#
  # Paso FINAL: Retornar resultados ----
  # -----------------------------------------------------------------------------#
  return (tibble::tibble(!!id_column := dplyr::pull(uvd, !!id_column),
                         variable = dplyr::pull(uvd, variable),
                         distribucion = ajuste.univariado$distribucion,
                         lmomentos = list(ajuste.univariado$lmomentos),
                         maxima.verosimilitud = list(ajuste.univariado$maxima.verosimilitud)))
}
