# ------------------------------------------------------------------------------#
# Definicion de funcion para ajustar y calcular bondad de ajuste para una    ----
# ÚNICA distribucion, utilizando metodo de L-Momentos y Maxima Verosimilitud ----
# ------------------------------------------------------------------------------#

AjusteUnivariadoConfig <- function(x, umbral.p.valor, configuracion, parametros.lmomentos, 
                                   parametros.maxima.verosimilitud) {
  
  # Ajustar por L-Momentos y determinar la bondad del ajuste
  ajuste.lmomentos <- NULL
  bondad.ajuste.lmomentos <- NULL
  tryCatch({
    ajuste.lmomentos <- do.call(what = configuracion$funcion_ajuste_lmomentos, 
                                args = parametros.lmomentos)
    bondad.ajuste.lmomentos <- TestearBondadAjuste(x = x, umbral.p.valor = umbral.p.valor, 
                                                   ajuste = ajuste.lmomentos)
  }, error = function(e) {
    warning(paste0("Error al ejecutar ajuste por L-Momentos de distribucion ",
                   configuracion$distribucion, ": ", as.character(e)))
  })
  
  # Ajustar por Maxima Verosimilitud y determinar la bondad del ajuste
  ajuste.maxima.verosimilitud <- NULL
  bondad.ajuste.maxima.verosimilitud <- NULL
  tryCatch({
    ajuste.maxima.verosimilitud <- do.call(what = configuracion$funcion_ajuste_maxima_verosimilitud, 
                                           args = parametros.maxima.verosimilitud)
    bondad.ajuste.maxima.verosimilitud <- TestearBondadAjuste(x = x, umbral.p.valor = umbral.p.valor, 
                                                              ajuste = ajuste.maxima.verosimilitud)  
  }, error = function(e) {
    warning(paste0("Error al ejecutar ajuste por Maxima Versomilitud de distribucion ",
                   configuracion$distribucion, ": ", as.character(e)))
  })
  
  # Devolver resultados
  resultados <- list(
    distribucion = configuracion$distribucion,
    lmomentos = list(
      parametros = ajuste.lmomentos$parametros,
      bondad = bondad.ajuste.lmomentos
    ),
    maxima.verosimilitud = list(
      parametros = ajuste.maxima.verosimilitud$parametros,
      bondad = bondad.ajuste.maxima.verosimilitud
    )
  )
  return (resultados)
}
# -------------------------------------------------------------------------------


# ------------------------------------------------------------------------------#
# Definicion de funcion para ajustar y calcular bondad de ajuste para TODAS  ----  
# las distribuciones, utilizando metodo de L-Momentos y Maxima Verosimilitud ----
# ------------------------------------------------------------------------------#

AjusteUnivariado <- function(x, p.valor, configuracion, parametros.lmomentos, 
                             parametros.maxima.verosimilitud) {
  
  ajustes.original <- purrr::pmap(
    .l = configuracion, 
    .f = function(distribucion, funcion_ajuste_lmomentos, funcion_ajuste_maxima_verosimilitud) {
      # Ajustar por L-Momentos y determinar la bondad del ajuste
      ajuste.lmomentos <- NULL
      bondad.ajuste.lmomentos <- NULL
      tryCatch({
        ajuste.lmomentos <- do.call(what = funcion_ajuste_lmomentos, 
          args = parametros.lmomentos)
        bondad.ajuste.lmomentos <- TestearBondadAjuste(x = x, 
          umbral.p.valor = umbral.p.valor, 
          ajuste = ajuste.lmomentos)  
      }, error = function(e) {
        warning(paste0("Error al ejecutar ajuste por L-Momentos de distribucion ",
          distribucion, ": ", as.character(e)))
      })
      
      # Ajustar por Maxima Verosimilitud y determinar la bondad del ajuste
      ajuste.maxima.verosimilitud <- NULL
      bondad.ajuste.maxima.verosimilitud <- NULL
      tryCatch({
        ajuste.maxima.verosimilitud <- do.call(what = funcion_ajuste_maxima_verosimilitud, 
          args = parametros.maxima.verosimilitud)
        bondad.ajuste.maxima.verosimilitud <- TestearBondadAjuste(x = x, 
          umbral.p.valor = umbral.p.valor, 
          ajuste = ajuste.maxima.verosimilitud)  
      }, error = function(e) {
        warning(paste0("Error al ejecutar ajuste por Maxima Versomilitud de distribucion ",
          distribucion, ": ", as.character(e)))
      })
      
      # Devolver lista con resultados
      resultados <- list(
        distribucion = distribucion,
        lmomentos = list(
          parametros = ajuste.lmomentos$parametros,
          bondad = bondad.ajuste.lmomentos
        ),
        maxima.verosimilitud = list(
          parametros = ajuste.maxima.verosimilitud$parametros,
          bondad = bondad.ajuste.maxima.verosimilitud
        )
      )
      return (resultados)
    }
  )
  
  return(ajustes.original)
  
}
# -------------------------------------------------------------------------------
