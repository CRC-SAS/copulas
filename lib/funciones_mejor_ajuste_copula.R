# Funcion para determinar cual familia de copula ajusta correctamente y
# luegos aquellas que ajusten positivamente son rankeadas en funcion del AIC, 
# BIC, RME, RMSE y una validacion cruzada

DeterminarMejorCopula <- function(ajustes.copula = NULL) {
  
  # Inicializar objeto para guardar resultados
  mejor.copula <- list(familia = NA, copula = NA, parametros = list(), dependencia = NA,
    aic = NA, bic = NA, mre = NA, rmse = NA, validacion.cruzada = NA)
  for (i in 1:length(ajustes.copula)) {
    familia <- ajustes.copula[[i]]$copula$familia
    copula <- ajustes.copula[[i]]$copula$copula
    parametros <- ajustes.copula[[i]]$copula$parametro
    dependencia <- ajustes.copula[[i]]$copula$dependencia
    bondad <- ajustes.copula[[i]]$bondad.ajuste
    
    if (!is.null(bondad) && !is.null(bondad$pasa.tests) && bondad$pasa.tests) {
      aic.ajuste <- bondad$estadisticos %>% 
        dplyr::filter(test == "AIC") %>%
        dplyr::pull(valor)
      bic.ajuste <- bondad$estadisticos %>% 
        dplyr::filter(test == "BIC") %>%
        dplyr::pull(valor)
      mre.ajuste <- bondad$estadisticos %>% 
        dplyr::filter(test == "MRE") %>%
        dplyr::pull(valor)       
      rmse.ajuste <- bondad$estadisticos %>% 
        dplyr::filter(test == "RMSE") %>%
        dplyr::pull(valor)
      xv.ajuste <- bondad$estadisticos %>% 
        dplyr::filter(test == "Validacion cruzada") %>%
        dplyr::pull(valor)
      
      if (!is.na(mejor.copula$familia)) {
        es.mejor.copula <- (aic.ajuste < mejor.copula$aic) ||
          ((aic.ajuste == mejor.copula$aic) && (bic.ajuste < mejor.copula$bic)) ||
          ((aic.ajuste == mejor.copula$aic) && (bic.ajuste == mejor.copula$ccc) && (mre.ajuste < mejor.copula$mre)) ||
          ((aic.ajuste == mejor.copula$aic) && (bic.ajuste == mejor.copula$ccc) && (mre.ajuste == mejor.copula$mre) && (rmse.ajuste < mejor.copula$rmse)) ||
          ((aic.ajuste == mejor.copula$aic) && (bic.ajuste == mejor.copula$ccc) && (mre.ajuste == mejor.copula$mre) && (rmse.ajuste == mejor.copula$rmse) && (xv.ajuste == mejor.copula$validacion.cruzada))
        
      } else {
        es.mejor.copula <- TRUE
      }
      if (es.mejor.copula) {
        mejor.copula$familia <- familia 
        mejor.copula$copula <- copula
        mejor.copula$parametros <- parametros 
        mejor.copula$dependencia <- dependencia
        mejor.copula$aic <- aic.ajuste
        mejor.copula$bic <- bic.ajuste
        mejor.copula$mre <- mre.ajuste
        mejor.copula$rmse <- rmse.ajuste 
        mejor.copula$validacion.cruzada <- xv.ajuste 
      }
    }
  }
  
  return(mejor.copula)
}
