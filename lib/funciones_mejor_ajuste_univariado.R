# Ordenar la combinación de distribución y paramétros en función del
# RMSE, CCC y el test de cuantiles

# Obtener resultados de test para ajuste con L-Momentos
DeterminarMejorAjusteUnivariado <- function(ajustes = NULL) {
  
  # Inicializar objeto para guardar resultados
  mejor.ajuste <- list(distribucion = NA, metodo_ajuste = NA, parametros = NA,
                       rmse = Inf, ccc = -Inf, cuantiles = 0)
  for (i in 1:length(ajustes)) {
    distribucion <- ajustes[[i]]$distribucion
    atributos <- names(ajustes[[i]])
    metodos_ajuste <- atributos[!atributos %in% "distribucion"]
    for (metodo_ajuste in metodos_ajuste) {
      parametros <- ajustes[[i]][[metodo_ajuste]]$parametros
      bondad <- ajustes[[i]][[metodo_ajuste]]$bondad
      if(!is.null(bondad) && !is.null(bondad$pasa.tests) && bondad$pasa.tests) {
        rmse.ajuste <- bondad$estadisticos %>% 
          dplyr::filter(test == "RMSEIQR") %>%
          dplyr::pull(valor)
        ccc.ajuste <- bondad$estadisticos %>% 
          dplyr::filter(test == "CCC") %>%
          dplyr::pull(valor)
        cuantiles.ajuste <- valor.cuantiles(bondad$estadisticos, cuantiles.evaluar = 0.5) 
        if (!is.na(mejor.ajuste$distribucion)) {
          es.mejor.ajuste <- (rmse.ajuste < mejor.ajuste$rmse) ||
            ((rmse.ajuste == mejor.ajuste$rmse) && (ccc.ajuste > mejor.ajuste$ccc)) ||
            ((rmse.ajuste == mejor.ajuste$rmse) && (ccc.ajuste == mejor.ajuste$ccc) && (cuantiles.ajuste > mejor.ajuste$cuantiles)) 
        } else {
          es.mejor.ajuste <- TRUE
        }
        if (es.mejor.ajuste) {
          mejor.ajuste$distribucion <- distribucion 
          mejor.ajuste$metodo_ajuste <- metodo_ajuste 
          mejor.ajuste$parametros <- parametros 
          mejor.ajuste$rmse <- rmse.ajuste 
          mejor.ajuste$ccc <- ccc.ajuste 
          mejor.ajuste$cuantiles <- cuantiles.ajuste 
        }
      }
    }
  }  
  return(mejor.ajuste)
}


# ---------------------------------------------------------------------------- #
# Función para normalizar cuantiles ----
# ---------------------------------------------------------------------------- #
# Función para extraer extraer un resumen del test de cuantiles
# A partir de los cuantiles elegido se hace un índice relativo en el que
# un valor de 1 indica que todos los cuantiles ajustan y un valor de 0 que ninguno
# no es estadísticamente diferente a los teóricos.
valor.cuantiles <- function(x, cuantiles.evaluar = 0.5) {
  purrr::map_df(x, .f = magrittr::extract) %>%
    dplyr::filter(., test == 'QCOMHD') %>%
    dplyr::select(., -test) %>%
    tidyr::separate(., col = parametro, into = c("parametros","cuantil"), sep = "-") %>%
    dplyr::filter(., cuantil >= cuantiles.evaluar) %>%
    tidyr::spread(., parametros, valor) %>%
    dplyr::mutate(., valor = if_else(p.value > p.crit, 1/nrow(.), 0)) %>%
    dplyr::summarise(., valor = sum(valor)) %>%
    dplyr::pull()
}


# ------------------------------------------------------------------------------
