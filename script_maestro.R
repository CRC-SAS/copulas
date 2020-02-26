# -----------------------------------------------------------------------------#
# Paso 1: Limpiar espacio de trabajo ----
rm(list = objects())
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# Paso 2: Cargar paquetes necesarios ----
require(dplyr)
require(purrr)
require(readxl)
require(lmomco)
require(goftest)
require(WRS2)
require(yardstick)
require(hydroGOF)
require(copula)
require(ggplot2)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# Paso 3: Cargar funciones necesarias ----
source('funciones_ajustes_marginales.R')
source('funciones_ajuste_distribuciones.R')
source('funciones_bondad_ajuste.R')
source('ajuste_distribuciones.R')
source('funciones_mejor_ajuste_univariado.R')
source('funciones_test_estacionaridad.R')
source('funciones_test_independencia.R')
source('funciones_ajuste_familias_copulas.R')
source('funciones_bondad_ajuste_copulas.R')
source('funciones_mejor_ajuste_copula.R')
source('funciones_auxiliares.R')
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# Paso 4: Leer archivo de configuracion ----

# Ajuste univariado 
# Leer archivo de configuracion de funciones de ajuste univariado
configuracion.ajuste.univariado <- readxl::read_excel(path = "funciones_ajuste.xls") %>%
  dplyr::filter(uso_sequias) %>%
  dplyr::select(distribucion, funcion_ajuste_lmomentos, funcion_ajuste_maxima_verosimilitud)

# Leer archivo de configuracion de funciones de ajuste multivariado
configuracion.ajuste.copula <- readxl::read_excel(path = "funciones_ajuste_copulas.xls") %>%
  dplyr::filter(uso_sequias) %>%
  dplyr::select(familia, funcion_ajuste)

umbral.p.valor  <- 0.05
n.realizaciones <- 1L
valor.minimo.deteccion.x <- 1
valor.minimo.deteccion.y <- 0.1
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# Paso 5: Conectar base de datos ----

# estacion.usar <- '87750'
estaciones <- c('87791')
escala <- 3L
variables <- c('duracion-minimo')

for (i in 1:length(estaciones)) {
  estacion.usar <- estaciones[i]

eventos.estacion <- eventos.completos %>%
  dplyr::filter(., estacion_id == estacion.usar & tipo_evento == 'seco')

x <- eventos.estacion$duracion
y <- eventos.estacion$minimo*-1
fechas <- eventos.estacion$fecha_inicio
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# Paso 6: Ajuste univariado ----

# Paso 6.a.: Ajuste de la variable a todas las distribuciones definidas

# Ajuste de la variable x
parametros.lmomentos            <- list(x = x, min.cantidad.valores = 50)
parametros.maxima.varosimilitud <- list(x = x, min.cantidad.valores = 50, numero.muestras = NULL)

ajuste.univariado.x <- AjusteUnivariado(x = x, p.valor = umbral.p.valor)



# Ajuste de la variable y
parametros.lmomentos            <- list(x = y, min.cantidad.valores = 50)
parametros.maxima.varosimilitud <- list(x = y, min.cantidad.valores = 50, numero.muestras = NULL)

ajuste.univariado.y <- AjusteUnivariado(x = y, p.valor = umbral.p.valor)

# ------------------------------------------------------------------------------

# Paso 6.b.: Determinar el mejor ajuste univariado ----

mejor.ajuste.x <- DeterminarMejorAjusteUnivariado(ajuste.univariado.x)

mejor.ajuste.y <- DeterminarMejorAjusteUnivariado(ajuste.univariado.y)
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# Paso 7. Graficos exploratorios ----

# Grafico de distribucion univariada de x
ggplot(data = as.data.frame(x), aes(x = x)) +
  stat_ecdf(geom = 'point') +
  stat_function(fun = paste0('p', mejor.ajuste.x$distribucion), args = mejor.ajuste.x$parametros) +
  theme_bw() +
  ggsave(paste0('distribucion.variable.x.', estacion.usar, '.', escala, '.', variables, '.png'), device = 'png', dpi = 600)

# Grafico de distribucion univariada de x
ggplot(data = as.data.frame(y), aes(x = y)) +
  stat_ecdf(geom = 'point') +
  stat_function(fun = paste0('p', mejor.ajuste.y$distribucion), args = mejor.ajuste.y$parametros) +
  theme_bw() +
  ggsave(paste0('distribucion.variable.y.', estacion.usar, '.', escala, '.', variables, '.png'), device = 'png', dpi = 600)
  
# Grafico exploratorio x e y
ggplot2::ggplot(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
  geom_point() +
  theme_bw() +
  ggsave(paste0('grafico.exploratorio.variables.originales.', estacion.usar, '.', escala, '.', variables, '.png'), device = 'png', dpi = 600)

# Grafico exploratorio x e y
ggplot2::ggplot(data = data.frame(x = copula::pobs(x), y = copula::pobs(y)), aes(x = x, y = y)) +
  geom_point() +
  theme_bw() +
  ggsave(paste0('grafico.exploratorio.variables.pesudoobservaciones.', estacion.usar, '.', escala, '.', variables, '.png'), device = 'png', dpi = 600)

# Grafico de Chi Plot
GraficarChiPlot(x = x, y = y)

# Grafico de Kendall
GraficarKPlot(x = x, y = y)

# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# Paso 8. Aleatorización de las variables y ajuste de copulas ----

# Paso 8.a.: Aleatorizacion y ajuste ----
# Inicializar objeto para guardar resultados
lista.ajuste.copulas <- list(realizacion = NA,
  variables.sinteticas = list(x = list(), y = list()),
  ajustes.univariados = list(x = list(), y = list()),
  son.dependientes = list(),
  son.estacionarias = list(), 
  ajustes.copulas = list())

configuracion.ajuste.univariado.x <- configuracion.ajuste.univariado %>%
  dplyr::filter(., distribucion == mejor.ajuste.x$distribucion) %>%
  dplyr::select(., distribucion, matches(mejor.ajuste.x$metodo_ajuste))

configuracion.ajuste.univariado.y <- configuracion.ajuste.univariado %>%
  dplyr::filter(., distribucion == mejor.ajuste.y$distribucion) %>%
  dplyr::select(., distribucion, matches(mejor.ajuste.y$metodo_ajuste))


# Ajuste univariado y evaluación de estacionaridad e independencia
for ( i in 1:n.realizaciones) {
  lista.ajuste.copulas$realizacion[[i]] <- i
  # Crear variable con ruido
  x.prima <- AgregarRuido(x, valor.minimo.deteccion.x)
  
  # Ajustar distribucion univariada a las n realizaciones
  # Ajuste de x.prima
  
  # Ajuste por L-momentos
  if (mejor.ajuste.x$metodo_ajuste == 'lmomentos') {
    ajuste <- do.call(what = configuracion.ajuste.univariado.x$funcion_ajuste_lmomentos, 
      args = list(x.prima, min.cantidad.valores = 50))
    
  } else if (mejor.ajuste.x$metodo_ajuste == 'maxima.verosimilitud')   { 
    # Ajustar por Maxima Verosimilitud 
    ajuste <- do.call(what = configuracion.ajuste.univariado.x$funcion_ajuste_maxima_verosimilitud, 
      args = list(x.prima, min.cantidad.valores = 50))
  }
  
  lista.ajuste.copulas$variables.sinteticas$x[[i]] <- x.prima
  lista.ajuste.copulas$ajustes.univariados$x[[i]] <- ajuste
  
  # Crear variable con ruido
  y.prima <- AgregarRuido(y, 0.1)
  
  # Ajustar distribucion univariada a las n realizaciones
  # Ajuste de y.prima
  
  # Ajuste por L-momentos
  if (mejor.ajuste.y$metodo_ajuste == 'lmomentos') {
    ajuste <- do.call(what = configuracion.ajuste.univariado.y$funcion_ajuste_lmomentos, 
      args = list(y.prima, min.cantidad.valores = 50))
    
  } else if (mejor.ajuste.y$metodo_ajuste == 'maxima.verosimilitud')   { 
    # Ajustar por Maxima Verosimilitud 
    ajuste <- do.call(what = configuracion.ajuste.univariado.y$funcion_ajuste_maxima_verosimilitud, 
      args = list(y.prima, min.cantidad.valores = 50))
  }
  
  lista.ajuste.copulas$variables.sinteticas$y[[i]] <- y.prima
  lista.ajuste.copulas$ajustes.univariados$y[[i]] <- ajuste
  
  # Evaluar dependencia de las variables sinteticas
  #resultados.dependencia <- SonDependientes(x.prima, y.prima)
  resultados.dependencia <- TRUE
  
  lista.ajuste.copulas$son.dependientes[i] <- resultados.dependencia
  
  # Evaluar estacionaridad de las variables sinteticas
  #resultados.estacionairdad.x <- EsEstacionaria(x = x.prima, fecha = fechas)
  #resultados.estacionairdad.y <- EsEstacionaria(x = y.prima, fecha = fechas)
  resultados.estacionairdad.x <- TRUE
  resultados.estacionairdad.y <- TRUE
  
  if (resultados.estacionairdad.x == TRUE && resultados.estacionairdad.y == TRUE) {
    lista.ajuste.copulas$son.estacionarias[[i]] <- TRUE
  } else {
    lista.ajuste.copulas$son.estacionarias[[i]] <- FALSE
  }
  
  if (lista.ajuste.copulas$son.dependientes[[i]] == TRUE && lista.ajuste.copulas$son.estacionarias[[i]] == TRUE) {
    
    # Crear parametros para el ajuste de las copulas
    parametros            <- list(x = lista.ajuste.copulas$variables.sinteticas$x[[i]], 
      y = lista.ajuste.copulas$variables.sinteticas$y[[i]])
    
    # Ajustar y calcular bondad de ajuste para cada familia
    ajustes.copula <- purrr::pmap(
      .l = configuracion.ajuste.copula, 
      .f = function(familia, funcion_ajuste) {
        # Ajustar y determinar la bondad del ajuste
        ajuste.copula <- NULL
        bondad.ajuste.copula <- NULL
        tryCatch({
          ajuste.copula <- do.call(what = funcion_ajuste, 
            args = parametros)
          
          bondad.ajuste.copula <- TestearBondadAjusteCopulas(x = parametros$x,
            y = parametros$y,
            umbral.p.valor = umbral.p.valor, 
            copula = ajuste.copula)
          
        }, error = function(e) {
          warning(paste0("Error al ejecutar de la familia ",
            familia, ": ", as.character(e)))
        })
        
        
        # Devolver lista con resultados
        resultados <- list(
          copula = ajuste.copula,
          bondad.ajuste = bondad.ajuste.copula)
        
        
        return (resultados)
      }
    )
    
    lista.ajuste.copulas$ajustes.copulas[[i]] <- ajustes.copula
  } else {
    lista.ajuste.copulas$ajustes.copulas[[i]] <- NA
  }
  
  print(paste0("Acaba de terminar la realización ", i, 'de la estacion ', estacion.usar))
}

save(lista.ajuste.copulas, file = paste0('lista.ajuste.copulas', '_estacion_', estacion.usar, '_', escala,'.', variables, '.rda'))

# Paso 8.b.: Extraer resultados por familia ----
# Inicializar objeto para guardar resultados
metricas.ajuste.copulas <- list(familia = NA, estadisticos = list())

for (j in 1:length(configuracion.ajuste.copula$familia)) {
  
  realizaciones.familia.j <- lista.ajuste.copulas$ajustes.copulas %>% rlist::list.map(.[[j]])
  
  # Inicializar objato para guardar resultados
  estadisticos <- NULL
  
  if (!is.null(realizaciones.familia.j[[i]]$bondad.ajuste$estadisticos)) {
    for (i in 1:n.realizaciones) {
      estadisticos <- rbind(estadisticos, realizaciones.familia.j[[i]]$bondad.ajuste$estadisticos)
    }
    medidas.resumen <- estadisticos %>%
      dplyr::filter(., test != 'estimate') %>%
      dplyr::group_by(test) %>%
      dplyr::summarise(., mediana = median(valor),
        media = mean(valor),
        sd = sd(valor)) 
    
    # guardar resultados
    metricas.ajuste.copulas$familia[j] <- configuracion.ajuste.copula$familia[j]
    metricas.ajuste.copulas$estadisticos[[j]] <- medidas.resumen
    
  } else {
    # guardar resultados
    metricas.ajuste.copulas$familia[j] <- configuracion.ajuste.copula$familia[j]
    metricas.ajuste.copulas$estadisticos[j] <- NA
    
  }
}

# Paso 8.c.: Graficos diagnosticos ----
# Crear data frame
metricas.ajuste.copulas.df <- purrr::pmap_dfr(
  .l = metricas.ajuste.copulas,
  .f = function(familia, estadisticos) {
    df <- NULL
    if (! is.na(estadisticos)) {
      df <- as.data.frame(estadisticos) %>%
        dplyr::mutate(familia := !! familia) %>%
        dplyr::select(familia, test, mediana, media, sd)
    }
    return (df)
  }
)

# Crear grafico diagnóstico de metricas

# Preparación del data frame
metricas.ajuste.copulas.df <- metricas.ajuste.copulas.df %>%
  dplyr::mutate(., familia = factor(metricas.ajuste.copulas.df$familia, 
    levels = c('clayton', 'frank', 'gumbel', 'joe', 'normal', 't'),
    labels = c('Clayton', 'Frank', 'Gumbel', 'Joe', 'Normal', 't'))) %>%
  dplyr::filter(., test != 'Sn')

ggplot2::ggplot(metricas.ajuste.copulas.df, aes(x = familia, y = mediana, fill = familia)) +
  ggplot2::geom_bar(stat = 'identity') +
  ggplot2::geom_errorbar(aes(ymin = mediana-sd, ymax = mediana + sd), width = .2,
    position = position_dodge(.9)) +
  scale_fill_brewer(palette = 'RdBu') +
  ggplot2::facet_wrap(~test, scales = 'free') +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none",
    axis.title.x=element_blank(),
    axis.title.y=element_blank()) +
  ggplot2::ggsave(filename = paste0('Graficos.diagnosticos.copulas.estacion.', estacion.usar, '.escala.', escala, '.', variables, '.png'), device = 'png',
    dpi = 600, width = 8, height = 5)

# Paso 8.c.: Determinar cual es la mejor familia de copulas
# Inicializar objeto para guardar resultados
mejor.ajuste <- list(familia = NA, AIC = Inf, BIC = Inf, MRE = Inf,
  RMSE = Inf, Sn = umbral.p.valor, Vx = Inf)

for (i in 1:length(configuracion.ajuste.copula$familia)) {
  familia <- metricas.ajuste.copulas$familia[i]
  bondad <- metricas.ajuste.copulas$estadisticos[[i]]
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
      mejor.ajuste$familia <- familia
      mejor.ajuste$AIC <- aic.ajuste
      mejor.ajuste$BIC <- bic.ajuste
      mejor.ajuste$MRE <- mre.ajuste
      mejor.ajuste$RMSE <- rmse.ajuste
      mejor.ajuste$Sn <- Sn.ajuste
      mejor.ajuste$Vx <- Vx.ajuste
      
    }
  }
}

# ---------------------------------------------------------------------------- #
# Pasoo 9: Evaluar efectos de la aleatorizacion ----
# ---------------------------------------------------------------------------- #
for (i in 1:length(configuracion.ajuste.copula$familia)) {
  if (configuracion.ajuste.copula$familia[i] == mejor.ajuste$familia) {
    familia <- i
  } 
}

# Extraer todas las realizaciones de la mejor familia
realizaciones.mejor.familia <- lista.ajuste.copulas$ajustes.copulas %>% rlist::list.map(.[[familia]])

# Distribución del p-valor de la mejor familia
aleatorizacion.p.valor <- purrr::map_dfr(
  .x = seq_len(length(realizaciones.mejor.familia)),
  .f = function(seq_index) {
    estadisticos <- realizaciones.mejor.familia[[seq_index]]$bondad.ajuste$estadisticos 
    return (estadisticos)
  }
) %>% dplyr::filter(., test == 'Sn' & parametro == 'p.value')

ggplot2::ggplot(aleatorizacion.p.valor, aes(x = valor)) +
  ggplot2::geom_histogram(bins = 11, fill="black", col="white") +
  labs(x = expression(p-valor), y = element_blank()) +
  theme_bw() + 
  ggsave(paste0("histograma.p.valor.mejor.familia.", estacion.usar, '.escala.', escala, '.', variables, ".png"),
    device = 'png', dpi = 600)

# ------------------------------------------------------------------------------

# Comprobar normalidad de las metricas
# parametro tau
# Extraer metricas de validacion
aleatorizacion.tau <- purrr::map_dfr(
  .x = seq_len(length(realizaciones.mejor.familia)),
  .f = function(seq_index) {
    tau.kendall <- realizaciones.mejor.familia[[seq_index]]$copula$dependencia$tau.kendall
    return (data.frame(tau.kendall = tau.kendall))
  }
)

# Parametro tau
ajuste.normal.tau <- AjustarMaximaVerosimilitudNormal(x = aleatorizacion.tau$tau.kendall, 
  min.cantidad.valores = 10)

resultados.test <- TestKS(aleatorizacion.tau$tau.kendall, ajuste.normal.tau)

media.tau <- round(ajuste.normal.tau$parametros$mu, 3)
media.anotacion <- paste("mu ==", media.tau)

sigma.tau <- round(ajuste.normal.tau$parametros$sigma, 3)
sigma.anotacion <- paste("sigma ==", sigma.tau)

p.tau <- round(resultados.test$p.value, 3)
p.anotacion <- paste("p-valor ==", p.tau)

# Grafico diagnóstico
ggplot2::ggplot(aleatorizacion.tau, aes(x = tau.kendall)) +
  stat_ecdf(geom = "point") +
  stat_function(fun = pnor, args = ajuste.normal.tau$parametros) +
  theme_bw() +
  labs(x = expression(tau), y = element_blank()) +
  annotate('text', x = min(aleatorizacion.tau$tau.kendall), y = 0.96, label = as.character(media.anotacion), parse = T, hjust = 0) +
  annotate('text', x = min(aleatorizacion.tau$tau.kendall), y = 0.9, label = as.character(sigma.anotacion), parse = T, hjust = 0) +
  annotate('text', x = min(aleatorizacion.tau$tau.kendall), y = 0.84, label = as.character(p.anotacion), parse = T, hjust = 0) +
  ggsave(paste0('distribucion.parametro.dependencia.tau.estacion', estacion.usar, '.escala.', '.', variables, '.png'),
    device = "png", dpi = 600)

# Parametro del ajuste de la copula
# Extraer metricas de validacion
aleatorizacion.theta <- purrr::map_dfr(
  .x = seq_len(length(realizaciones.mejor.familia)),
  .f = function(seq_index) {
    theta <- realizaciones.mejor.familia[[seq_index]]$copula$parametro$theta
    return (data.frame(theta = theta))
  }
)

ggplot2::ggplot(aleatorizacion.theta, aes(x = theta)) +
  ggplot2::geom_histogram(bins = 11, fill="black", col="white") +
  labs(x = expression(theta), y = element_blank()) +
  theme_bw() +
  ggsave(paste0('histograma.parametro.theta.estacion.', estacion.usar, '.escala.', escala, '.', variables, '.png'),
    device = 'png', dpi = 600)

# --------------------------------------------------------------------------------

# Prueba de los parametros de la distribucion univariada de x

# Extraer las realizaciones de la distribucipn de la variable x
realizaciones.univariadas.x <- lista.ajuste.copulas$ajustes.univariados$x 

parametros.univariados.x <- purrr::map_dfr(
  .x = seq_len(length(realizaciones.univariadas.x)),
  .f = function(seq_index) {
    parametros <- realizaciones.univariadas.x[[seq_index]]$parametros
    return (data.frame(parametros))
  }
) 


# Graficar la distribución de cada parametro de la disttribución unnivariada de x
parametros = names(parametros.univariados.x) 

for (i in 1:length(parametros)) {
  
  # Ajustar distribución del paramétro a una normal
  ajuste.normal.parametro <- AjustarMaximaVerosimilitudNormal(x = parametros.univariados.x[,i], 
    min.cantidad.valores = 10)
  
  # Evaluar bondad del ajuste a la normal
  resultados.test <- TestKS(parametros.univariados.x[,i], ajuste.normal.parametro)
  
  # Extraer resultados para la anotación del grafico
  media.parametro <- round(ajuste.normal.parametro$parametros$mu, 3)
  media.anotacion <- paste("mu ==", media.parametro)
  
  sigma.parametro <- round(ajuste.normal.parametro$parametros$sigma, 3)
  sigma.anotacion <- paste("sigma ==", sigma.parametro)
  
  p.parametro <- round(resultados.test$p.value, 3)
  p.anotacion <- paste("p-valor ==", p.parametro)
  
  # Graficar distribucion del parametro
  ggplot2::ggplot(parametros.univariados.x, aes(x = !!sym(parametros[i]))) +
    stat_ecdf(geom = "point") +
    stat_function(fun = pnor, args = ajuste.normal.parametro$parametros) +
    theme_bw() +
    labs(x = element_blank(), y = element_blank()) +
    annotate('text', x = min(parametros.univariados.x[,i]), y = 0.96, label = as.character(media.anotacion), parse = T, hjust = 0) +
    annotate('text', x = min(parametros.univariados.x[,i]), y = 0.9, label = as.character(sigma.anotacion), parse = T, hjust = 0) +
    annotate('text', x = min(parametros.univariados.x[,i]), y = 0.84, label = as.character(p.anotacion), parse = T, hjust = 0) +
    ggsave(paste0('distribucion.parametro.univariado.x.', i,'.estacion.', estacion.usar, '.escala.', escala, '.', variables, '.png'), device = "png", dpi = 600)

}
# ------------------------------------------------------------------------------

# Prueba de los parametros de la distribucion univariada de y

# Extraer las realizaciones de la distribucipn de la variable y
realizaciones.univariadas.y <- lista.ajuste.copulas$ajustes.univariados$y 

parametros.univariados.y <- purrr::map_dfr(
  .x = seq_len(length(realizaciones.univariadas.y)),
  .f = function(seq_index) {
    parametros <- realizaciones.univariadas.y[[seq_index]]$parametros
    return (data.frame(parametros))
  }
) 

# Graficar la distribución de cada parametro de la disttribución unnivariada de y
parametros = names(parametros.univariados.y) 

for (i in 1:length(parametros)) {
  
  # Ajustar distribución del paramétro a una normal
  ajuste.normal.parametro <- AjustarMaximaVerosimilitudNormal(x = parametros.univariados.y[,i], 
    min.cantidad.valores = 10)
  
  # Evaluar bondad del ajuste a la normal
  resultados.test <- TestKS(parametros.univariados.y[,i], ajuste.normal.parametro)
  
  # Extraer resultados para la anotación del grafico
  media.parametro <- round(ajuste.normal.parametro$parametros$mu, 3)
  media.anotacion <- paste("mu ==", media.parametro)
  
  sigma.parametro <- round(ajuste.normal.parametro$parametros$sigma, 3)
  sigma.anotacion <- paste("sigma ==", sigma.parametro)
  
  p.parametro <- round(resultados.test$p.value, 3)
  p.anotacion <- paste("p-valor ==", p.parametro)
  
  # Graficar distribucion del parametro
  ggplot2::ggplot(parametros.univariados.y, aes(x = !!sym(parametros[i]))) +
    stat_ecdf(geom = "point") +
    stat_function(fun = pnor, args = ajuste.normal.parametro$parametros) +
    theme_bw() +
    labs(x = element_blank(), y = element_blank()) +
    annotate('text', x = min(parametros.univariados.y[,i]), y = 0.96, label = as.character(media.anotacion), parse = T, hjust = 0) +
    annotate('text', x = min(parametros.univariados.y[,i]), y = 0.9, label = as.character(sigma.anotacion), parse = T, hjust = 0) +
    annotate('text', x = min(parametros.univariados.y[,i]), y = 0.84, label = as.character(p.anotacion), parse = T, hjust = 0) +
    ggsave(paste0('distribucion.parametro.univariado.y.', i,'.estacion.', estacion.usar, '.escala.', escala, '.', variables, '.png'), device = "png", dpi = 600)
  
}

# ----------------------------------------------------------------------------- 
}
# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Paso 10: Analisis multivariado ----
# ---------------------------------------------------------------------------- #

realizaciones.mejor.familia[1]

mejor.copula <- realizaciones.mejor.familia[[1]]$copula$copula

variable.x <- lista.ajuste.copulas$variables.sinteticas$x[[1]]

mejor.ajuste.dsitribucion.x <- lista.ajuste.copulas$ajustes.univariados$x[[1]]

variable.y <- lista.ajuste.copulas$variables.sinteticas$y[[1]]

mejor.ajuste.dsitribucion.y <- lista.ajuste.copulas$ajustes.univariados$y[[1]]

myMvd1<-mvdc(copula=mejor.copula,margins=c(mejor.ajuste.dsitribucion.x$distribucion, mejor.ajuste.dsitribucion.y$distribucion),
  paramMargins=list(mejor.ajuste.dsitribucion.x$parametros, mejor.ajuste.dsitribucion.y$parametros))

contour(myMvd1,dMvdc,xlim=c(0,30),ylim=c(0,3),col="blue",
  main="Gumbel-copula",xlab="duration",ylab="intensity")
points(variable.x,variable.y,col="red")


# Crear valores de la copula empirica
u <- seq(0, 1, length.out = length(variable.x))
grid <- as.matrix(expand.grid(u1 = u, u2 = u))
val <- cbind(grid, z = C.n(grid, X = matrix(cbind(variable.x, variable.y), ncol = 2)))

# Crear valores de la copula teorica
copula.teorica <- as.data.frame(grid) %>%
  dplyr::mutate(., z = pCopula(grid, copula = mejor.copula))

# Graficar comparacion de la copula empirica y teórica
ggplot(data = as.data.frame(copula.teorica), aes(x = u1, y = u2, z = z, color = 'black')) +
  stat_contour() +
  geom_dl(aes(label=..level..), method="bottom.pieces", stat="contour") +
  geom_contour(data = as.data.frame(val), aes(x = u1, y = u2, z = z,  color = 'red')) +
  theme_bw() +
  scale_colour_manual(guide = 'legend', values = c('black' = 'black', 'red' = 'red'), 
    labels = c('Cópula teórica', 'Cópula empírica')) +
  theme(legend.position = 'bottom', legend.title = element_blank(), legend.text = element_text(size = 15)) +
  xlab('x') + ylab('y') 
  ggsave(paste0('Copula.teorica.empirica.estacion.', estacion.usar, '.escala.', escala, '.png'),
    device = 'png', dpi = 600)

datos <- as.data.frame(cbind(variable.x, variable.y))

ggplot(as.data.frame(val), aes(u1, u2, z = z)) +
  stat_contour()

