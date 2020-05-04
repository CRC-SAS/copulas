# -----------------------------------------------------------------------------#
# --- PASO 1. Cargar paquetes necesarios                                    ----
# -----------------------------------------------------------------------------#

rm(list = ls()); gc()
Sys.setenv(TZ = "UTC")
list.of.packages <- c("dplyr", "purrr", "lubridate", "magrittr", 
                      "lmomco", "stringr", "utils", "yaml", "goftest",
                      "WRS2", "futile.logger", "doSNOW", "foreach", 
                      "iterators", "snow", "yardstick", "hydroGOF", 
                      "copula", "ggplot2", "R6", "RPostgres")
for (pack in list.of.packages) {
  if (!require(pack, character.only = TRUE)) {
    stop(paste0("Paquete no encontrado: ", pack))
  }
}
rm(pack); gc()

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 2. Leer archivo de configuracion                                 ----
# -----------------------------------------------------------------------------#

normalize_dirnames <- function(dirnames) {
  if (is.atomic(dirnames)) 
    dirnames <- base::sub('/$', '', dirnames)
  if (!is.atomic(dirnames))
    for (nm in names(dirnames)) 
      dirnames[[nm]] <- normalize_dirnames(dirnames[[nm]])
    return (dirnames)
}

# a) YAML de configuracion del generador de estadísticas móviles
args <- base::commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  archivo.config <- args[1]
} else {
  # No vino el archivo de configuracion por linea de comandos. Utilizo un archivo default
  archivo.config <- paste0(getwd(), "/configuracion_copulas.yml")
}
if (! file.exists(archivo.config)) {
  stop(paste0("El archivo de configuración ", archivo.config, " no existe\n"))
} else {
  cat(paste0("Leyendo archivo de configuración ", archivo.config, "...\n"))
  config <- yaml::yaml.load_file(archivo.config)
  config$dir <- normalize_dirnames(config$dir)
}

# b) YAML de parametros del generador de estadísticas móviles
if (length(args) > 1) {
  archivo.params <- args[2]
} else {
  # No vino el archivo de configuracion por linea de comandos. Utilizo un archivo default
  archivo.params <- paste0(getwd(), "/parametros_copulas.yml")
}
if (! file.exists(archivo.params)) {
  stop(paste0("El archivo de parámetros ", archivo.params, " no existe\n"))
} else {
  cat(paste0("Leyendo archivo de parámetros ", archivo.params, "...\n"))
  config$params <- yaml::yaml.load_file(archivo.params)
}

replace_run_identifier <- function(filenames, identifier) {
  if (is.atomic(filenames)) 
    filenames <- base::sub('<\\*idc>', identifier, filenames)
  if (!is.atomic(filenames))
    for (nm in names(filenames)) 
      filenames[[nm]] <- replace_run_identifier(filenames[[nm]], identifier)
    return (filenames)
}

# c) YAML de configuración del intercambio de archivos del proceso de generación de índices
if (length(args) > 1) {
  archivo.nombres <- args[3]
} else {
  # No vino el archivo de configuracion por linea de comandos. Utilizo un archivo default
  archivo.nombres <- paste0(config$dir$data, "/configuracion_archivos_utilizados.yml")
}
if (! file.exists(archivo.nombres)) {
  stop(paste0("El archivo de configuración ", archivo.nombres, " no existe\n"))
} else {
  cat(paste0("Leyendo archivo de configuración ", archivo.nombres, "...\n"))
  config$files <- yaml::yaml.load_file(archivo.nombres)
  config$files <- replace_run_identifier(config$files, config$files$identificador_corrida)
}

rm(archivo.config, archivo.params, archivo.nombres, args); gc()

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 3. Cargar librerias propias e iniciar script                     ----
# -----------------------------------------------------------------------------#

# a) Cargar librerias
source(glue::glue("{config$dir$lib}/Script.R"), echo = FALSE)
source(glue::glue("{config$dir$lib}/Task.R"), echo = FALSE)
source(glue::glue("{config$dir$lib}/Helpers.R"), echo = FALSE)

# b) Cargar funciones necesarias 
source(glue::glue("{config$dir$base}/lib/funciones_ajustes_marginales.R"), echo = FALSE)
source(glue::glue("{config$dir$base}/lib/funciones_ajuste_distribuciones.R"), echo = FALSE)
source(glue::glue("{config$dir$base}/lib/funciones_bondad_ajuste.R"), echo = FALSE)
source(glue::glue("{config$dir$base}/lib/ajuste_distribuciones.R"), echo = FALSE)
source(glue::glue("{config$dir$base}/lib/funciones_mejor_ajuste_univariado.R"), echo = FALSE)
source(glue::glue("{config$dir$base}/lib/funciones_test_estacionaridad.R"), echo = FALSE)
source(glue::glue("{config$dir$base}/lib/funciones_test_independencia.R"), echo = FALSE)
source(glue::glue("{config$dir$base}/lib/funciones_ajuste_familias_copulas.R"), echo = FALSE)
source(glue::glue("{config$dir$base}/lib/funciones_bondad_ajuste_copulas.R"), echo = FALSE)
source(glue::glue("{config$dir$base}/lib/funciones_mejor_ajuste_copula.R"), echo = FALSE)
source(glue::glue("{config$dir$base}/lib/funciones_auxiliares.R"), echo = FALSE)
source(glue::glue("{config$dir$base}/lib/funciones_worker.R"), echo = FALSE)


# c.1) Definir nombre del script
script_name <- "CalcCopulas"
script_logfile <- glue::glue("{config$dir$run}/{script_name}.log")

# c.2) borrar archivo .log de corridas anteriores
if (file.exists(script_logfile))
  file.remove(script_logfile)

# c.3) Iniciar script
script <- Script$new(run.dir = config$dir$run, name = script_name, create.appender = T)
script$start()

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 4. Leer/procesar variables de entrada                            ----
# -----------------------------------------------------------------------------#

## CONFIGURACIÓN AJUSTES

# Ajuste univariado 
# Leer archivo de configuracion de funciones de ajuste univariado
configuracion.ajuste.univariado.completo <- purrr::pmap_dfr(
  .l = utils::tail(config$params$configuracion.ajuste.univariado, -1) %>% purrr::transpose(), 
  .f = function(..., nombres_columnas) { 
    tibble::tibble(..., .name_repair = "minimal") %>% setNames(nombres_columnas) },
  nombres_columnas = config$params$configuracion.ajuste.univariado[[1]]
)
# Tomar solo configuraciones con uso_sequias igual a True
configuracion.ajuste.univariado <- configuracion.ajuste.univariado.completo %>% dplyr::filter(uso_sequias) %>%
  dplyr::select(distribucion, funcion_ajuste_lmomentos, funcion_ajuste_maxima_verosimilitud)

# Ajuste multivariado 
# Leer archivo de configuracion de funciones de ajuste multivariado
configuracion.ajuste.copula.completo <- purrr::pmap_dfr(
  .l = utils::tail(config$params$configuracion.ajuste.copula, -1) %>% purrr::transpose(), 
  .f = function(..., nombres_columnas) { 
    tibble::tibble(..., .name_repair = "minimal") %>% setNames(nombres_columnas) },
  nombres_columnas = config$params$configuracion.ajuste.copula[[1]]
)
# Tomar solo configuraciones con uso_sequias igual a True
configuracion.ajuste.copula <- configuracion.ajuste.copula.completo %>% dplyr::filter(uso_sequias) %>%
  dplyr::select(familia, funcion_ajuste)


## PARÁMETROS RESTANTES

# Obtener la configuración a ser analizada
configuraciones_eventos <- purrr::pmap_dfr(
  .l = utils::tail(config$params$configuraciones.eventos, -1) %>% purrr::transpose(), 
  .f = function(..., nombres_columnas) { 
    tibble::tibble(..., .name_repair = "minimal") %>% setNames(nombres_columnas) },
  nombres_columnas = config$params$configuraciones.eventos[[1]]
)
if (nrow(configuraciones_eventos) > 1)
  stop("Aún no se implementó la funcionalidad que permite procesar varias configuraciones en paralelo.")

# Obtener los eventos a ser analizados!
eventos <- feather::read_feather(glue::glue("{config$dir$data}/{config$files$eventos_identificados}")) %>%
  dplyr::filter(conf_id %in% configuraciones_eventos$conf_id) %>% 
  dplyr::filter(realizacion <= config$params$n.realizaciones) %>%
  dplyr::filter(tipo_evento == config$params$eventos$tipo, 
                duracion >= config$params$eventos$duracion_minima)
id_column <- IdentificarIdColumn(eventos)
# Controlar que eventos tenga un id identificable
if (! any(c("station_id", "point_id") %in% colnames(eventos)))
  stop("El archivo con los eventos a ser identificados (indicado por el parámetro: eventos_identificados) ",
       "debe contnener un tibble con alguna de las siguientes columnas como id: \"station_id\" o \"point_id\"!!")

# Obtener ubicaciones a ser consideradas (deben existir en eventos!)
ubicaciones <- purrr::pmap_dfr(
  .l = config$params$ubicaciones %>% purrr::transpose(), 
  .f = function(...) { tibble::tibble(..., .name_repair = "minimal") %>% 
      setNames(c(id_column, 'nombre')) }
)
# Controlar que las estaciones seleccionadas estén entre aquellas 
# para las que se identificaron eventos en los pasos previos!
ids_in_params <- sapply(config$params$ubicaciones, function(u) {u$id})
ids_in_events <- eventos %>% dplyr::pull(id_column) %>% as.character()
if (!all(ids_in_params %in% ids_in_events))
  stop("No hay datos, en eventos, para todas las ubicaciones en parametros_copulas.yml")

# Obtener variables a ser consideradas (deben existir en eventos!)
variables_copulas <- purrr::pmap_dfr(
  .l = config$params$variables_copulas %>% purrr::transpose(), 
  .f = function(...) { tibble::tibble(..., .name_repair = "minimal") %>% 
      setNames(c('variable_x', 'variable_y')) }
)
# Controlar que las variables seleccionadas sean columnas de eventos!
if (!all(union(variables_copulas$variable_x, variables_copulas$variable_y) %in% colnames(eventos)))
  stop("No hay datos, en eventos, para todas las variables en variables_copulas en parametros_copulas.yml")

valores_minimos <- purrr::pmap_dfr(
  .l = config$params$valores_minimos %>% purrr::transpose(), 
  .f = function(...) { tibble::tibble(..., .name_repair = "minimal") %>% 
      setNames(c('variable', 'valor_minimo_deteccion')) }
)
# Controlar que hayan valores mínimos para todas las variables!
if (!all(union(variables_copulas$variable_x, variables_copulas$variable_y) %in% valores_minimos$variable))
  stop("No hay datos, en eventos, para todas las variables en valores_minimos en parametros_copulas.yml")

# A partir de los eventos y variables_copulas se genera eventos_completos. Este tibble, a diferencia 
# de eventos, tiene el valor absoluto de las variables y , además, todas las series perturbadas a ser  
# utilizadas en los análisis, de manera a optimizar el uso de memoria para grandes cantidades de eventos. 
serie_observada <- eventos %>%
  dplyr::mutate(intensidad = abs(intensidad), magnitud = abs(magnitud),
                duracion = abs(duracion), minimo = abs(minimo), maximo = abs(maximo)) %>%
  dplyr::select(realizacion, !!id_column, tipo_evento, conf_id, numero_evento, 
                fecha_inicio, intensidad, magnitud, duracion, minimo, maximo) %>%
  tidyr::pivot_longer(cols = c(intensidad, magnitud, duracion, minimo, maximo),
                      names_to = "variable", values_to = "valor") %>%
  dplyr::mutate(tipo_serie = "observada", n_serie = 0) %>% 
  dplyr::arrange(variable) 

# OBS: para el cálculo de las series perturbadas se utiliza la semilla definida en el 
# archivo de configuración de parámetros.
set.seed(config$params$s.series.perturbadas)

series_perturbadas <- purrr::map_dfr(
  .x = 1:config$params$n.series.perturbadas,
  .f = function(n_serie_perturbada, valores_minimos, serie_observada) {
    serie_perturbada <- serie_observada %>%
      dplyr::mutate(tipo_serie = "perturbada", n_serie = n_serie_perturbada) %>%
      dplyr::inner_join(valores_minimos, by = "variable") %>%
      dplyr::group_by(!!rlang::sym(id_column), tipo_evento, variable) %>%
      dplyr::mutate(valor = AgregarRuido(valor, unique(valor_minimo_deteccion))) %>%
      dplyr::ungroup() %>% dplyr::select(-valor_minimo_deteccion) %>%
      dplyr::arrange(realizacion, !!rlang::sym(id_column), variable, desc(tipo_evento), conf_id, numero_evento)
    return(serie_perturbada)
  }, valores_minimos, serie_observada)

# Finalmente, se crea el tibble eventos_completos!!
eventos_completos <- dplyr::bind_rows(serie_observada, series_perturbadas)

# Se guarda una copia de los eventos completos!
eventos_filename <- glue::glue("{config$dir$data}/{config$files$copulas$eventos_completos}")
if (file.exists(eventos_filename))
  file.remove(eventos_filename)
feather::write_feather(eventos_completos, eventos_filename)

# Se guarda la semilla, para poder generar las mismas series perturbadas posteriormente
seed_filename <- glue::glue("{config$dir$data}/{config$files$copulas$semilla_utilizada}")
if (file.exists(seed_filename))
  file.remove(seed_filename)
saveRDS(config$params$s.series.perturbadas, seed_filename)

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 5. Ajustar, por separado, utilizando las series observadas, cada una 
# --- de las variables que conforman las cópulas (parámetro: variables_copulas) 
# --- a las distribuciones señaladas en el archivo de parámetros (parámetro: 
# --- configuracion.ajuste.univariado). 
# --- Clave primaria: ubicación, variable, distribución ----
# -----------------------------------------------------------------------------#


# Definir el objetos sobre los cuales iterar
ubicacion_x_variable_x_distribucion <- 
  tidyr::crossing(ubicaciones, variable = union(variables_copulas$variable_x, variables_copulas$variable_y),
                  configuracion.ajuste.univariado)


# Definir el nombre de la función a ser paralelizada
function_name <- "AjusteUnivariadoUVD"

# Definir nombre de archivos .log y .out de corridas anteriores
task_logfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.log")
task_outfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.out")

# Borrar archivos .log y .out de corridas anteriores
if (file.exists(task_logfile))
  file.remove(task_logfile)
if (file.exists(task_outfile))
  file.remove(task_outfile)

# Definir nombre del archivo donde se van a guardar los resultados
results_filename <- glue::glue("{config$dir$data}/{config$files$copulas$ajustes_univariados}")

# Borrar archivo de resultado de corridas anteriores
if (file.exists(results_filename))
  file.remove(results_filename)

# Crear tarea distribuida y ejecutarla
task <- Task$new(parent.script = script,
                 func.name = function_name,
                 packages = list.of.packages)

# Informar inicio de ejecución 
script$info("Computando ajuste univariado para cada combinación de ubicación, variable, distribución")
# Ejecutar tarea distribuida
ajuste.univariado.x.ubic.var.dist <- task$run(number.of.processes = config$max.procesos,
                                              input.values = ubicacion_x_variable_x_distribucion,  
                                              serie.observada = serie_observada,
                                              umbral.p.valor = config$params$umbral.p.valor)

# Transformar resultados a un objeto de tipo tibble
ajuste.univariado.x.ubic.var.dist <- ajuste.univariado.x.ubic.var.dist %>% purrr::map_dfr(~.x)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.errors <- task$getErrors()
if (length(task.errors) > 0) {
  for (error.obj in task.errors) {
    id_column <- IdentificarIdColumn(ubicacion_x_variable_x_distribucion[1,])
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}, ",
                           "variable=\"{error.obj$input.value[['variable']]}\", ",
                           "distribucion=\"{error.obj$input.value[['distribucion']]}\")",
                           ": {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
} else {
  # Guardar resultados en un archivo fácil de compartir
  script$info(glue::glue("Guardando ajustes univariados en el archivo {results_filename}"))
  saveRDS(ajuste.univariado.x.ubic.var.dist, results_filename)
}

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 6. Determinar la distribución que mejor se ajusta a cada una de las
# --- variables (utilizando series observadas) en cada una de las ubicaciones.
# --- Clave primaria: ubicación, variable  ----
# --- OBS: a partir de acá para cada variable, en cada ubicación, ya se sabe 
# --- cual es la distribución que mejor ajusta y el mejor método de ajuste
# -----------------------------------------------------------------------------#


# Definir el objetos sobre los cuales iterar
ubicacion_x_variable <- ajuste.univariado.x.ubic.var.dist %>%
  dplyr::select(!!id_column, variable) %>% dplyr::distinct()


# Definir el nombre de la función a ser paralelizada
function_name <- "MejorAjusteUnivariadoUV"

# Definir nombre de archivos .log y .out de corridas anteriores
task_logfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.log")
task_outfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.out")

# Borrar archivos .log y .out de corridas anteriores
if (file.exists(task_logfile))
  file.remove(task_logfile)
if (file.exists(task_outfile))
  file.remove(task_outfile)

# Definir nombre del archivo donde se van a guardar los resultados
results_filename <- glue::glue("{config$dir$data}/{config$files$copulas$mejores_ajustes_univariados}")

# Borrar archivo de resultado de corridas anteriores
if (file.exists(results_filename))
  file.remove(results_filename)

# Crear tarea distribuida y ejecutarla
task.mejor.ajuste <- Task$new(parent.script = script,
                              func.name = function_name,
                              packages = list.of.packages)

# Informar inicio de ejecución 
script$info("Determinando mejor ajuste univariado para cada combinación de ubicación, variable")
# Ejecutar tarea distribuida
mejor.ajuste.univariado.x.ubic.var <- task.mejor.ajuste$run(number.of.processes = config$max.procesos,
                                                            input.values = ubicacion_x_variable,  
                                                            ajustes.univariados = ajuste.univariado.x.ubic.var.dist)

# Transformar resultados a un objeto de tipo tibble
mejor.ajuste.univariado.x.ubic.var <- mejor.ajuste.univariado.x.ubic.var %>% purrr::map_dfr(~.x)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.mejor.ajuste.errors <- task.mejor.ajuste$getErrors()
if (length(task.mejor.ajuste.errors) > 0) {
  for (error.obj in task.mejor.ajuste.errors) {
    id_column <- IdentificarIdColumn(ubicacion_x_variable[1,])
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}, ",
                           "variable=\"{error.obj$input.value[['variable']]}\")",
                           ": {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
} else {
  # Guardar resultados en un archivo fácil de compartir
  script$info(glue::glue("Guardando mejores ajustes univariados en el archivo {results_filename}"))
  saveRDS(mejor.ajuste.univariado.x.ubic.var, results_filename)
}

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 7. Ajustar las series perturbadas utilizando la distribución y el 
# --- método de ajuste que mejor ajustaron las series observadas. 
# --- Clave primaria: ubicación, variable, serie_perturbada ----
# --- OBS: la distribución ya no forma parte de la clave primaria porque cada
# --- par ubicación, variable, ya tiene definida la mejor distribución y el 
# --- mejor método de ajuste!!
# -----------------------------------------------------------------------------#


# Definir el objetos sobre los cuales iterar
mejor_ajuste_x_serie_perturbada <- mejor.ajuste.univariado.x.ubic.var %>% 
  dplyr::select(-parametros, -rmse, -ccc, -cuantiles) %>%
  dplyr::left_join(configuracion.ajuste.univariado, by = "distribucion") %>%
  dplyr::mutate(funcion_mejor_ajuste = ifelse(mejor_ajuste == 'lmomentos', 
                                              funcion_ajuste_lmomentos, 
                                              funcion_ajuste_maxima_verosimilitud)) %>%
  dplyr::select(-funcion_ajuste_lmomentos, -funcion_ajuste_maxima_verosimilitud) %>%
  dplyr::mutate(tipo_serie = "perturbada") %>%
  tidyr::crossing(n_serie = 1:config$params$n.series.perturbadas) %>%
  dplyr::left_join(valores_minimos, by = "variable")


# Definir el nombre de la función a ser paralelizada
function_name <- "AplicarMejorAjusteASeriesPerturbadas"

# Definir nombre de archivos .log y .out de corridas anteriores
task_logfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.log")
task_outfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.out")

# Borrar archivos .log y .out de corridas anteriores
if (file.exists(task_logfile))
  file.remove(task_logfile)
if (file.exists(task_outfile))
  file.remove(task_outfile)

# Definir nombre del archivo donde se van a guardar los resultados
results_filename <- glue::glue("{config$dir$data}/{config$files$copulas$series_perturbadas_ajustadas}")

# Borrar archivo de resultado de corridas anteriores
if (file.exists(results_filename))
  file.remove(results_filename)

# Crear tarea distribuida y ejecutarla
task <- Task$new(parent.script = script,
                 func.name = function_name,
                 packages = list.of.packages)

# Informar inicio de ejecución 
script$info("Generar series perturbadas y aplicarles el mejor ajuste")
# Ejecutar tarea distribuida
series.perturbadas.ajustadas <- task$run(number.of.processes = config$max.procesos,
                                         input.values = mejor_ajuste_x_serie_perturbada,  
                                         series.perturbadas = series_perturbadas)

# Transformar resultados a un objeto de tipo tibble
series.perturbadas.ajustadas <- series.perturbadas.ajustadas %>% purrr::map_dfr(~.x)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.errors <- task$getErrors()
if (length(task.errors) > 0) {
  for (error.obj in task.errors) {
    id_column <- IdentificarIdColumn(mejor_ajuste_x_serie_perturbada[1,])
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}, ",
                           "variable=\"{error.obj$input.value[['variable']]}\", ",
                           "distribucion=\"{error.obj$input.value[['distribucion']]}\", ",
                           "tipo_serie=\"{error.obj$input.value[['tipo_serie']]}\", ",
                           "n_serie={error.obj$input.value[['n_serie']]})",
                           ": {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
} else {
  # Guardar resultados en un archivo fácil de compartir
  script$info(glue::glue("Guardando las series perturbadas y sus ajustes univariados en el archivo {results_filename}"))
  saveRDS(series.perturbadas.ajustadas, results_filename)
}

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 8. Determinar la estacionariedad de las series asociadas a las
# --- variables que conforman cada una de las cópulas.
# --- Clave primaria: ubicación, variable, serie_perturbada ----
# -----------------------------------------------------------------------------#


# Definir el objetos sobre los cuales iterar
ubicacion_x_variable <- ubicaciones %>%
  tidyr::crossing(variable = union(variables_copulas$variable_x, variables_copulas$variable_y)) 
series_observadas_a_computar <- ubicacion_x_variable %>% 
  dplyr::mutate(tipo_serie = "observada", n_serie = 0) %>%
  tidyr::crossing(realizacion = unique(eventos$realizacion))
series_perturbadas_a_computar <- ubicacion_x_variable %>% 
  dplyr::mutate(tipo_serie = "perturbada") %>%
  tidyr::crossing(n_serie = 1:config$params$n.series.perturbadas) %>%
  tidyr::crossing(realizacion = unique(eventos$realizacion))
series_a_computar <- 
  dplyr::bind_rows(series_observadas_a_computar, series_perturbadas_a_computar)


# Definir el nombre de la función a ser paralelizada
function_name <- "DeterminarEstacionariedad"

# Definir nombre de archivos .log y .out de corridas anteriores
task_logfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.log")
task_outfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.out")

# Borrar archivos .log y .out de corridas anteriores
if (file.exists(task_logfile))
  file.remove(task_logfile)
if (file.exists(task_outfile))
  file.remove(task_outfile)

# Definir nombre del archivo donde se van a guardar los resultados
results_filename <- glue::glue("{config$dir$data}/{config$files$copulas$estacionariedad}")

# Borrar archivo de resultado de corridas anteriores
if (file.exists(results_filename))
  file.remove(results_filename)

# Crear tarea distribuida y ejecutarla
task <- Task$new(parent.script = script,
                 func.name = function_name,
                 packages = list.of.packages)

# Informar inicio de ejecución 
script$info("Determinar estacionariedad, tanto de la serie observada como de las series perturbadas")
# Ejecutar tarea distribuida
estacionariedad <- task$run(number.of.processes = config$max.procesos,
                            input.values = series_a_computar, 
                            eventos.completos = eventos_completos,
                            tests.estacionariedad = config$params$tests.estacionariedad,
                            umbral.p.valor = config$params$umbral.p.valor)

# Transformar resultados a un objeto de tipo tibble
estacionariedad <- estacionariedad %>% purrr::map_dfr(~.x)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.errors <- task$getErrors()
if (length(task.errors) > 0) {
  for (error.obj in task.errors) {
    id_column <- IdentificarIdColumn(series_a_computar[1,])
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}, ",
                           "variable=\"{error.obj$input.value[['variable']]}\", ",
                           "tipo_serie=\"{error.obj$input.value[['tipo_serie']]}\", ",
                           "n_serie={error.obj$input.value[['n_serie']]}, ",
                           "realizacion={error.obj$input.value[['realizacion']]})",
                           ": {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
} else {
  # Guardar resultados en un archivo fácil de compartir
  script$info(glue::glue("Guardando estacionariedad de las series en el archivo {results_filename}"))
  saveRDS(estacionariedad, results_filename)
}

# Consolidar resultados. La estacionariedad es la única métrica en la que se consideran
# separadamente cada una de las realizaciones del generador, por lo tanto, esto debe
# consolidarse de manera a reagrupar las realizaciones. EL criterio, por el momento,
# es tomar la estacionaridad de la mayoría de las realizaciones.
estacionariedad.consolidada <- estacionariedad %>%
  dplyr::group_by(!!rlang::sym(id_column), nombre, variable, tipo_serie, n_serie) %>%
  dplyr::summarise(es_estacionaria = sum(es_estacionaria) >= round(n()/2))

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 9. Determinar la dependencia entre las series asociadas a las
# --- variables que conforman cada una de las cópulas.
# --- Clave primaria: ubicación, serie, variable_x, variable_y ----
# -----------------------------------------------------------------------------#

# Definir el objetos sobre los cuales iterar
input_copulas <- variables_copulas %>%
  tidyr::crossing(ubicaciones, n_serie = 0:config$params$n.series.perturbadas) %>%
  dplyr::mutate(tipo_serie = ifelse(n_serie == 0, 'observada', 'perturbada')) %>%
  dplyr::select(!!id_column, nombre, tipo_serie, n_serie, variable_x, variable_y) %>%
  dplyr::inner_join(estacionariedad.consolidada %>% dplyr::rename(x_es_estacionaria = es_estacionaria), 
                    by = c(id_column, 'nombre', 'variable_x' = 'variable', 'tipo_serie', 'n_serie')) %>%
  dplyr::inner_join(estacionariedad.consolidada %>% dplyr::rename(y_es_estacionaria = es_estacionaria), 
                    by = c(id_column, 'nombre', 'variable_y' = 'variable', 'tipo_serie', 'n_serie'))


# Definir el nombre de la función a ser paralelizada
function_name <- "DeterminarDependencia"

# Definir nombre de archivos .log y .out de corridas anteriores
task_logfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.log")
task_outfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.out")

# Borrar archivos .log y .out de corridas anteriores
if (file.exists(task_logfile))
  file.remove(task_logfile)
if (file.exists(task_outfile))
  file.remove(task_outfile)

# Definir nombre del archivo donde se van a guardar los resultados
results_filename <- glue::glue("{config$dir$data}/{config$files$copulas$dependencia}")

# Borrar archivo de resultado de corridas anteriores
if (file.exists(results_filename))
  file.remove(results_filename)

# Crear tarea distribuida y ejecutarla
task <- Task$new(parent.script = script,
                 func.name = function_name,
                 packages = list.of.packages)

# Informar inicio de ejecución 
script$info("Determinando la dependencia entre las series perturbadas que conforman las cópulas")
# Ejecutar tarea distribuida
dependencia <- task$run(number.of.processes = config$max.procesos,
                        input.values = input_copulas, 
                        eventos.completos = eventos_completos,
                        tests.dependencia = config$params$tests.dependencia,
                        umbral.p.valor =  config$params$umbral.p.valor)

# Transformar resultados a un objeto de tipo tibble
dependencia <- dependencia %>% purrr::map_dfr(~.x)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.errors <- task$getErrors()
if (length(task.errors) > 0) {
  for (error.obj in task.errors) {
    id_column <- IdentificarIdColumn(input_copulas[1,])
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}, ",
                           "copula=\"{error.obj$input.value[['variable_x']]}-{error.obj$input.value[['variable_y']]}\" ",
                           "tipo_serie=\"{error.obj$input.value[['tipo_serie']]}\", ",
                           "n_serie={error.obj$input.value[['n_serie']]})",
                           ": {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
} else {
  # Guardar resultados en un archivo fácil de compartir
  script$info(glue::glue("Guardando dependencia de las series en el archivo {results_filename}"))
  saveRDS(dependencia, results_filename)
}

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 10. Ajustar cópulas usando cada una de las series perturbadas. El
# --- ajuste se realiza únicamente cuando las dos series perturbas son estacionarias 
# --- y, además, ambas series perturbadas son dependientes ----
# -----------------------------------------------------------------------------#

# Definir el objetos sobre los cuales iterar
copulas_a_ajustar <- dependencia %>%
  dplyr::select(-detalles_dependencia, -segundos_calc_dependencia) %>%
  tidyr::crossing(configuracion.ajuste.copula)


# Definir el nombre de la función a ser paralelizada
function_name <- "AjustarCopulas"

# Definir nombre de archivos .log y .out de corridas anteriores
task_logfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.log")
task_outfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.out")

# Borrar archivos .log y .out de corridas anteriores
if (file.exists(task_logfile))
  file.remove(task_logfile)
if (file.exists(task_outfile))
  file.remove(task_outfile)

# Definir nombre del archivo donde se van a guardar los resultados
results_filename <- glue::glue("{config$dir$data}/{config$files$copulas$ajustes_multivariados}")

# Borrar archivo de resultado de corridas anteriores
if (file.exists(results_filename))
  file.remove(results_filename)

# Crear tarea distribuida y ejecutarla
task <- Task$new(parent.script = script,
                 func.name = function_name,
                 packages = list.of.packages)

# Informar inicio de ejecución 
script$info("Ajustar las cópulas usando las series perturbadas estacionarias y dependientes")
# Ejecutar tarea distribuida
copulas.ajustadas <- task$run(number.of.processes = config$max.procesos,
                              input.values = copulas_a_ajustar,
                              eventos.completos = eventos_completos,
                              umbral.p.valor = config$params$umbral.p.valor)

# Transformar resultados a un objeto de tipo tibble
copulas.ajustadas <- copulas.ajustadas %>% purrr::map_dfr(~.x)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.errors <- task$getErrors()
if (length(task.errors) > 0) {
  for (error.obj in task.errors) {
    id_column <- IdentificarIdColumn(copulas_a_ajustar[1,])
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}, ",
                           "copula=\"{error.obj$input.value[['variable_x']]}-{error.obj$input.value[['variable_y']]}\")",
                           "tipo_serie=\"{error.obj$input.value[['tipo_serie']]}\", ",
                           "n_serie={error.obj$input.value[['n_serie']]})",
                           ": {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
} else {
  # Guardar resultados en un archivo fácil de compartir
  script$info(glue::glue("Guardando los ajustes multivariados en el archivo {results_filename}"))
  saveRDS(copulas.ajustadas, results_filename)
}

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 11. Determinar la familia que mejor ajusta cada una de las cópulas
# --- consideradas en cada una de las ubicaciones abordadas ---- 
# -----------------------------------------------------------------------------#

# Definir el objetos sobre los cuales iterar
copulas_x_ubicacion <- copulas.ajustadas %>% 
  dplyr::select(!!id_column, nombre, variable_x, variable_y) %>% 
  dplyr::distinct()


# Definir el nombre de la función a ser paralelizada
function_name <- "MejorAjusteMultivariadoUC"

# Definir nombre de archivos .log y .out de corridas anteriores
task_logfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.log")
task_outfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.out")

# Borrar archivos .log y .out de corridas anteriores
if (file.exists(task_logfile))
  file.remove(task_logfile)
if (file.exists(task_outfile))
  file.remove(task_outfile)

# Definir nombre del archivo donde se van a guardar los resultados
results_filename <- glue::glue("{config$dir$data}/{config$files$copulas$mejores_ajustes_multivariados}")

# Borrar archivo de resultado de corridas anteriores
if (file.exists(results_filename))
  file.remove(results_filename)

# Crear tarea distribuida y ejecutarla
task <- Task$new(parent.script = script,
                 func.name = function_name,
                 packages = list.of.packages)

# Informar inicio de ejecución 
script$info("Determinar la familia que mejor ajusta a cada una de las cópulas en cada ubicación")
# Ejecutar tarea distribuida
mejor.ajuste.copulas <- task$run(number.of.processes = config$max.procesos,
                                 input.values = copulas_x_ubicacion,
                                 copulas.ajustadas = copulas.ajustadas,
                                 umbral.p.valor = config$params$umbral.p.valor)

# Transformar resultados a un objeto de tipo tibble
mejor.ajuste.copulas <- mejor.ajuste.copulas %>% purrr::map_dfr(~.x)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.errors <- task$getErrors()
if (length(task.errors) > 0) {
  for (error.obj in task.errors) {
    id_column <- IdentificarIdColumn(copulas_x_ubicacion[1,])
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}, ",
                           "copula=\"{error.obj$input.value[['variable_x']]}-{error.obj$input.value[['variable_y']]}\")",
                           "familia=\"{error.obj$input.value[['familia']]}\", funcion_ajuste=\"{error.obj$input.value[['funcion_ajuste']]}\")",
                           ": {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
} else {
  # Guardar resultados en un archivo fácil de compartir
  script$info(glue::glue("Guardando mejores ajustes multivariados en el archivo {results_filename}"))
  saveRDS(mejor.ajuste.copulas, results_filename)
}

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 12. Finalmente, realizar el ajuste multivariado ---- 
# --- OBS: como el ajuste multivariado también debe realizarse sobre la serie 
# --- observada, se combina combina ésta con las series perturbadas en 
# --- un solo input, pero con n_serie_perturbada igual a cero!!
# -----------------------------------------------------------------------------#

# Definir el objetos sobre los cuales iterar
copulas_x_ubicacion_x_serie <- mejor.ajuste.copulas %>%
  dplyr::select(!!id_column, nombre, variable_x, variable_y) %>%
  tidyr::crossing(n_serie = 0:config$params$n.series.perturbadas) %>%
  dplyr::mutate(tipo_serie = ifelse(n_serie == 0, 'observada', 'perturbada')) %>%
  dplyr::select(!!id_column, nombre, variable_x, variable_y, tipo_serie, n_serie)


# Definir el nombre de la función a ser paralelizada
function_name <- "AplicarMejorAjusteACopulas"

# Definir nombre de archivos .log y .out de corridas anteriores
task_logfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.log")
task_outfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.out")

# Borrar archivos .log y .out de corridas anteriores
if (file.exists(task_logfile))
  file.remove(task_logfile)
if (file.exists(task_outfile))
  file.remove(task_outfile)

# Definir nombre del archivo donde se van a guardar los resultados
results_filename <- glue::glue("{config$dir$data}/{config$files$copulas$resultado_final}")

# Borrar archivo de resultado de corridas anteriores
if (file.exists(results_filename))
  file.remove(results_filename)

# Crear tarea distribuida y ejecutarla
task <- Task$new(parent.script = script,
                 func.name = function_name,
                 packages = list.of.packages)

# Informar inicio de ejecución 
script$info("Aplicar el mejor ajuste multivariado a cada una de las cópulas en cada ubicación")
# Ejecutar tarea distribuida
ajuste.multivariado.final <- task$run(number.of.processes = config$max.procesos,
                                      input.values = copulas_x_ubicacion_x_serie,
                                      copulas.ajustadas = copulas.ajustadas,
                                      mejor.ajuste.univariado = mejor.ajuste.univariado.x.ubic.var,
                                      mejor.ajuste.multivariado = mejor.ajuste.copulas,
                                      eventos.completos = eventos_completos)

# Transformar resultados a un objeto de tipo tibble
ajuste.multivariado.final <- ajuste.multivariado.final %>% purrr::map_dfr(~.x)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.errors <- task$getErrors()
if (length(task.errors) > 0) {
  for (error.obj in task.errors) {
    id_column <- IdentificarIdColumn(copulas_x_ubicacion_x_serie[1,])
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}, ",
                           "copula=\"{error.obj$input.value[['variable_x']]}-{error.obj$input.value[['variable_y']]}\", ",
                           "tipo_serie=\"{error.obj$input.value[['tipo_serie']]}\", ",
                           "n_serie={error.obj$input.value[['n_serie']]})",
                           ": {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
} else {
  # Guardar resultados en un archivo fácil de compartir
  script$info(glue::glue("Guardando ajuste multivariado final en el archivo {results_filename}"))
  saveRDS(ajuste.multivariado.final, results_filename)
}

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO XX. Finalizar script ----
# -----------------------------------------------------------------------------#

# a) Finalizar script
script$stop()

# b) Crear archivo .info
info_filename <- glue::glue("{config$dir$data}/{config$files$copulas$info_corrida}")
if (file.exists(info_filename))
  file.remove(info_filename)
file.copy(from = script_logfile, to = info_filename)

# ------------------------------------------------------------------------------