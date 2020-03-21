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
  .f = function(..., nombres_columnas) { tibble::tibble(...) %>% setNames(nombres_columnas) },
  nombres_columnas = config$params$configuracion.ajuste.univariado[[1]]
)
# Tomar solo configuraciones con uso_sequias igual a True
configuracion.ajuste.univariado <- configuracion.ajuste.univariado.completo %>% dplyr::filter(uso_sequias) %>%
  dplyr::select(distribucion, funcion_ajuste_lmomentos, funcion_ajuste_maxima_verosimilitud)

# Ajuste multivariado 
# Leer archivo de configuracion de funciones de ajuste multivariado
configuracion.ajuste.copula.completo <- purrr::pmap_dfr(
  .l = utils::tail(config$params$configuracion.ajuste.copula, -1) %>% purrr::transpose(), 
  .f = function(..., nombres_columnas) { tibble::tibble(...) %>% setNames(nombres_columnas) },
  nombres_columnas = config$params$configuracion.ajuste.copula[[1]]
)
# Tomar solo configuraciones con uso_sequias igual a True
configuracion.ajuste.copula <- configuracion.ajuste.copula.completo %>% dplyr::filter(uso_sequias) %>%
  dplyr::select(familia, funcion_ajuste)


## PARÁMETROS RESTANTES

# Obtener los eventos a ser analizados!
eventos <- feather::read_feather(glue::glue("{config$dir$data}/{config$files$eventos.identificados}")) 
id_column <- IdentificarIdColumn(eventos)

# Obtener ubicaciones a ser consideradas (deben existir en eventos!)
ubicaciones <- purrr::pmap_dfr(
  .l = config$params$ubicaciones %>% purrr::transpose(), 
  .f = function(...) { tibble::tibble(...) %>% setNames(c(id_column, 'nombre')) }
)
# Controlar que las estaciones seleccionadas estén entre aquellas 
# para las que se identificaron eventos en los pasos previos!
ids_in_params <- sapply(config$params$ubicaciones, function(u) {u$uid})
ids_in_events <- eventos %>% dplyr::pull(id_column) %>% as.character()
if (!all(ids_in_params %in% ids_in_events))
  stop("No hay datos, en eventos, para todas las ubicaciones en parametros_copulas.yml")

# Obtener variables a ser consideradas (deben existir en eventos!)
variables <- purrr::pmap_dfr(
  .l = config$params$variables_copulas %>% purrr::transpose(), 
  .f = function(...) { tibble::tibble(...) %>% setNames(c('variable_x', 'variable_y')) }
)
# Controlar que las variables seleccionadas sean columnas de eventos!
if (!all(union(variables$variable_x, variables$variable_y) %in% colnames(eventos)))
  stop("No hay datos, en eventos, para todas las variables en parametros_copulas.yml")

valores_minimos <- purrr::pmap_dfr(
  .l = config$params$valores_minimos %>% purrr::transpose(), 
  .f = function(...) { tibble::tibble(...) %>% setNames(c('variable', 'valor_minimo_deteccion')) }
)
# Controlar que hayan valores mínimos para todas las variables!
if (!all(union(variables$variable_x, variables$variable_y) %in% valores_minimos$variable))
  stop("No hay datos, en eventos, para todas las variables en parametros_copulas.yml")

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 5. Ajustar, por separado, cada una de las variables que conforman las
# --- cópulas, parámetro: variables_copulas, a las distribuciones señaladas en el  
# --- archivo de parámetros, parámetro: configuracion.ajuste.univariado. ----
# --- Clave primaria: ubicación, variable, distribución
# -----------------------------------------------------------------------------#


# Definir el objetos sobre los cuales iterar
ubicaciones_x_variables_x_distribuciones <- 
  tidyr::crossing(ubicaciones, variable = union(variables$variable_x, variables$variable_y),
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

# Crear tarea distribuida y ejecutarla
task.ajustar <- Task$new(parent.script = script,
                         func.name = function_name,
                         packages = list.of.packages)

# Informar inicio de ejecución 
script$info("Computando ajuste univariado para cada combinación de ubicación, variable, configuración")
# Ejecutar tarea distribuida
ajustes.x.ubic.var.conf <- task.ajustar$run(number.of.processes = config$max.procesos,
                                            input.values = ubicaciones_x_variables_x_distribuciones[c(1:4,40:43),],  
                                            config = config, eventos.completos = eventos)

# Transformar resultados a un objeto de tipo tibble
ajustes.x.ubic.var.conf.tibble <- ajustes.x.ubic.var.conf %>% purrr::map_dfr(~.x)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.ajustar.errors <- task.ajustar$getErrors()
if (length(task.ajustar.errors) > 0) {
  for (error.obj in task.ajustar.errors) {
    id_column <- IdentificarIdColumn(ubicaciones_x_variables_x_distribuciones[1,])
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}, ",
                           "variable=\"{error.obj$input.value[['variable']]}\", ",
                           "distribucion=\"{error.obj$input.value[['distribucion']]}\")",
                           ": {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
}

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 6. Determinar la distribución que mejor se ajusta a cada una de las
# --- variables en cada una de las ubicaciones. ----
# --- Clave primaria: ubicación, variable
# --- OBS: a partir de acá para cada variable, en cada ubicación, ya se sabe 
# --- cual es la distribución que mejor ajusta y el mejor método de ajuste
# -----------------------------------------------------------------------------#


# Definir el objetos sobre los cuales iterar
ubicaciones_x_variables <- 
  tidyr::crossing(ubicaciones, variable = union(variables$variable_x, variables$variable_y))


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

# Crear tarea distribuida y ejecutarla
task.mejor.ajuste <- Task$new(parent.script = script,
                              func.name = function_name,
                              packages = list.of.packages)

# Informar inicio de ejecución 
script$info("Determinando mejor ajuste univariado para cada combinación de ubicación, variable")
# Ejecutar tarea distribuida
mejor.ajuste.x.ubic.var <- task.mejor.ajuste$run(number.of.processes = config$max.procesos,
                                                 input.values = ubicaciones_x_variables,  
                                                 ajustes_univariados = ajustes.x.ubic.var.conf.tibble)

# Transformar resultados a un objeto de tipo tibble
mejor.ajuste.x.ubic.var.tibble <- mejor.ajuste.x.ubic.var %>% purrr::map_dfr(~.x)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.mejor.ajuste.errors <- task.mejor.ajuste$getErrors()
if (length(task.mejor.ajuste.errors) > 0) {
  for (error.obj in task.mejor.ajuste.errors) {
    id_column <- IdentificarIdColumn(ubicaciones_x_variables[1,])
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}, ",
                           "variable=\"{error.obj$input.value[['variable']]}\")",
                           ": {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
}

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 7. Generar series perturbadas para la distribución que mejor ajusta
# --- en cada una de las ubicaciones, y además, determinar los parámetros de
# --- ajuste de cada una de esas series perturbadas. ----
# --- Clave primaria: ubicación, variable, serie_perturbada
# --- OBS: la distribución ya no forma parte de la clave primaria porque cada
# --- par ubicación, variable, ya tiene definida la mejor distribución y el 
# --- mejor método de ajuste
# -----------------------------------------------------------------------------#


# Definir el objetos sobre los cuales iterar
mejor_ajuste_x_n_series_perturbadas <- mejor.ajuste.x.ubic.var.tibble %>% 
  dplyr::select(-parametros, -rmse, -ccc, -cuantiles) %>%
  dplyr::left_join(configuracion.ajuste.univariado, by = "distribucion") %>%
  dplyr::mutate(funcion_mejor_ajuste = ifelse(mejor_ajuste == 'lmomentos', 
                                              funcion_ajuste_lmomentos, 
                                              funcion_ajuste_maxima_verosimilitud)) %>%
  dplyr::select(-funcion_ajuste_lmomentos, -funcion_ajuste_maxima_verosimilitud) %>%
  tidyr::crossing(n_serie_perturbada = 1:config$params$n.series.perturbadas) %>%
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

# Crear tarea distribuida y ejecutarla
task.series.perturbadas <- Task$new(parent.script = script,
                                    func.name = function_name,
                                    packages = list.of.packages)

# Informar inicio de ejecución 
script$info("Generar series perturbadas y aplicarles el mejor ajuste")
# Ejecutar tarea distribuida
series.perturbadas.ajustadas <- task.series.perturbadas$run(number.of.processes = config$max.procesos,
                                                            input.values = mejor_ajuste_x_n_series_perturbadas,  
                                                            eventos.completos = eventos.completos)

# Transformar resultados a un objeto de tipo tibble
series.perturbadas.ajustadas.tibble <- series.perturbadas.ajustadas %>% purrr::map_dfr(~.x)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.series.perturbadas.errors <- task.series.perturbadas$getErrors()
if (length(task.series.perturbadas.errors) > 0) {
  for (error.obj in task.series.perturbadas.errors) {
    id_column <- IdentificarIdColumn(mejor_ajuste_x_n_series_perturbadas[1,])
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}, ",
                           "variable=\"{error.obj$input.value[['variable']]}\", ",
                           "distribucion=\"{error.obj$input.value[['distribucion']]}\", ",
                           "n_serie_perturbada={error.obj$input.value[['n_serie_perturbada']]})",
                           ": {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
}

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 8. Determinar la estacionariedad de la serie perturbada correspondiente
# --- a cada variable en cada ubicación ----
# --- Clave primaria: ubicación, variable, serie_perturbada
# -----------------------------------------------------------------------------#


# Definir el objetos sobre los cuales iterar
series_perturbadas_ajustadas <- series.perturbadas.ajustadas.tibble


# Definir el nombre de la función a ser paralelizada
function_name <- "DeterminarEstacionariedadDeSeriesPerturbadas"

# Definir nombre de archivos .log y .out de corridas anteriores
task_logfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.log")
task_outfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.out")

# Borrar archivos .log y .out de corridas anteriores
if (file.exists(task_logfile))
  file.remove(task_logfile)
if (file.exists(task_outfile))
  file.remove(task_outfile)

# Crear tarea distribuida y ejecutarla
task.estacionariedad <- Task$new(parent.script = script,
                                 func.name = function_name,
                                 packages = list.of.packages)

# Informar inicio de ejecución 
script$info("Generar series perturbadas y aplicarles el mejor ajuste")
# Ejecutar tarea distribuida
series.perturbadas.estacionariedad <- task.estacionariedad$run(number.of.processes = config$max.procesos,
                                                               input.values = series_perturbadas_ajustadas,  
                                                               umbral.p.valor = config$params$umbral.p.valor)

# Transformar resultados a un objeto de tipo tibble
series.perturbadas.estacionariedad.tibble <- series.perturbadas.estacionariedad %>% purrr::map_dfr(~.x)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.estacionariedad.errors <- task.estacionariedad$getErrors()
if (length(task.estacionariedad.errors) > 0) {
  for (error.obj in task.estacionariedad.errors) {
    id_column <- IdentificarIdColumn(mejor_ajuste_x_n_series_perturbadas[1,])
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}, ",
                           "variable=\"{error.obj$input.value[['variable']]}\", ",
                           "distribucion=\"{error.obj$input.value[['distribucion']]}\", ",
                           "n_serie_perturbada={error.obj$input.value[['n_serie_perturbada']]})",
                           ": {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
}

# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 9. Determinar la dependencia entre las series perturbadas que    ----
# --- conforman cada una de los pares de variables de las cópulas           ----
# -----------------------------------------------------------------------------#
# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- PASO 8. Finalizar script ----
# -----------------------------------------------------------------------------#

# a) Guardar resultados en un archivo fácil de compartir
results_filename <- glue::glue("{config$dir$data}/{config$files$copulas$resultados}")
script$info(glue::glue("Guardando resultados en el archivo {results_filename}"))
feather::write_feather(resultados.copulas.tibble, results_filename)

# b) Finalizar script
script$stop()

# c) Crear archivo .info
info_filename <- glue::glue("{config$dir$data}/{config$files$copulas$info_corrida}")
if (file.exists(info_filename))
  file.remove(info_filename)
file.copy(from = script_logfile, to = info_filename)

# ------------------------------------------------------------------------------