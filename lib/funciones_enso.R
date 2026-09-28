# Definicion de funcion para leer el archivo de eventos ENSO (RONI trimestral
# solapado) y calcular, para cada temporada, el mes/anio calendario real de
# inicio y fin de la ventana de 3 meses que representa ----
LeerEventosEnso <- function(archivo) {
  # Offset de anio (respecto de la columna Year) que corresponde al mes de
  # inicio/fin de cada temporada solapada de 3 meses. Solo DJF y NDJ cruzan
  # el limite de anio calendario (convencion estandar CPC/ONI): DJF(Year) es
  # Dic(Year-1)-Ene(Year)-Feb(Year); NDJ(Year) es Nov(Year)-Dic(Year)-Ene(Year+1).
  # El resto de las temporadas (JFM...OND) no tienen offset.
  offsets_temporada <- tibble::tribble(
    ~Season, ~offset_inicio, ~offset_fin,
    "DJF",   -1L,             0L,
    "NDJ",    0L,             1L
  )

  eventos_enso <- data.table::fread(archivo, sep = "\t", encoding = "UTF-8") %>%
    dplyr::as_tibble() %>%
    dplyr::left_join(offsets_temporada, by = "Season") %>%
    dplyr::mutate(offset_inicio = dplyr::coalesce(offset_inicio, 0L),
                  offset_fin    = dplyr::coalesce(offset_fin, 0L),
                  ym_inicio     = (Year + offset_inicio) * 12L + Start_Month,
                  ym_fin        = (Year + offset_fin)    * 12L + End_Month,
                  orden         = dplyr::row_number())

  # Devolver objeto resultado
  return(eventos_enso)
}
#-------------------------------------------------------------------------------

# Definicion de funcion para aplicar la regla de desambiguacion Nino/Nina/Neutro
# a las (hasta 3) clasificaciones ENSO de temporadas solapadas que coinciden
# con un mismo mes calendario. tipos.evento debe venir ordenado cronologicamente ----
ClasificarGrupoEnso <- function(tipos.evento) {
  es.nino <- startsWith(tipos.evento, "Niño")
  es.nina <- startsWith(tipos.evento, "Niña")

  if (any(es.nino) && ! any(es.nina)) {
    resultado <- tipos.evento[es.nino][1]
  } else if (any(es.nina) && ! any(es.nino)) {
    resultado <- tipos.evento[es.nina][1]
  } else {
    resultado <- "Neutro"
  }

  # Devolver objeto resultado
  return(resultado)
}
#-------------------------------------------------------------------------------

# Definicion de funcion para, a partir del archivo de eventos ENSO, construir
# una tabla que indica para cada mes calendario cubierto por la serie (Year*12+mes)
# la clasificacion ENSO resultante de combinar las temporadas solapadas que lo incluyen ----
ClasificarMesesEnso <- function(eventos_enso) {
  meses.expandidos <- purrr::pmap_dfr(
    eventos_enso,
    function(...) {
      fila <- list(...)
      tibble::tibble(ym = fila$ym_inicio:fila$ym_fin,
                      tipo_evento = fila$Tipo_evento,
                      orden = fila$orden)
    }
  )

  clasificacion_enso <- meses.expandidos %>%
    dplyr::arrange(ym, orden) %>%
    dplyr::group_by(ym) %>%
    dplyr::summarise(evento_enso = ClasificarGrupoEnso(tipo_evento), .groups = "drop")

  # Devolver objeto resultado
  return(clasificacion_enso)
}
#-------------------------------------------------------------------------------

# Definicion de funcion para asignar a cada fecha de inicio de evento (ola de
# frio/calor, periodo calido/frio) la clasificacion ENSO correspondiente al
# periodo (anio, mes) en el que comenzo. Devuelve NA para fechas fuera del
# rango cubierto por el archivo de eventos ENSO (0 coincidencias) ----
AsignarEventoEnso <- function(fecha_inicio, clasificacion_enso) {
  ym <- lubridate::year(fecha_inicio) * 12L + lubridate::month(fecha_inicio)
  indices <- match(ym, clasificacion_enso$ym)

  # Devolver objeto resultado
  return(clasificacion_enso$evento_enso[indices])
}
#-------------------------------------------------------------------------------
