# -----------------------------------------------------------------------------#
# ---- Funciona de ajuste de copulas ----

# --- Definicion de funciones para copulas Arquimedianas ----

# Definicion de funcion para el ajuste de una copula de la familia Gumbel
AjustarCopulaGumbel <- function(x, y) {
  
  # Inicializar objeto a devolver
  copula <- list(familia = 'gumbel', copula = NA, 
    dependencia = list(tau.kendall = NA),
    parametro = list(theta = NA),
    bondad.ajuste = list(varianza = NA, loglike = NA))
  
  # Creacion de matriz de eventos en valores absolutos para evitar dependencias negativas
  matriz.eventos <- as.matrix(cbind(abs(x), abs(y)))
  # Convertir la matriz de eventos en pseudoobservaciones
  matriz.eventos.pobs <- copula::pobs(matriz.eventos, ties.method = 'random')
  
  # La matriz debe tener al menos dos variables
  base::stopifnot(ncol(matriz.eventos.pobs) >= 2)
  
  # Estimacion del tau de Kendall
  tau <- cor(matriz.eventos.pobs[,1], matriz.eventos.pobs[,2], method='kendall')
  
  # Detener la estimacion si tau es negativo porque está fuera del dominio de la copula
  base::stopifnot(tau > 0)
  
  # Inicializar parametros para la familia Gumbel a partir del tau de Kendall
  theta.0 <- copula::initOpt("Gumbel",  interval = FALSE,  u = matriz.eventos.pobs, method="tau.mean")
  
  # Crear objeto copula de la familia Gumbel
  # Las dimensiones de la copula son iguales a la cantidad de variables
  gumbel.copula.inicial <- copula::gumbelCopula(param = theta.0,
    dim = ncol(matriz.eventos.pobs))  
  
  # Ajuste de copula Gumbel
  ajuste.gumbel.copula.mpl <- copula::fitCopula(copula = gumbel.copula.inicial,
    data = matriz.eventos.pobs, start = theta.0, estimate.variance = T,
    lower = 1.0001, upper = Inf)
  
  # Comprobar que el parametro ajustado se encuentra dentro del rango válido para la familia
  CI <- stats::confint(ajuste.gumbel.copula.mpl, quietly = TRUE)
  stopifnot(CI[1] <= ajuste.gumbel.copula.mpl@estimate, ajuste.gumbel.copula.mpl@estimate <= CI[2])
  
  # Devolver objeto resultado
  copula$copula <- gumbelCopula(param = ajuste.gumbel.copula.mpl@estimate, dim = ncol(matriz.eventos.pobs))
  copula$dependencia$tau.kendall <- tau
  copula$parametro$theta <- unname(ajuste.gumbel.copula.mpl@estimate)
  copula$bondad.ajuste$varianza <- as.vector(unname(ajuste.gumbel.copula.mpl@var.est))
  copula$bondad.ajuste$loglike <- unname(ajuste.gumbel.copula.mpl@loglik)
  
  return(copula)
}

# Definicion de funcion para el ajuste de una copula de la familia Frank
AjustarCopulaFrank <- function(x, y) {
  
  # Inicializar objeto a devolver
  copula <- list(familia = 'frank', copula = NA, 
    dependencia = list(tau.kendall = NA),
    parametro = list(theta = NA),
    bondad.ajuste = list(varianza = NA, loglike = NA))
  
  # Creacion de matriz de eventos en valores absolutos para evitar dependencias negativas
  matriz.eventos <- as.matrix(cbind(abs(x), abs(y)))
  # Convertir la matriz de eventos en pseudoobservaciones
  matriz.eventos.pobs <- copula::pobs(matriz.eventos, ties.method = 'random')
  
  # La matriz debe tener al menos dos variables
  base::stopifnot(ncol(matriz.eventos.pobs) >= 2)
  
  # Estimacion del tau de Kendall
  tau <- cor(matriz.eventos.pobs[,1], matriz.eventos.pobs[,2], method='kendall')
  
  # Detener la estimacion si tau es negativo porque está fuera del dominio de la copula
  base::stopifnot(tau > 0)
  
  # Inicializar parametros para la familia Gumbel a partir del tau de Kendall
  theta.0 <- copula::initOpt("Frank",  interval = FALSE,  u = matriz.eventos.pobs, method="tau.mean")
  
  # Crear objeto copula de la familia Gumbel
  # Las dimensiones de la copula son iguales a la cantidad de variables
  frank.copula.inicial <- copula::frankCopula(param = theta.0,
    dim = ncol(matriz.eventos.pobs))  
  
  # Ajuste de copula 
  ajuste.frank.copula.mpl <- copula::fitCopula(copula = frank.copula.inicial,
    data = matriz.eventos.pobs, start = theta.0, estimate.variance = T,
    lower = -Inf, upper = Inf)
  
  # Comprobar que el parametro ajustado se encuentra dentro del rango válido para la familia
  CI <- stats::confint(ajuste.frank.copula.mpl, quietly = TRUE)
  stopifnot(CI[1] <= ajuste.frank.copula.mpl@estimate, ajuste.frank.copula.mpl@estimate <= CI[2])
  
  # Devolver objeto resultado
  copula$copula <- copula::frankCopula(param = ajuste.frank.copula.mpl@estimate, dim = ncol(matriz.eventos.pobs))
  copula$dependencia$tau.kendall <- tau
  copula$parametro$theta <- unname(ajuste.frank.copula.mpl@estimate)
  copula$bondad.ajuste$varianza <- as.vector(unname(ajuste.frank.copula.mpl@var.est))
  copula$bondad.ajuste$loglike <- unname(ajuste.frank.copula.mpl@loglik)
  
  return(copula)
}

# Definicion de funcion para el ajuste de una copula de la familia AMH
AjustarCopulaAMH <- function(x, y) {
  
  # Inicializar objeto a devolver
  copula <- list(familia = 'amh', copula = NA, 
    dependencia = list(tau.kendall = NA),
    parametro = list(theta = NA),
    bondad.ajuste = list(varianza = NA, loglike = NA))
  
  # Creacion de matriz de eventos en valores absolutos para evitar dependencias negativas
  matriz.eventos <- as.matrix(cbind(abs(x), abs(y)))
  # Convertir la matriz de eventos en pseudoobservaciones
  matriz.eventos.pobs <- copula::pobs(matriz.eventos, ties.method = 'random')
  
  # La matriz debe tener al menos dos variables
  base::stopifnot(ncol(matriz.eventos.pobs) >= 2)
  
  # Estimacion del tau de Kendall
  tau <- cor(matriz.eventos.pobs[,1], matriz.eventos.pobs[,2], method='kendall')
  
  # Detener la estimacion si tau es negativo porque está fuera del dominio de la copula
  base::stopifnot(tau > 0)
  
  # Inicializar parametros para la familia Gumbel a partir del tau de Kendall
  theta.0 <- copula::initOpt("AMH",  interval = FALSE,  u = matriz.eventos.pobs, method="tau.mean")
  
  # Crear objeto copula de la familia Gumbel
  # Las dimensiones de la copula son iguales a la cantidad de variables
  amh.copula.inicial <- copula::amhCopula(param = theta.0,
    dim = ncol(matriz.eventos.pobs))  
  
  # Ajuste de copula Gumbel
  ajuste.amh.copula.mpl <- copula::fitCopula(copula = amh.copula.inicial,
    data = matriz.eventos.pobs, start = theta.0)
  
  # Comprobar que el parametro ajustado se encuentra dentro del rango válido para la familia
  CI <- stats::confint(ajuste.amh.copula.mpl, quietly = TRUE)
  stopifnot(CI[1] <= ajuste.amh.copula.mpl@estimate, ajuste.amh.copula.mpl@estimate <= CI[2])
  
  # Devolver objeto resultado
  copula$copula <- copula::amhCopula(param = ajuste.amh.copula.mpl@estimate, dim = ncol(matriz.eventos.pobs))
  copula$dependencia$tau.kendall <- tau
  copula$parametro$theta <- unname(ajuste.amh.copula.mpl@estimate)
  copula$bondad.ajuste$varianza <- as.vector(unname(ajuste.amh.copula.mpl@var.est))
  copula$bondad.ajuste$loglike <- unname(ajuste.amh.copula.mpl@loglik)
  
  return(copula)
}

# Definicion de funcion para el ajuste de una copula de la familia Joe
AjustarCopulaJoe <- function(x, y) {
  
  # Inicializar objeto a devolver
  copula <- list(familia = 'joe', copula = NA, 
    dependencia = list(tau.kendall = NA),
    parametro = list(theta = NA),
    bondad.ajuste = list(varianza = NA, loglike = NA))
  
  # Creacion de matriz de eventos en valores absolutos para evitar dependencias negativas
  matriz.eventos <- as.matrix(cbind(abs(x), abs(y)))
  # Convertir la matriz de eventos en pseudoobservaciones
  matriz.eventos.pobs <- copula::pobs(matriz.eventos, ties.method = 'random')
  
  # La matriz debe tener al menos dos variables
  base::stopifnot(ncol(matriz.eventos.pobs) >= 2)
  
  # Estimacion del tau de Kendall
  tau <- cor(matriz.eventos.pobs[,1], matriz.eventos.pobs[,2], method = 'kendall')
  
  # Detener la estimacion si tau es negativo porque está fuera del dominio de la copula
  base::stopifnot(tau > 0)
  
  # Inicializar parametros para la familia Gumbel a partir del tau de Kendall
  theta.0 <- copula::initOpt("Joe",  interval = FALSE,  u = matriz.eventos.pobs, method="tau.mean")
  
  # Crear objeto copula de la familia Gumbel
  # Las dimensiones de la copula son iguales a la cantidad de variables
  joe.copula.inicial <- copula::joeCopula(param = theta.0,
    dim = ncol(matriz.eventos.pobs))  
  
  # Ajuste de copula Gumbel
  ajuste.joe.copula.mpl <- copula::fitCopula(copula = joe.copula.inicial,
    data = matriz.eventos.pobs, start = theta.0)
  
  # Comprobar que el parametro ajustado se encuentra dentro del rango válido para la familia
  CI <- stats::confint(ajuste.joe.copula.mpl, quietly = TRUE)
  stopifnot(CI[1] <= ajuste.joe.copula.mpl@estimate, ajuste.joe.copula.mpl@estimate <= CI[2])
  
  # Devolver objeto resultado
  copula$copula <- copula::joeCopula(param = ajuste.joe.copula.mpl@estimate, dim = ncol(matriz.eventos.pobs))
  copula$dependencia$tau.kendall <- tau
  copula$parametro$theta <- unname(ajuste.joe.copula.mpl@estimate)
  copula$bondad.ajuste$varianza <- as.vector(unname(ajuste.joe.copula.mpl@var.est))
  copula$bondad.ajuste$loglike <- unname(ajuste.joe.copula.mpl@loglik)
  
  return(copula)
}

# Definicion de funcion para el ajuste de una copula de la familia Clayton
AjustarCopulaClayton <- function(x, y) {
  
  # Inicializar objeto a devolver
  copula <- list(familia = 'clayton', copula = NA, 
    dependencia = list(tau.kendall = NA),
    parametro = list(theta = NA),
    bondad.ajuste = list(varianza = NA, loglike = NA))
  
  # Creacion de matriz de eventos en valores absolutos para evitar dependencias negativas
  matriz.eventos <- as.matrix(cbind(abs(x), abs(y)))
  # Convertir la matriz de eventos en pseudoobservaciones
  matriz.eventos.pobs <- copula::pobs(matriz.eventos, ties.method = 'random')
  
  # La matriz debe tener al menos dos variables
  base::stopifnot(ncol(matriz.eventos.pobs) >= 2)
  
  # Estimacion del tau de Kendall
  tau <- cor(matriz.eventos.pobs[,1], matriz.eventos.pobs[,2], method = 'kendall')
  
  # Detener la estimacion si tau es negativo porque está fuera del dominio de la copula
  base::stopifnot(tau > 0)
  
  # Inicializar parametros para la familia Gumbel a partir del tau de Kendall
  theta.0 <- copula::initOpt("Clayton",  interval = FALSE,  u = matriz.eventos.pobs, method="tau.mean")
  
  # Crear objeto copula de la familia Gumbel
  # Las dimensiones de la copula son iguales a la cantidad de variables
  clayton.copula.inicial <- copula::claytonCopula(param = theta.0,
    dim = ncol(matriz.eventos.pobs))  
  
  # Ajuste de copula Gumbel
  ajuste.clayton.copula.mpl <- copula::fitCopula(copula = clayton.copula.inicial,
    data = matriz.eventos.pobs, start = theta.0)
  
  # Comprobar que el parametro ajustado se encuentra dentro del rango válido para la familia
  CI <- stats::confint(ajuste.clayton.copula.mpl, quietly = TRUE)
  stopifnot(CI[1] <= ajuste.clayton.copula.mpl@estimate, ajuste.clayton.copula.mpl@estimate <= CI[2])
  
  # Devolver objeto resultado
  copula$copula <- copula::claytonCopula(param = ajuste.clayton.copula.mpl@estimate, dim = ncol(matriz.eventos.pobs))
  copula$dependencia$tau.kendall <- tau
  copula$parametro$theta <- unname(ajuste.clayton.copula.mpl@estimate)
  copula$bondad.ajuste$varianza <- as.vector(unname(ajuste.clayton.copula.mpl@var.est))
  copula$bondad.ajuste$loglike <- unname(ajuste.clayton.copula.mpl@loglik)
  
  return(copula)
}

# ------------------------------------------------------------------------------

# --- Definicion de funciones para copulas Elípticas ----

# Definicion de funcion para el ajuste de una copula de la familia Normal
AjustarCopulaNormal <- function(x, y) {
  
  # Inicializar objeto a devolver
  copula <- list(familia = 'normal', copula = NA, 
    dependencia = list(tau.kendall = NA),
    parametro = list(rho = NA),
    bondad.ajuste = list(varianza = NA, loglike = NA))
  
  # Creacion de matriz de eventos en valores absolutos para evitar dependencias negativas
  matriz.eventos <- as.matrix(cbind(abs(x), abs(y)))
  # Convertir la matriz de eventos en pseudoobservaciones
  matriz.eventos.pobs <- copula::pobs(matriz.eventos, ties.method = 'random')
  
  # La matriz debe tener al menos dos variables
  base::stopifnot(ncol(matriz.eventos.pobs) >= 2)
  
  # Estimacion del tau de Kendall
  tau <- cor(matriz.eventos.pobs[,1], matriz.eventos.pobs[,2], method='kendall')
  
  # Detener la estimacion si tau es negativo porque está fuera del dominio de la copula
  base::stopifnot(tau > 0)
  
  # Crear objeto copula de la familia Gumbel
  # Las dimensiones de la copula son iguales a la cantidad de variables
  normal.copula.inicial <- copula::normalCopula(param = NA_real_,
    dim = ncol(matriz.eventos.pobs), dispstr = 'un')  
  
  # Ajuste de copula Gumbel
  ajuste.normal.copula.mpl <- copula::fitCopula(copula = normal.copula.inicial,
    data = matriz.eventos.pobs, estimate.variance = T)
  
  # Comprobar que el parametro ajustado se encuentra dentro del rango válido para la familia
  CI <- stats::confint(ajuste.normal.copula.mpl, quietly = TRUE)
  stopifnot(CI[1] <= ajuste.normal.copula.mpl@estimate, ajuste.normal.copula.mpl@estimate <= CI[2])
  
  # Extraer parametro inicial de la copula Gumbel
  theta <- unname(coef(ajuste.normal.copula.mpl))
  
  # Devolver objeto resultado
  copula$copula <- copula::normalCopula(param = ajuste.normal.copula.mpl@estimate, dim = ncol(matriz.eventos.pobs))
  copula$dependencia$tau.kendall <- tau
  copula$parametro$rho <- unname(ajuste.normal.copula.mpl@estimate)
  copula$bondad.ajuste$varianza <- as.vector(unname(ajuste.normal.copula.mpl@var.est))
  copula$bondad.ajuste$loglike <- unname(ajuste.normal.copula.mpl@loglik)
  
  return(copula)
}

# Definicion de funcion para el ajuste de una copula de la familia t-Student
AjustarCopulatstudent <- function(x, y) {
  
  # Inicializar objeto a devolver
  copula <- list(familia = 't', copula = NA, 
    dependencia = list(tau.kendall = NA),
    parametro = list(rho = NA),
    bondad.ajuste = list(varianza = NA, loglike = NA))
  
  # Creacion de matriz de eventos en valores absolutos para evitar dependencias negativas
  matriz.eventos <- as.matrix(cbind(abs(x), abs(y)))
  # Convertir la matriz de eventos en pseudoobservaciones
  matriz.eventos.pobs <- copula::pobs(matriz.eventos, ties.method = 'random')
  
  # La matriz debe tener al menos dos variables
  base::stopifnot(ncol(matriz.eventos.pobs) >= 2)
  
  # Estimacion del tau de Kendall
  tau <- cor(matriz.eventos.pobs[,1], matriz.eventos.pobs[,2], method = 'kendall')
  
  # Detener la estimacion si tau es negativo porque está fuera del dominio de la copula
  base::stopifnot(tau > 0)
  
  # Crear objeto copula de la familia Gumbel
  # Las dimensiones de la copula son iguales a la cantidad de variables
  t.copula.inicial <- copula::tCopula(param = NA_real_,
    dim = ncol(matriz.eventos.pobs), dispstr = 'un', df = 4, df.fixed = T)  
  
  # Ajuste de copula Gumbel
  ajuste.t.copula.mpl <- copula::fitCopula(copula = t.copula.inicial,
    data = matriz.eventos.pobs,  method = 'mpl')
  
  # Devolver objeto resultado
  copula$copula <- copula::tCopula(param = ajuste.t.copula.mpl@estimate, dim = ncol(matriz.eventos.pobs), df.fixed = T)
  copula$dependencia$tau.kendall <- tau
  copula$parametro$rho <- unname(ajuste.t.copula.mpl@estimate)
  copula$bondad.ajuste$varianza <- as.vector(unname(ajuste.t.copula.mpl@var.est))
  copula$bondad.ajuste$loglike <- unname(ajuste.t.copula.mpl@loglik)
  
  return(copula)
  
}

# ------------------------------------------------------------------------------



