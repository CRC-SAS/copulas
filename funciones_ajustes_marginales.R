# ---------------------------------------------------------------------------- #
# Ajustar 4-Parameter Asymmetric Exponential Power ----
# ---------------------------------------------------------------------------- #

# Basado en Asquith, W. H. (2014). "Parameter estimation for the 4-parameter 
# Asymmetric Exponential Power distribution by the method of L-moments using R."
#C omputational Statistics & Data Analysis 71: 955-970.

# Ejemplo
# xi <- 0.002
# alpha <- 0.05
# kappa <- 1
# h <- 1.5
# set.seed(123)
# x <- lmomco::rlmomco(100, vec2par(c(xi, alpha, kappa, h), type = 'aep4'))
# y <- lmomco::rlmomco(100, vec2par(c(0.002150876, 0.04788925, 1.006429, 1.472859), type = 'aep4'))

# Ajuste por L-Momentos
AjustarLMomentosAsymetricExponentialPower <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'aep4',
    parametros = list(xi = NA, alpha = NA, kappa = NA, h = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (lmomco::are.lmom.valid(mom) && !is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      aep4.pars <- lmomco::lmom2par(mom, type='aep4')
      if (are.paraep4.valid(aep4.pars)) {
        ajuste$parametros$xi <- unname(aep4.pars$para[1])
        ajuste$parametros$alpha <- unname(aep4.pars$para[2])
        ajuste$parametros$kappa <- unname(aep4.pars$para[3])
        ajuste$parametros$h <- unname(aep4.pars$para[4])
      }
    } 
  }
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudAsymetricExponentialPower <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'aep4',
    parametros = list(xi = NA, alpha = NA, kappa = NA, h = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    aep4.pars.guess <- lmomco::paraep4(lmomco.fit, checklmom = TRUE)
    start.params <- list(xi = unname(aep4.pars.guess$para[1]), 
      alpha = unname(aep4.pars.guess$para[2]),
      kappa = unname(aep4.pars.guess$para[3]),
      h = unname(aep4.pars.guess$para[4]))
    
    # Realizar estimacion por metodo ML
    aep4.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'aep4',
      method = 'mle',
      keepdata = FALSE,
      start = start.params)
    ajuste <- list(xi = NA, alpha = NA, kappa = NA, h = NA)
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$xi <- unname(aep4.mle.fit$estimate[1])
      ajuste$alpha <- unname(aep4.mle.fit$estimate[2])
      ajuste$kappa <- unname(aep4.mle.fit$estimate[3])
      ajuste$h <- unname(aep4.mle.fit$estimate[4])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = aep4.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$xi <- unname(parametros.remuestreo[1])
      ajuste$alpha <- unname(parametros.remuestreo[2])
      ajuste$kappa <- unname(parametros.remuestreo[3])
      ajuste$h <- unname(parametros.remuestreo[4])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Cauchy ----
# ---------------------------------------------------------------------------- #

# Ejemplo
# xi <- 0 # location
# alpha <- 0.5 # scale
# set.seed(123)
# x <- lmomco::rlmomco(100, vec2par(c(xi, alpha), type = 'cau'))

# Ajuste por L-Momentos
AjustarLMomentosCauchy <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'cau',
    parametros = list(xi = NA, alpha = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::TLmoms(x, trim = 1, nm = 5)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (!is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      cau.pars <- lmomco::lmom2par(mom, type='cau')
      if (are.parcau.valid(cau.pars)) {
        ajuste$parametros$xi <- unname(cau.pars$para[1])
        ajuste$parametros$alpha <- unname(cau.pars$para[2])
      }
    } 
  }
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudCauchy <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'cau',
    parametros = list(xi = NA, alpha = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::TLmoms(x, trim = 1, nm = 5)
    cau.pars.guess <- lmomco::parcau(lmomco.fit, checklmom = TRUE)
    start.params <- list(xi = unname(cau.pars.guess$para[1]), 
      alpha = unname(cau.pars.guess$para[2]))
    
    # Realizar estimacion por metodo ML
    cau.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'cau',
      method = 'mle',
      keepdata = FALSE,
      start = start.params)
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$xi <- unname(cau.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(cau.mle.fit$estimate[2])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = cau.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$xi <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
    }
  }
  
  return (ajuste)
}
# ------------------------------------------------------------------------------ 

# ---------------------------------------------------------------------------- #
# Ajustar Exponencial ----
# ---------------------------------------------------------------------------- #


# alpha <- 0.5
# xi <- 10 
# x <- x <- lmomco::rlmomco(100, vec2par(c(xi, alpha), type = 'exp'))

# Ajuste por L-Momentos
AjustarLMomentosExponencial <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'exl',
    parametros = list(xi = NA, alpha = NA))
  if (length(x) >= min.cantidad.valores) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (lmomco::are.lmom.valid(mom) && !is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      exp.pars <- lmomco::parexp(mom, checklmom = TRUE)
      if (are.parexp.valid(exp.pars)) {
        ajuste$parametros$xi <- unname(exp.pars$para[1])
        ajuste$parametros$alpha <- unname(exp.pars$para[2])
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudExponencial <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Ajuste de la distribución Gamma por metodo de maxima verosimilitud
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'exl',
    parametros = list(xi = NA, alpha = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    exp.pars.guess <- lmomco::parexp(lmomco.fit, checklmom = TRUE)
    start.params <- list(xi = unname(exp.pars.guess$para[1]), 
      alpha = unname(exp.pars.guess$para[2]))
    
    # Realizar estimacion por metodo ML
    exp.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'exl',
      method = 'mle',
      keepdata = FALSE,
      start = start.params,
      lower = c(0, 0))
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$xi <- unname(exp.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(exp.mle.fit$estimate[2])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = exp.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$xi <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Gamma ----
# ---------------------------------------------------------------------------- #


# alpha = 1.2
# beta = 0.5
# x <- lmomco::rlmomco(100, vec2par(c(alpha, beta), type = 'gam'))

# Ajuste por L-Momentos
AjustarLMomentosGamma <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'gam',
    parametros = list(alpha = NA, beta = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (lmomco::are.lmom.valid(mom) && !is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      gam.pars <- lmomco::lmom2par(mom, type='gam')
      if (are.pargam.valid(gam.pars)) {
        ajuste$parametros$alpha <- unname(gam.pars$para[1])
        ajuste$parametros$beta <- unname(gam.pars$para[2])
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudGamma <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'gam',
    parametros = list(alpha = NA, beta = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    gam.pars.guess <- lmomco::pargam(lmomco.fit, checklmom = TRUE)
    start.params <- list(shape = unname(gam.pars.guess$para[1]), 
      scale = unname(gam.pars.guess$para[2]))
    
    # Realizar estimacion por metodo ML
    gam.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'gamma',
      method = 'mle',
      keepdata = FALSE,
      start = start.params,
      lower = c(0,0))
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$alpha <- unname(gam.mle.fit$estimate[1])
      ajuste$parametros$beta <- unname(gam.mle.fit$estimate[2])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = gam.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$alpha <- unname(parametros.remuestreo[1])
      ajuste$parametros$beta <- unname(parametros.remuestreo[2])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Generalized Exponential Poisson ----
# ---------------------------------------------------------------------------- #

# beta = 2700
# kappa = 3
# h = 5
# x <- lmomco::rlmomco(100, vec2par(c(beta, kappa, h), type = 'gep'))

# Ajuste por L-Momentos
AjustarLMomentosGeneralizedExponentialPoisson  <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'gep',
    parametros = list(beta = NA, kappa = NA, h = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (lmomco::are.lmom.valid(mom) && !is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      gep.pars <- lmomco::lmom2par(mom, type='gep')
      if (are.pargep.valid(gep.pars)) {
        ajuste$parametros$beta <- unname(gep.pars$para[1])
        ajuste$parametros$kappa <- unname(gep.pars$para[2])
        ajuste$parametros$h <- unname(gep.pars$para[3])
      }
    } 
  }
  
  return (ajuste)
}

AjustarMaximaVerosimilitudGeneralizedExponentialPoisson  <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(type = list(distribucion = 'gep'),
    parametros = list(beta = NA, kappa = NA, h = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    gep.pars.guess <- lmomco::pargep(lmomco.fit, checklmom = TRUE)
    start.params <- list(beta = unname(gep.pars.guess$para[1]), 
      kappa = unname(gep.pars.guess$para[2]),
      h = unname(gep.pars.guess$para[3]))
    
    # Realizar estimacion por metodo ML
    gep.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'gep',
      method = 'mle',
      keepdata = FALSE,
      start = start.params)
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$beta <- unname(gep.mle.fit$estimate[1])
      ajuste$parametros$kappa <- unname(gep.mle.fit$estimate[2])
      ajuste$parametros$h <- unname(gep.mle.fit$estimate[3])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = gep.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$beta <- unname(parametros.remuestreo[1])
      ajuste$parametros$kappa <- unname(parametros.remuestreo[2])
      ajuste$parametros$h <- unname(parametros.remuestreo[3])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Generalized Extreme Value ----
# ---------------------------------------------------------------------------- #

# xi = 0
# alpha = 1
# kappa = 0.5
# x <- lmomco::rlmomco(100, vec2par(c(xi, alpha, kappa), type = 'gev'))

# Ajuste por L-Momentos
AjustarLMomentosGeneralizedExtremeValue  <- function(x, min.cantidad.valores) {
 # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'gev',
    parametros = list(xi = NA, alpha = NA, kappa = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (lmomco::are.lmom.valid(mom) && !is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      gev.pars <- lmomco::lmom2par(mom, type='gev')
      if (are.pargev.valid(gev.pars)) {
        ajuste$parametros$xi <- unname(gev.pars$para[1])
        ajuste$parametros$alpha <- unname(gev.pars$para[2])
        ajuste$parametros$kappa <- unname(gev.pars$para[3])
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudGeneralizedExtremeValue  <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'gev',
    parametros = list(xi = NA, alpha = NA, kappa = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    gev.pars.guess <- lmomco::pargev(lmomco.fit, checklmom = TRUE)
    start.params <- list(xi = unname(gev.pars.guess$para[1]), 
      alpha = unname(gev.pars.guess$para[2]),
      kappa = unname(gev.pars.guess$para[3]))
    
    # Realizar estimacion por metodo ML
    gev.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'gev',
      method = 'mle',
      keepdata = FALSE,
      start = start.params)
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$xi <- unname(gev.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(gev.mle.fit$estimate[2])
      ajuste$parametros$kappa <- unname(gev.mle.fit$estimate[3])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = gev.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$xi <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
      ajuste$parametros$kappa <- unname(parametros.remuestreo[3])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Generalized Lambda ----
# ---------------------------------------------------------------------------- #

# xi = 0.0305000
# alpha = 0.7313684
# kappa = 0.0045810
# h <- 0.0102000
# x <- lmomco::rlmomco(100000, vec2par(c(xi, alpha, kappa, h), type = 'gld'))

# Ajuste por L-Momentos
AjustarLMomentosGeneralizedLambda  <- function(x, min.cantidad.valores) {
 # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'gld',
    parametros = list(xi = NA, alpha = NA, kappa = NA, h = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (lmomco::are.lmom.valid(mom) && !is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      gld.pars <- lmomco::lmom2par(mom, type='gld')
      if (are.pargld.valid(gld.pars)) {
        ajuste$parametros$xi <- unname(gld.pars$para[1])
        ajuste$parametros$alpha <- unname(gld.pars$para[2])
        ajuste$parametros$kappa <- unname(gld.pars$para[3])
        ajuste$parametros$h <- unname(gld.pars$para[4])
        
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudGeneralizedLambda  <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'gld',
    parametros = list(xi = NA, alpha = NA, kappa = NA, h = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    gld.pars.guess <- lmomco::pargld(lmomco.fit, checklmom = TRUE)
    start.params <- list(xi = unname(gld.pars.guess$para[1]), 
      alpha = unname(gld.pars.guess$para[2]),
      kappa = unname(gld.pars.guess$para[3]),
      h = unname(gld.pars.guess$para[4]))
    
    # Realizar estimacion por metodo ML
    gld.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'gld',
      method = 'mle',
      keepdata = FALSE,
      start = start.params)
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$xi <- unname(gld.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(gld.mle.fit$estimate[2])
      ajuste$parametros$kappa <- unname(gld.mle.fit$estimate[3])
      ajuste$parametros$h <- unname(gld.mle.fit$estimate[4])
      
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = gld.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$xi <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
      ajuste$parametros$kappa <- unname(parametros.remuestreo[3])
      ajuste$parametros$h <- unname(parametros.remuestreo[4])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Generalized Logistic ----
# ---------------------------------------------------------------------------- #

# xi = 0
# alpha = 1
# kappa = 2
# x <- lmomco::rlmomco(100, vec2par(c(xi, alpha, kappa), type = 'glo'))

# Ajuste por L-Momentos
AjustarLMomentosGeneralizedLogistic  <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'glo',
    parametros = list(xi = NA, alpha = NA, kappa = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (lmomco::are.lmom.valid(mom) && !is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      glo.pars <- lmomco::lmom2par(mom, type='glo')
      if (are.parglo.valid(glo.pars)) {
        ajuste$parametros$xi <- unname(glo.pars$para[1])
        ajuste$parametros$alpha <- unname(glo.pars$para[2])
        ajuste$parametros$kappa <- unname(glo.pars$para[3])
        
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudGeneralizedLogistic  <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'glo',
    parametros = list(xi = NA, alpha = NA, kappa = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    glo.pars.guess <- lmomco::parglo(lmomco.fit, checklmom = TRUE)
    start.params <- list(xi = unname(glo.pars.guess$para[1]), 
      alpha = unname(glo.pars.guess$para[2]),
      kappa = unname(glo.pars.guess$para[3]))
    
    # Realizar estimacion por metodo ML
    glo.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'glo',
      method = 'mle',
      keepdata = FALSE,
      start = start.params,
      lower = c(-Inf, 0, -Inf))
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$xi <- unname(glo.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(glo.mle.fit$estimate[2])
      ajuste$parametros$kappa <- unname(glo.mle.fit$estimate[3])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = glo.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo   <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$xi <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
      ajuste$parametros$kappa <- unname(parametros.remuestreo[3])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Generalized Normal -----
# ---------------------------------------------------------------------------- #

# xi = 0
# alpha = 1
# kappa = 0
# x <- lmomco::rlmomco(100, vec2par(c(xi, alpha, kappa), type = 'gno'))

# Ajuste por L-Momentos
AjustarLMomentosGeneralizedNormal  <- function(x, min.cantidad.valores) {
 # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'gno',
    parametros = list(xi = NA, alpha = NA, kappa = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (lmomco::are.lmom.valid(mom) && !is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      gno.pars <- lmomco::lmom2par(mom, type='gno')
      if (are.pargno.valid(gno.pars)) {
        ajuste$parametros$xi <- unname(gno.pars$para[1])
        ajuste$parametros$alpha <- unname(gno.pars$para[2])
        ajuste$parametros$kappa <- unname(gno.pars$para[3])
        
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudGeneralizedNormal  <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'gno',
    parametros = list(xi = NA, alpha = NA, kappa = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    gno.pars.guess <- lmomco::pargno(lmomco.fit, checklmom = TRUE)
    start.params   <- list(xi = unname(gno.pars.guess$para[1]), 
      alpha = unname(gno.pars.guess$para[2]),
      kappa = unname(gno.pars.guess$para[3]))
    
    # Realizar estimacion por metodo ML
    gno.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'gno',
      method = 'mle',
      keepdata = FALSE,
      start = start.params,
      lower = c(-Inf, 0, -1))
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$xi <- unname(gno.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(gno.mle.fit$estimate[2])
      ajuste$parametros$kappa <- unname(gno.mle.fit$estimate[3])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = gno.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$xi <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
      ajuste$parametros$kappa <- unname(parametros.remuestreo[3])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Govindarajulu -----
# ---------------------------------------------------------------------------- #

# xi = -3
# alpha = -2.5
# beta = 1
# x <- lmomco::rlmomco(100, vec2par(c(xi, alpha, kappa), type = 'gov'))

# Ajuste por L-Momentos
AjustarLMomentosGovindarajulu  <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'gov',
    parametros = list(xi = NA, alpha = NA, beta = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (lmomco::are.lmom.valid(mom) && !is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      gov.pars <- lmomco::lmom2par(mom, type='gov')
      if (are.pargov.valid(gov.pars)) {
        ajuste$parametros$xi <- unname(gov.pars$para[1])
        ajuste$parametros$alpha <- unname(gov.pars$para[2])
        ajuste$parametros$beta <- unname(gov.pars$para[3])
        
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudGovindarajulu  <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'gov',
    parametros = list(xi = NA, alpha = NA, beta = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    gov.pars.guess <- lmomco::pargov(lmomco.fit, checklmom = TRUE)
    start.params <- list(xi = unname(gov.pars.guess$para[1]), 
      alpha = unname(gov.pars.guess$para[2]),
      beta = unname(gov.pars.guess$para[3]))
    
    # Realizar estimacion por metodo ML
    gov.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'gov',
      method = 'mle',
      keepdata = FALSE,
      start = start.params)
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$xi <- unname(gov.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(gov.mle.fit$estimate[2])
      ajuste$parametros$beta <- unname(gov.mle.fit$estimate[3])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = gov.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$xi <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
      ajuste$parametros$beta <- unname(parametros.remuestreo[3])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Generalized Pareto ----
# ---------------------------------------------------------------------------- #

# xi = -0.07128146
# alpha = 1.22485715
# kappa = 0.02792960
# x <- lmomco::rlmomco(100, vec2par(c(xi, alpha, kappa), type = 'gpa'))

# Ajuste por L-Momentos
AjustarLMomentosGeneralizedPareto  <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'gpa',
    parametros = list(xi = NA, alpha = NA, kappa = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (lmomco::are.lmom.valid(mom) && !is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      gpa.pars <- lmomco::lmom2par(mom, type='gpa')
      if (are.pargpa.valid(gpa.pars)) {
        ajuste$parametros$xi <- unname(gpa.pars$para[1])
        ajuste$parametros$alpha <- unname(gpa.pars$para[2])
        ajuste$parametros$kappa <- unname(gpa.pars$para[3])
        
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudGeneralizedPareto  <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'gpa',
    parametros = list(xi = NA, alpha = NA, kappa = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    gpa.pars.guess <- lmomco::pargpa(lmomco.fit, checklmom = TRUE)
    start.params <- list(xi = unname(gpa.pars.guess$para[1]), 
      alpha = unname(gpa.pars.guess$para[2]),
      kappa = unname(gpa.pars.guess$para[3]))
    
    # Realizar estimacion por metodo ML
    gpa.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'gpa',
      method = 'mle',
      keepdata = FALSE,
      start = start.params,
      lower = c(-Inf, 0, -1))
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$xi <- unname(gpa.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(gpa.mle.fit$estimate[2])
      ajuste$parametros$kappa <- unname(gpa.mle.fit$estimate[3])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = gpa.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$xi <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
      ajuste$parametros$kappa <- unname(parametros.remuestreo[3])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Gumbel -----
# ---------------------------------------------------------------------------- #

# xi = -0.2975743
# alpha = 0.8613353
# x <- lmomco::rlmomco(100, vec2par(c(xi, alpha), type = 'gum'))

# Ajuste por L-Momentos
AjustarLMomentosGumbel  <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'gum',
    parametros = list(xi = NA, alpha = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (lmomco::are.lmom.valid(mom) && !is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      gum.pars <- lmomco::lmom2par(mom, type='gum')
      if (are.pargum.valid(gum.pars)) {
        ajuste$parametros$xi <- unname(gum.pars$para[1])
        ajuste$parametros$alpha <- unname(gum.pars$para[2])
        
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudGumbel  <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'gum',
    parametros = list(xi = NA, alpha = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    gum.pars.guess <- lmomco::pargum(lmomco.fit, checklmom = TRUE)
    start.params <- list(xi = unname(gum.pars.guess$para[1]), 
      alpha = unname(gum.pars.guess$para[2]))
    
    # Realizar estimacion por metodo ML
    gum.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'gum',
      method = 'mle',
      keepdata = FALSE,
      start = start.params,
      lower = c(-Inf, 0))
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$xi <- unname(gum.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(gum.mle.fit$estimate[2])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = gum.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$xi <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Kappa ----
# ---------------------------------------------------------------------------- #

# xi = -0.9133408
# alpha = 4.2596767
# kappa = 2.0039776
# h = 0.8568096
# x <- lmomco::rlmomco(100, vec2par(c(xi, alpha, kappa, h), type = 'kap'))

# Ajuste por L-Momentos
AjustarLMomentosKappa  <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'kap',
    parametros = list(xi = NA, alpha = NA, kappa = NA, h = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (lmomco::are.lmom.valid(mom) && !is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      kap.pars <- lmomco::lmom2par(mom, type='kap')
      if (are.parkap.valid(kap.pars)) {
        ajuste$parametros$xi <- unname(kap.pars$para[1])
        ajuste$parametros$alpha <- unname(kap.pars$para[2])
        ajuste$parametros$kappa <- unname(kap.pars$para[3])
        ajuste$parametros$h <- unname(kap.pars$para[4])
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudKappa  <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'kap',
    parametros = list(xi = NA, alpha = NA, kappa = NA, h = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    kap.pars.guess <- lmomco::parkap(lmomco.fit, checklmom = TRUE)
    start.params <- list(xi = unname(kap.pars.guess$para[1]), 
      alpha = unname(kap.pars.guess$para[2]),
      kappa = unname(kap.pars.guess$para[3]),
      h = unname(kap.pars.guess$para[4]))
    
    # Realizar estimacion por metodo ML
    kap.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'kap',
      method = 'mle',
      keepdata = FALSE,
      start = start.params,
      lower = c(-Inf, -Inf, -1, 0))
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$xi <- unname(kap.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(kap.mle.fit$estimate[2])
      ajuste$parametros$kappa <- unname(kap.mle.fit$estimate[3])
      ajuste$parametros$h <- unname(kap.mle.fit$estimate[4])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = kap.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$xi <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
      ajuste$parametros$kappa <- unname(parametros.remuestreo[3])
      ajuste$parametros$h <- unname(parametros.remuestreo[4])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Kumaraswamy ----
# ---------------------------------------------------------------------------- #

# alpha = 0.5
# beta = 0.5
# x <- lmomco::rlmomco(100, vec2par(c(alpha, beta), type = 'kur'))

# Ajuste por L-Momentos
AjustarLMomentosKumaraswamy   <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'kur',
    parametros = list(alpha = NA, beta = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (lmomco::are.lmom.valid(mom) && !is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      kur.pars <- lmomco::lmom2par(mom, type='kur')
      if (are.parkur.valid(kur.pars)) {
        ajuste$parametros$alpha <- unname(kur.pars$para[1])
        ajuste$parametros$beta <- unname(kur.pars$para[2])
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudKumaraswamy  <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'kur',
    parametros = list(alpha = NA, beta = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    kur.pars.guess <- lmomco::parkur(lmomco.fit, checklmom = TRUE)
    start.params <- list(alpha = unname(kur.pars.guess$para[1]), 
      beta = unname(kur.pars.guess$para[2]))
    
    # Realizar estimacion por metodo ML
    kur.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'kur',
      method = 'mle',
      keepdata = FALSE,
      start = start.params,
      lower = c(0, 0))
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$alpha <- unname(kur.mle.fit$estimate[1])
      ajuste$parametros$beta <- unname(kur.mle.fit$estimate[2])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = kur.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$alpha <- unname(parametros.remuestreo[1])
      ajuste$parametros$beta <- unname(parametros.remuestreo[2])
      
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Laplace ----
# xi = 0.4921780
# alpha = 0.7485578
# x <- lmomco::rlmomco(100, vec2par(c(xi, alpha), type = 'lap'))

# Ajuste por L-Momentos
AjustarLMomentosLaplace   <- function(x, min.cantidad.valores) {
 # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'lap',
    parametros = list(xi = NA, alpha = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (lmomco::are.lmom.valid(mom) && !is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      lap.pars <- lmomco::lmom2par(mom, type='lap')
      if (are.parlap.valid(lap.pars)) {
        ajuste$parametros$xi <- unname(lap.pars$para[1])
        ajuste$parametros$alpha <- unname(lap.pars$para[2])
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudLaplace  <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'lap',
    parametros = list(xi = NA, alpha = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    lap.pars.guess <- lmomco::parlap(lmomco.fit, checklmom = TRUE)
    start.params <- list(xi = unname(lap.pars.guess$para[1]), 
      alpha = unname(lap.pars.guess$para[2]))
    
    # Realizar estimacion por metodo ML
    lap.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'lap',
      method = 'mle',
      keepdata = FALSE,
      start = start.params)
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$xi <- unname(lap.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(lap.mle.fit$estimate[2])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = lap.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$xi <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
      
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Linear Mean Residual Quantile ----
# mu = 1.30714286
# alpha = -0.04714286
# x <- lmomco::rlmomco(100, vec2par(c(mu, alpha), type = 'lmrq'))

# Ajuste por L-Momentos
AjustarLMomentosLinearMeanResidualQuantile    <- function(x, min.cantidad.valores) {
 # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'lmrq',
    parametros = list(mu = NA, alpha = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (lmomco::are.lmom.valid(mom) && !is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      lmrq.pars <- lmomco::lmom2par(mom, type='lmrq')
      if (are.parlmrq.valid(lmrq.pars)) {
        ajuste$parametros$mu <- unname(lmrq.pars$para[1])
        ajuste$parametros$alpha <- unname(lmrq.pars$para[2])
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudLinearMeanResidualQuantile   <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'lmrq',
    parametros = list(mu = NA, alpha = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    lmrq.pars.guess <- lmomco::parlmrq(lmomco.fit, checklmom = TRUE)
    start.params <- list(mu = unname(lmrq.pars.guess$para[1]), 
      alpha = unname(lmrq.pars.guess$para[2]))
    
    # Realizar estimacion por metodo ML
    lmqr.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'lmrq',
      method = 'mle',
      keepdata = FALSE,
      start = start.params)
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$mu <- unname(lmqr.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(lmqr.mle.fit$estimate[2])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = lmqr.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$mu <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
      
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar 3-Parameter Log-Normal ----
# ---------------------------------------------------------------------------- #

# zeta <- 0
# meanlog <- 1
# sdlog <- 0.7
# set.seed(128)
# x <- lmomco::rlmomco(n = 100, vec2par(c(zeta, meanlog, sdlog), type = 'ln3'))

# Ajuste por L-Momentos
AjustarLMomentosLogNormal <- function(x, min.cantidad.valores) {
 # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'ln3',
    parametros = list(zeta = 0, mulog = NA, sigmalog = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (!is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      ln3.pars <- lmomco::lmom2par(mom, type='ln3', zeta = 0)
      if (are.parln3.valid(ln3.pars)) {
        ajuste$parametros$zeta <- 0
        ajuste$parametros$mulog <- unname(ln3.pars$para[2])
        ajuste$parametros$sigmalog <- unname(ln3.pars$para[3])
        
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudLogNormal   <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
 # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'ln3',
    parametros = list(zeta = NA, mulog = NA, sigmalog = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    ln3.pars.guess <- lmomco::parln3(lmomco.fit, checklmom = TRUE, zeta = 0)
    start.params <- list(mulog = unname(ln3.pars.guess$para[2]),
      sigmalog = unname(ln3.pars.guess$para[3]))
    
    # Realizar estimacion por metodo ML
    ln3.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'ln3',
      method = 'mle',
      keepdata = FALSE,
      start = start.params)
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$zeta <- 0
      ajuste$parametros$mulog <- unname(ln3.mle.fit$estimate[1])
      ajuste$parametros$sigmalog <- unname(ln3.mle.fit$estimate[2])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = ln3.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$zeta <- 0
      ajuste$parametros$mulog <- unname(parametros.remuestreo[1])
      ajuste$parametros$sigmalog <- unname(parametros.remuestreo[2])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Normal ----
# ---------------------------------------------------------------------------- #

# mean <- 0
# sd <- 1
# set.seed(128)
# x <- lmomco::rlmomco(n = 1000, vec2par(c(mean, sd), type = 'nor'))

# Ajuste por L-Momentos
AjustarLMomentosNormal <- function(x, min.cantidad.valores) {
  # Ajuste de la distribución Gamma por metodo de L-Momentos
  
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'nor',
    parametros = list(mu = NA, sigma = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (!is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      nor.pars <- lmomco::lmom2par(mom, type='nor')
      if (are.parnor.valid(nor.pars)) {
        ajuste$parametros$mu <- unname(nor.pars$para[1])
        ajuste$parametros$sigma <- unname(nor.pars$para[2])
        
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudNormal   <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Ajuste de la distribución Gamma por metodo de maxima verosimilitud
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'nor',
    parametros = list(mu = NA, sigma = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    nor.pars.guess <- lmomco::parnor(lmomco.fit, checklmom = TRUE)
    start.params <- list(mean = unname(nor.pars.guess$para[1]), 
      sd = unname(nor.pars.guess$para[2]))
    
    # Realizar estimacion por metodo ML
    nor.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'norm',
      method = 'mle',
      keepdata = FALSE,
      start = start.params)
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$mu <- unname(nor.mle.fit$estimate[1])
      ajuste$parametros$sigma <- unname(nor.mle.fit$estimate[2])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = nor.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo   <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$mu <- unname(parametros.remuestreo[1])
      ajuste$parametros$sigma <- unname(parametros.remuestreo[2])
      
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Ajustar Pearson Type III ----
# ---------------------------------------------------------------------------- #

# mu <- 0.1136436
# sigma <- 1.2043534
# gamma <- 0.6730858 
# set.seed(128)
# x <- lmomco::rlmomco(n = 100, vec2par(c(mu, sigma, gamma), type = 'pe3'))

# Ajuste por L-Momentos
AjustarLMomentosPearsonIII <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'pe3',
    parametros = list(mu = NA, sigma = NA, gamma = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (!is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      pe3.pars <- lmomco::lmom2par(mom, type='pe3')
      if (are.parpe3.valid(pe3.pars)) {
        ajuste$parametros$mu <- unname(pe3.pars$para[1])
        ajuste$parametros$sigma <- unname(pe3.pars$para[2])
        ajuste$parametros$gamma <- unname(pe3.pars$para[3])
        
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudPearsonIII    <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'pe3',
    parametros = list(mu = NA, sigma = NA, gamma = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    pe3.pars.guess <- lmomco::parpe3(lmomco.fit, checklmom = TRUE)
    start.params <- list(mu = unname(pe3.pars.guess$para[1]), 
      sigma = unname(pe3.pars.guess$para[2]),
      gamma = unname(pe3.pars.guess$para[3]))
    
    # Realizar estimacion por metodo ML
    pe3.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'pe3',
      method = 'mle',
      keepdata = FALSE,
      start = start.params,
      lower = c(-Inf, 0, 0))
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$mu <- unname(pe3.mle.fit$estimate[1])
      ajuste$parametros$sigma <- unname(pe3.mle.fit$estimate[2])
      ajuste$parametros$gamma <- unname(pe3.mle.fit$estimate[3])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = pe3.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$mu <- unname(parametros.remuestreo[1])
      ajuste$parametros$sigma <- unname(parametros.remuestreo[2])
      ajuste$parametros$gamma <- unname(parametros.remuestreo[3])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------ 

# ---------------------------------------------------------------------------- #
# Ajustar Rayleigh ----
# ---------------------------------------------------------------------------- #

# xi <- -2.060577
# alpha <- 1.824173
# set.seed(128)
# x <- lmomco::rlmomco(n = 100, vec2par(c(xi, alpha), type = 'ray'))

# Ajuste por L-Momentos
AjustarLMomentosRayleigh  <- function(x, min.cantidad.valores) {
 # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'ray',
    parametros = list(xi = NA, alpha = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (!is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      ray.pars <- lmomco::lmom2par(mom, type='ray')
      if (are.parray.valid(ray.pars)) {
        ajuste$parametros$xi <- unname(ray.pars$para[1])
        ajuste$parametros$alpha <- unname(ray.pars$para[2])
        
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudRayleigh   <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'ray',
    parametros = list(xi = NA, alpha = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    ray.pars.guess <- lmomco::parray(lmomco.fit, checklmom = TRUE)
    start.params <- list(xi = unname(ray.pars.guess$para[1]), 
      alpha = unname(ray.pars.guess$para[2]))
    
    # Realizar estimacion por metodo ML
    ray.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'ray',
      method = 'mle',
      keepdata = FALSE,
      start = start.params,
      lower = c(-Inf, 0))
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$xi <- unname(ray.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(ray.mle.fit$estimate[2])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = ray.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$xi <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
      
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------ 

# ---------------------------------------------------------------------------- #
# Ajustar Rice ----
# ---------------------------------------------------------------------------- #

# nu <- 0
# alpha <- 10
# set.seed(128)
# x <- lmomco::rlmomco(n = 100, vec2par(c(nu, alpha), type = 'rice'))

# Ajuste por L-Momentos
AjustarLMomentosRice  <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'rice',
    parametros = list(nu = NA, alpha = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (!is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      rice.pars <- lmomco::lmom2par(mom, type='rice')
      if (are.parrice.valid(rice.pars)) {
        ajuste$parametros$nu <- unname(rice.pars$para[1])
        ajuste$parametros$alpha <- unname(rice.pars$para[2])
        
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudRice   <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'rice',
    parametros = list(nu = NA, alpha = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    rice.pars.guess <- lmomco::parrice(lmomco.fit, checklmom = TRUE)
    start.params <- list(nu = unname(rice.pars.guess$para[1]), 
      alpha = unname(rice.pars.guess$para[2]))
    
    # Realizar estimacion por metodo ML
    rice.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'rice',
      method = 'mle',
      keepdata = FALSE,
      start = start.params,
      lower = c(0, 0))
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$nu <- unname(rice.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(rice.mle.fit$estimate[2])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = rice.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$nu <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
      
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------ 

# ---------------------------------------------------------------------------- #
# Ajustar Slash ----
# ---------------------------------------------------------------------------- #

# xi <- -100
# alpha <- 30
# set.seed(128)
# x <- lmomco::rlmomco(n = 100, vec2par(c(xi, alpha), type = 'sla'))

# Ajuste por L-Momentos
AjustarLMomentosSlash  <- function(x, min.cantidad.valores) {
  # Ajuste de la distribución Gamma por metodo de L-Momentos
  
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'sla',
    parametros = list(xi = NA, alpha = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::TLmoms(x, trim = 1)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (!is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      sla.pars <- lmomco::lmom2par(mom, type='sla')
      if (are.parsla.valid(sla.pars)) {
        ajuste$parametros$xi <- unname(sla.pars$para[1])
        ajuste$parametros$alpha <- unname(sla.pars$para[2])
        
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudSlash   <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'sla',
    parametros = list(xi = NA, alpha = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::TLmoms(x, trim = 1)
    sla.pars.guess <- lmomco::parsla(lmomco.fit, checklmom = TRUE)
    start.params <- list(xi = unname(sla.pars.guess$para[1]), 
      alpha = unname(sla.pars.guess$para[2]))
    
    # Realizar estimacion por metodo ML
    sla.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'sla',
      method = 'mle',
      keepdata = FALSE,
      start = start.params)
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$xi <- unname(sla.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(sla.mle.fit$estimate[2])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = sla.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$xi <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
      
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------ 

# ---------------------------------------------------------------------------- #
# Ajustar 3-Parameter Student t ----
# ---------------------------------------------------------------------------- #

# xi <- 10
# alpha <- 1.5
# nu <- 2
# set.seed(128)
# x <- lmomco::rlmomco(n = 100, vec2par(c(xi, alpha, nu), type = 'st3'))

# Ajuste por L-Momentos
AjustarLMomentosStudent_t3  <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'st3',
    parametros = list(xi = NA, alpha = NA, nu = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmoms(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (!is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      st3.pars <- lmomco::lmom2par(mom, type='st3')
      if (are.parst3.valid(st3.pars)) {
        ajuste$parametros$xi <- unname(st3.pars$para[1])
        ajuste$parametros$alpha <- unname(st3.pars$para[2])
        ajuste$parametros$nu <- unname(st3.pars$para[3])
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudStudent_t3  <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Ajuste de la distribución Gamma por metodo de maxima verosimilitud
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'st3',
    parametros = list(xi = NA, alpha = NA, nu = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    st3.pars.guess <- lmomco::parst3(lmomco.fit, checklmom = TRUE)
    start.params <- list(xi = unname(st3.pars.guess$para[1]), 
      alpha = unname(st3.pars.guess$para[2]),
      nu = unname(st3.pars.guess$para[3]))
    
    # Realizar estimacion por metodo ML
    st3.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'st3',
      method = 'mle',
      keepdata = FALSE,
      start = start.params,
      lower = c(-Inf, 0, 1))
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$xi <- unname(st3.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(st3.mle.fit$estimate[2])
      ajuste$parametros$nu <- unname(st3.mle.fit$estimate[3])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = st3.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$xi <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
      ajuste$parametros$nu <- unname(parametros.remuestreo[3])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------ 

# ---------------------------------------------------------------------------- #
# Ajustar Asymmetric Triangular ----
# ---------------------------------------------------------------------------- #

# min <- 10
# mode <- 90
# max <- 100
# set.seed(128)
# x <- lmomco::rlmomco(n = 100, vec2par(c(min, mode, max), type = 'tri'))

# Ajuste por L-Momentos
AjustarLMomentosAsymmetricTriangular  <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'tri',
    parametros = list(min = NA, mode = NA, max = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmoms(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (!is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      tri.pars <- lmomco::lmom2par(mom, type='tri')
      if (are.partri.valid(tri.pars)) {
        ajuste$parametros$min <- unname(tri.pars$para[1])
        ajuste$parametros$mode <- unname(tri.pars$para[2])
        ajuste$parametros$max <- unname(tri.pars$para[3])
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudAsymmetricTriangular  <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'tri',
    parametros = list(min = NA, mode = NA, max = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    tri.pars.guess <- lmomco::partri(lmomco.fit, checklmom = TRUE)
    start.params <- list(min = unname(tri.pars.guess$para[1]), 
      mode = unname(tri.pars.guess$para[2]),
      max = unname(tri.pars.guess$para[3]))
    
    # Realizar estimacion por metodo ML
    tri.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'tri',
      method = 'mle',
      keepdata = FALSE,
      start = start.params)
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$min <- unname(tri.mle.fit$estimate[1])
      ajuste$parametros$mode <- unname(tri.mle.fit$estimate[2])
      ajuste$parametros$max <- unname(tri.mle.fit$estimate[3])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = tri.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$min <- unname(parametros.remuestreo[1])
      ajuste$parametros$mode <- unname(parametros.remuestreo[2])
      ajuste$parametros$max <- unname(parametros.remuestreo[3])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------ 

# ---------------------------------------------------------------------------- #
# Ajustar Wakeby ----
# ---------------------------------------------------------------------------- #

# xi <- -2
# alpha <- 4
# beta <- 6
# gamma <- 1.5
# delta <- -0.25
# set.seed(128)
# x <- lmomco::rlmomco(n = 100, vec2par(c(xi, alpha, beta, gamma, delta), type = 'wak'))

# Ajuste por L-Momentos
AjustarLMomentosWakeby  <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'wak',
    parametros = list(xi = NA, alpha = NA, beta = NA,
      gamma = NA, delta = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (!is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      wak.pars <- lmomco::lmom2par(mom, type='wak')
      if (are.parwak.valid(wak.pars)) {
        ajuste$parametros$xi <- unname(wak.pars$para[1])
        ajuste$parametros$alpha <- unname(wak.pars$para[2])
        ajuste$parametros$beta <- unname(wak.pars$para[3])
        ajuste$parametros$gamma <- unname(wak.pars$para[4])
        ajuste$parametros$delta <- unname(wak.pars$para[5])
        
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudWakeby  <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'wak',
    parametros = list(xi = NA, alpha = NA, beta = NA,
      gamma = NA, delta = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    wak.pars.guess <- lmomco::parwak(lmomco.fit, checklmom = TRUE)
    start.params <- list(xi = unname(wak.pars.guess$para[1]), 
      alpha = unname(wak.pars.guess$para[2]),
      beta = unname(wak.pars.guess$para[3]),
      gamma = unname(wak.pars.guess$para[4]),
      delta = unname(wak.pars.guess$para[5]))
    
    # Realizar estimacion por metodo ML
    wak.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'wak',
      method = 'mle',
      keepdata = FALSE,
      start = start.params)
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$xi <- unname(wak.mle.fit$estimate[1])
      ajuste$parametros$alpha <- unname(wak.mle.fit$estimate[2])
      ajuste$parametros$beta <- unname(wak.mle.fit$estimate[3])
      ajuste$parametros$gamma <- unname(wak.mle.fit$estimate[4])
      ajuste$parametros$delta <- unname(wak.mle.fit$estimate[5])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = wak.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$xi <- unname(parametros.remuestreo[1])
      ajuste$parametros$alpha <- unname(parametros.remuestreo[2])
      ajuste$parametros$beta <- unname(parametros.remuestreo[3])
      ajuste$parametros$gamma <- unname(parametros.remuestreo[4])
      ajuste$parametros$delta <- unname(parametros.remuestreo[5])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------ 

# ---------------------------------------------------------------------------- #
# Ajustar Weibull ----
# ---------------------------------------------------------------------------- #

# zeta <- 3
# beta <- 2.5
# alpha <- 2.5
# set.seed(128)
# x <- lmomco::rlmomco(n = 100, vec2par(c(zeta, beta, alpha), type = 'wei'))


# Ajuste por L-Momentos
AjustarLMomentosWeibull  <- function(x, min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'wei',
    parametros = list(zeta = NA, beta = NA,  delta = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar Unbiased Sample L-moments
    mom <- lmomco::lmom.ub(x)
    
    # Convertir Probability-Weighted Moments (PWM) a L-moments
    if (!is.na(sum(mom[[1]])) && !is.nan(sum(mom[[1]]))) {
      # Extraer parametros
      wei.pars <- lmomco::lmom2par(mom, type='wei')
      if (are.parwei.valid(wei.pars)) {
        ajuste$parametros$zeta <- unname(wei.pars$para[1])
        ajuste$parametros$beta <- unname(wei.pars$para[2])
        ajuste$parametros$delta <- unname(wei.pars$para[3])
        
      }
    } 
  }
  
  return (ajuste)
}

# Ajuste por Maxima verosimilitud
AjustarMaximaVerosimilitudWeibull  <- function(x, numero.muestras = NULL, 
  min.cantidad.valores) {
  # Verificar si hay un minimo numero de valores para estimar parametros
  ajuste <- list(distribucion = 'wei',
    parametros = list(zeta = NA, beta = NA,  delta = NA))
  if (length(x) >= (min.cantidad.valores)) {
    # Estimar L-momentos para usar como valores iniciales
    lmomco.fit <- lmomco::pwm2lmom(lmomco::pwm.ub(x))
    wei.pars.guess <- lmomco::parwei(lmomco.fit, checklmom = TRUE)
    start.params <- list(zeta = unname(wei.pars.guess$para[1]), 
      beta = unname(wei.pars.guess$para[2]),
      delta = unname(wei.pars.guess$para[3]))
    
    # Realizar estimacion por metodo ML
    wei.mle.fit <- fitdistrplus::fitdist(data = x,
      distr = 'wei',
      method = 'mle',
      keepdata = FALSE,
      start = start.params,
      lower = c(-Inf, 0, 0))
    if (is.null(numero.muestras)) {
      # Estimar parametros sin remuestreo
      ajuste$parametros$zeta <- unname(wei.mle.fit$estimate[1])
      ajuste$parametros$beta <- unname(wei.mle.fit$estimate[2])
      ajuste$parametros$delta <- unname(wei.mle.fit$estimate[3])
    } else {
      # Realizar bootstrap parametrico para los parametros estimados
      bootstrap.obj <- fitdistrplus::bootdist(f = wei.mle.fit,
        bootmethod = "param",
        niter = numero.muestras,
        silent = TRUE)
      
      # Calcular la mediana de cada parametro
      parametros.remuestreo <- apply(X = bootstrap.obj$estim, MARGIN = 2, FUN = median, na.rm = TRUE)
      ajuste$parametros$zeta <- unname(parametros.remuestreo[1])
      ajuste$parametros$beta <- unname(parametros.remuestreo[2])
      ajuste$parametros$delta <- unname(parametros.remuestreo[3])
    }
  }
  
  return (ajuste)
}

# ------------------------------------------------------------------------------