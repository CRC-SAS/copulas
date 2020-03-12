# Definicion de funcion para aleatorizar las variables ----
AgregarRuido <- function(x, valor.minimo.deteccion) {
  # Agregar ruido con una distribución uniforme cuyos umbrales son la 
  # mínima unidad de medida de la variable.
  minimo <- min(0, valor.minimo.deteccion)
  maximo <- max(0, valor.minimo.deteccion)
  
  repeticiones.ruido <- x + valor.minimo.deteccion * runif(length(x), min  = minimo, max = maximo)
  
  # Devolver objeto resultado
  return(repeticiones.ruido)
  
}
#-------------------------------------------------------------------------------

# Definicion de funcion para graficar Chi plot ----
GraficarChiPlot <- function(x, y) {
  # Crear lista con las variables a graficar
  datos.graficar <- CDVine::BiCopChiPlot(pobs(x), pobs(y), PLOT = F)
  
  # Convertir variables a puntos
  lambda.chi <- data.frame(lambda = datos.graficar$lambda,
    chi = datos.graficar$chi)
  # Graficar Chi plot
  ggplot2::ggplot(data = lambda.chi, aes(x = lambda, y = chi)) +
    ggplot2::geom_point() +
    geom_hline(yintercept = datos.graficar$control.bounds[1], linetype="dashed") +
    geom_hline(yintercept = datos.graficar$control.bounds[2], linetype="dashed") +
    geom_vline(xintercept = 0, linetype="dashed") +
    theme_bw() +
    ylab(expression(chi)) + xlab(expression(lambda)) +
    ggsave("Chiplot.png", device = 'png', dpi = 600)
  
}
#-------------------------------------------------------------------------------

# Definicion de funcion para graficar Kendall plot ----
GraficarKPlot <- function(x, y) {
  # Crear lista con las variables a graficar
  datos.graficar <- CDVine::BiCopKPlot(copula::pobs(x), copula::pobs(y), PLOT = F)
  
  # Convertir variables a puntos
  W.Hisort <- data.frame(W = datos.graficar$W.in,
    H = datos.graficar$Hi.sort)
  # Graficar Chi plot
  ggplot2::ggplot(data = W.Hisort, aes(x = W, y = H)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    stat_function(fun = function(x) x - x * log(x), geom = 'line') +
    theme_bw() +
    xlim(0, 1) + ylim(0, 1) +
    xlab(expression(W[1:n])) + ylab('H') +
    ggsave("Kplot.png", device = 'png', dpi = 600)
  
}
#-------------------------------------------------------------------------------

