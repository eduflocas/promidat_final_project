
# RSS: Residual Sum of Squares
RSS <- function(Pred, Real) {
  return(sum((Real - Pred)^2))
}

# MSE: Mean Squared Error

MSE <- function(Pred, Real) {
  N <- length(Real)
  rss <- sum((Real - Pred)^2)
  return((1/N) * rss)
}

# RMSE: Root Mean Squared Error

RMSE <- function(Pred, Real) {
  N <- length(Real)
  rss <- sum((Real - Pred)^2)
  return(sqrt((1/N) * rss))
}

# PFA: Porcentaje de veces en las que el pronóstico fue mayor o igual a la realidad

PFA <- function(Pred, Real) {
  Total <- 0
  N <- length(Pred)
  for(i in 1:N) {
    if(Pred[i] > Real[i])
      Total <- Total + 1      
  }
  return(Total/N)
}

# PTFA: Porcentaje de fallos hacia arriba en términos absolutos

PTFA <- function(Pred, Real) {
  Total <- 0
  SReal <- 0
  N <- length(Pred)
  for(i in 1:N) {
    if(Pred[i] > Real[i]) {
      Total <- Total + (Pred[i] - Real[i])
      SReal <- SReal + abs(Real[i])
    }
  }
  if(Total == 0)
    SReal = 1
  return(Total/SReal)
}

# Obtención de todos los errores

tabla.errores <- function(predicciones, real, nombres = NULL) {
  r <- data.frame()
  for (pred in predicciones) {
    r <- rbind(r, data.frame(
      'MSE'    = MSE(pred, real), 'RMSE'   = RMSE(pred, real),
      'PFA'    = PFA(pred, real), 'PTFA'   = PTFA(pred, real)
    )
    )
  }
  row.names(r) <- nombres
  return(r)
}

# tabla.errores <- function(predicciones, real) {
#   r <- data.frame()
#   for (pred in predicciones) {
#     r <- rbind(r, data.frame(MSE = MSE(pred, real), RMSE = RMSE(pred, real), 
#                              PFA = PFA(pred, real), PTFA = PTFA(pred, real)))
#   }
#   return(r)
# }


#Gráfico radar

grafico.errores <- function (errores) {
  library(ggplot2)
  library(reshape)
  library(scales)
  
  centros <- as.data.frame(apply(errores, 2, function(i)
    scales::rescale(i, to = c(0, 100))))
  
  res <- melt(t(centros), varnames = c("E", "Modelos"))
  res <- res[order(res$E, decreasing = F), ]
  res$M <- as.character(res$M)
  y = c(0, 25, 50, 75, 100)
  
  ggplot(res, aes(x = E, y = value, group = Modelos, color = Modelos, fill = Modelos)) +
    geom_polygon(alpha = 0.3, size = 1) + geom_point(size = 3) + 
    theme_minimal() + theme(axis.text.y = element_blank()) + xlab("") + 
    ylab("") + scale_y_continuous(limits = c(-10, 100), breaks = y) + 
    annotate("text", x = 0.5, y = y, label = paste0(y, "%"), color = "gray60") +
    ggproto("CordRadar", CoordPolar, theta = "x", r = "y", 
            start = 0, direction = sign(1))
}


# Calibración Holt-Winters

calibrar.HW <- function(entrenamiento, prueba, paso = 0.1) {
  # se calculan todas las combinaciones para los parametros
  params <- purrr::cross(list(a = seq(0, 1, by = paso), b = seq(0, 1, by = paso), 
                              g = seq(0, 1, by = paso)))
  
  # se calcula un modelos para cada combinacion de parametros
  hw_secure <- purrr::possibly(stats::HoltWinters, otherwise = NULL)
  models <- purrr::map(params, ~suppressWarnings(hw_secure(entrenamiento, 
                                                           alpha = ifelse(.$a == 0, F, .$a), beta = ifelse(.$b == 0, F, .$b), gamma = ifelse(.$g == 
                                                                                                                                               0, F, .$g))))
  
  # se realiza la prediccion con cada uno de los modelos
  predictions <- purrr::map(models, ~{
    if (is.null(.)) {
      return(NULL)
    }
    forecast(., h = length(prueba))
  })
  
  # se calcula el error para cada prediccion
  error <- purrr::map_dbl(predictions, ~{
    if (is.null(.)) {
      return(Inf)
    }
    sum((as.numeric(prueba) - as.numeric(.$mean))^2)
  })
  
  # se retorna el modelo con el menor error
  best_model <- models[[which.min(error)]]
  p <- params[[which.min(error)]]
  best_model$call <- call("HoltWinters", x = quote(datos), alpha = ifelse(p$a == 
                                                                            0, F, p$a), beta = ifelse(p$b == 0, F, p$b), gamma = ifelse(p$g == 0, 
                                                                                                                                        F, p$g))
  return(best_model)
}

#Calibrar Arima

calibrar.arima <- function(entrenamiento = NULL, prueba = NULL, periodo = NA_integer_, ar = 0:2, es = 0:1) {
  # se calculan todas las combinaciones para los parametros
  params <- purrr::cross(list(a = ar, b = ar, c = ar, d = es, e = es, f = es))
  
  # se calcula un modelos para cada combinacion de parametros
  arima_secure <- purrr::possibly(stats::arima, otherwise = NULL)
  models <- purrr::map(params, ~ suppressWarnings(arima_secure(
    entrenamiento, order = c(.$a,.$b,.$c), 
    seasonal = list(order = c(.$d,.$e,.$f), period = periodo))))
  
  # se realiza la prediccion con cada uno de los modelos
  predictions <- purrr::map(models, ~{
    if (is.null(.)) {return(NULL)}
    forecast(., h = length(prueba))
  })
  
  # se calcula el error para cada prediccion
  error <- purrr::map_dbl(predictions, ~{
    if(is.null(.)) {return(Inf)}
    sum((as.numeric(prueba) - as.numeric(.$mean))^2)
  })
  
  # se retorna el modelo con el menor error
  best_model <- models[[which.min(error)]]
  p <- params[[which.min(error)]]
  best_model$call <- call(
    "arima", x = quote(datos), order = as.numeric(c(p$a, p$b, p$c)),
    seasonal = list(order = as.numeric(c(p$d, p$e, p$f)), period = periodo))
  return(best_model)
}
