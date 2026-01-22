#'@title Anomaly detector using ARIMA.
#'@description Anomaly detection using ARIMA
#'The ARIMA model adjusts to the time series. Observations distant from the model are labeled as anomalies.
#'It wraps the ARIMA model presented in the forecast library.
#'@return `arima_filter` object
#'@export
arima_filter <- function(filtro) {
  obj <- harbinger()
  obj$filtro <- filtro
  hutils <- harutils()
  obj$har_distance <- hutils$har_distance_l1
  obj$sw_size <- NULL
  
  class(obj) <- append("arima_filter", class(obj))
  return(obj)
}

#'@importFrom forecast auto.arima
#'@importFrom stats residuals
#'@importFrom stats na.omit
#'@export
fit.arima_filter <- function(obj, serie, ...) {
  if (is.null(serie)) stop("No data was provided for computation", call. = FALSE)
  
  serie <- stats::na.omit(serie)
  
  # Aplicar o filtro armazenado no objeto
  if (!is.null(obj$filtro)) {
    serie <- transform(obj$filtro, serie)
  }
  
  obj$model <- forecast::auto.arima(serie, allowdrift = TRUE, allowmean = TRUE)
  order <- obj$model$arma[c(1, 6, 2, 3, 7, 4, 5)]
  obj$p <- order[1]
  obj$d <- order[2]
  obj$q <- order[3]
  obj$drift <- (NCOL(obj$model$xreg) == 1) && is.element("drift", names(obj$model$coef))
  params <- list(p = obj$p, d = obj$d, q = obj$q, drift = obj$drift)
  attr(obj, "params") <- params
  
  if (is.null(obj$sw_size))
    obj$sw_size <- max(obj$p, obj$d + 1, obj$q)
  
  return(obj)
}

#'@importFrom forecast auto.arima
#'@importFrom stats residuals
#'@importFrom stats na.omit
#'@export
detect.arima_filter <- function(obj, serie, ...) {
  if (is.null(serie)) stop("No data was provided for computation", call. = FALSE)
  
  obj <- obj$har_store_refs(obj, serie)
  
  # Aplicar o filtro armazenado no objeto
  
  model_serie <- obj$serie
  if (!is.null(obj$filtro)) {
    model_serie_f <- try(transform(obj$filtro, model_serie), silent = TRUE)
    if (!inherits(model_serie_f, "try-error")) {
      model_serie <- model_serie_f
    }
  }

  # Garantias defensivas: alguns filtros podem inserir NA nas bordas ou alterar comprimento
  model_serie <- as.numeric(model_serie)
  ref_serie <- as.numeric(obj$serie)
  if (length(model_serie) != length(ref_serie)) {
    model_serie <- ref_serie
  }
  if (anyNA(model_serie)) {
    model_serie[is.na(model_serie)] <- ref_serie[is.na(model_serie)]
  }
  
  # Ajustando o modelo para a série completa
  model <- tryCatch(
    {
      forecast::Arima(model_serie, order = c(obj$p, obj$d, obj$q), include.drift = obj$drift)
    },
    error = function(cond) {
      forecast::auto.arima(model_serie, allowdrift = TRUE, allowmean = TRUE)
    }
  )
  
  pred <- model$x - model$residuals

  pred <- as.numeric(pred)
  # Alinhar comprimentos para evitar falha em séries específicas
  if (length(pred) != length(ref_serie)) {
    n <- min(length(pred), length(ref_serie))
    pred_tail <- tail(pred, n)
    serie_tail <- tail(ref_serie, n)
    res_tail <- serie_tail - pred_tail
    if (n < length(ref_serie)) {
      res <- c(rep(NA_real_, length(ref_serie) - n), res_tail)
    } else {
      res <- res_tail
    }
  } else {
    res <- ref_serie - pred
  }
  
  res <- obj$har_distance(res)
  anomalies <- obj$har_outliers(res)
  anomalies <- obj$har_outliers_check(anomalies, res)
  
  anomalies[1:obj$sw_size] <- FALSE
  
  detection <- obj$har_restore_refs(obj, anomalies = anomalies, res = res)
  
  return(detection)
}

