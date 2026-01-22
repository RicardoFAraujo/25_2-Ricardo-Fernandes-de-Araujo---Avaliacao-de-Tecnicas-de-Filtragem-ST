############################################################
# Análise Benchmarks Yahoo (A1..A4) com aplicação de filtros e ARIMA
# Objetivo: testar múltiplos filtros TS e detectar anomalias via ARIMA
# Autor: Ricardo (assistido)
# Data: Sys.Date()
############################################################

suppressPackageStartupMessages({
  library(united)
  library(daltoolbox)
  library(harbinger)
  library(dplyr)
})

source("arima_filter.R")

YAHOO_BENCHMARKS <- c("A1","A2","A3","A4")
yahoo_obj_name <- function(code) paste0(code, "Benchmark")

load_yahoo_benchmarks <- function(benchmarks = YAHOO_BENCHMARKS) {
  for (b in benchmarks) {
    nm <- yahoo_obj_name(b)
    if (!exists(nm, envir = .GlobalEnv, inherits = FALSE)) {
      tryCatch({
        data(list = nm, package = "united", envir = .GlobalEnv)
        if (exists(nm, envir = .GlobalEnv, inherits = FALSE)) {
          message("Carregado: ", nm)
        } else {
          warning("Falha ao carregar ", nm, ": não encontrado no pacote 'united'.")
        }
      }, error = function(e) warning("Erro ao carregar ", nm, ": ", conditionMessage(e)))
    } else {
      message("Já existente: ", nm)
    }
  }
  invisible(TRUE)
}

find_label_column <- function(df) {
  candidates <- c("anomaly","is_anomaly","label","IsAnomaly","anomalous")
  hits <- candidates[candidates %in% names(df)]
  if (!length(hits)) return(NA_character_)
  if (length(hits) > 1) {
    for (h in hits) {
      col <- df[[h]]
      if (is.logical(col)) return(h)
      if (is.numeric(col) && all(na.omit(unique(col)) %in% c(0,1))) return(h)
    }
  }
  hits[1]
}

prepare_yahoo_series <- function(benchmark = "A1", series_index = 1, value_col = NA) {
  load_yahoo_benchmarks(benchmark)
  obj_name <- yahoo_obj_name(benchmark)
  if (!exists(obj_name)) stop("Benchmark não carregado: ", benchmark)
  lst <- get(obj_name)
  if (series_index < 1 || series_index > length(lst)) stop("Índice fora da faixa.")
  df <- lst[[series_index]]
  if (!is.data.frame(df)) stop("Elemento não é data.frame.")
  if (is.na(value_col)) {
    candidates_val <- c("value","metric","y","val","calue")
    hit_val <- candidates_val[candidates_val %in% names(df)]
    if (!length(hit_val)) stop("Não encontrei coluna de valor; informe value_col.")
    value_col <- hit_val[1]
  }
  if (!value_col %in% names(df)) stop("Coluna de valor ausente: ", value_col)
  label_col <- find_label_column(df)
  truth_source <- NULL
  if ("event" %in% names(df)) {
    evt <- df$event
    if (is.logical(evt) || (is.numeric(evt) && all(na.omit(unique(evt)) %in% c(0,1)))) {
      anomalies_truth <- (evt == 1 | evt == TRUE)
      truth_source <- "event"
    } else {
      anomalies_truth <- if (!is.na(label_col)) (df[[label_col]] == 1 | df[[label_col]] == TRUE) else rep(NA, nrow(df))
      truth_source <- if (!is.na(label_col)) "label" else NA
    }
  } else {
    anomalies_truth <- if (!is.na(label_col)) (df[[label_col]] == 1 | df[[label_col]] == TRUE) else rep(NA, nrow(df))
    truth_source <- if (!is.na(label_col)) "label" else NA
  }
  serie <- df[[value_col]]
  list(df = df, serie = serie, benchmark = benchmark, series_index = series_index,
       value_col = value_col, label_col = label_col, truth = anomalies_truth, truth_source = truth_source)
}

compute_confusion <- function(pred_events, truth_events) {
  if (all(is.na(truth_events))) return(c(TP=NA, FP=NA, FN=NA, TN=NA))
  truth_bin <- truth_events == TRUE | truth_events == 1
  pred_bin <- pred_events == TRUE | pred_events == 1
  if (length(pred_bin) != length(truth_bin)) stop("Comprimentos divergentes em pred vs truth.")
  TP <- sum(pred_bin & truth_bin)
  FP <- sum(pred_bin & !truth_bin)
  FN <- sum(!pred_bin & truth_bin)
  TN <- sum(!pred_bin & !truth_bin)
  c(TP=TP, FP=FP, FN=FN, TN=TN)
}

compute_scores <- function(TP, FP, FN, TN) {
  # Coerção defensiva
  TP <- as.numeric(TP); FP <- as.numeric(FP); FN <- as.numeric(FN); TN <- as.numeric(TN)
  if (any(is.na(c(TP,FP,FN,TN)))) {
    return(c(precision=NA, recall=NA, f1=NA, specificity=NA))
  }
  # Precision: se há pelo menos uma predição positiva, calcula; se não há mas existe verdade positiva, define 0.
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else if ((TP + FN) > 0) 0 else NA
  # Recall: se existe verdade positiva, calcula; caso contrário NA.
  recall <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  # F1: somente se precision e recall não NA e >0 no denominador
  f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) 2 * precision * recall / (precision + recall) else if (!is.na(precision) && !is.na(recall) && precision==0 && recall==0) 0 else NA
  # Specificity: se existe TN+FP
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA
  c(precision=precision, recall=recall, f1=f1, specificity=specificity)
}

compute_errors <- function(series, fitted_vals) {
  n <- length(fitted_vals)
  series_adj <- tail(series, n)
  mae <- mean(abs(series_adj - fitted_vals))
  rmse <- sqrt(mean((series_adj - fitted_vals)^2))
  mape <- ifelse(any(series_adj==0), NA, mean(abs((series_adj - fitted_vals)/series_adj)))
  c(mae=mae, rmse=rmse, mape=mape)
}

get_default_filters <- function() {
  list(
    MA = ts_fil_ma(),
    Kalman = ts_fil_kalman(),
    Lowess = ts_fil_lowess(),
    Spline = ts_fil_spline(),
    Winsor = ts_fil_winsor()
  )
}

# Seção de definição de filtros custom (garantia de que todas as funções estejam presentes)
if (!exists("make_custom_filter")) {
  make_custom_filter <- function(name, fun) structure(list(name = name, fun = fun), class = "custom_filter")
  fit.custom_filter <- function(obj, serie) obj
  transform.custom_filter <- function(obj, serie) {
    out <- try(obj$fun(serie), silent = TRUE)
    if (inherits(out, "try-error") || length(out) != length(serie)) serie else out
  }
  custom_filter_identity <- function() make_custom_filter("None", function(x) x)
  custom_filter_ema <- function(alpha = 0.2) {
    make_custom_filter(paste0("EMA", "_a", alpha), function(x) {
      if (!length(x)) return(x)
      y <- numeric(length(x)); y[1] <- x[1]
      for (i in 2:length(x)) y[i] <- alpha * x[i] + (1 - alpha) * y[i-1]
      y
    })
  }
  custom_filter_smooth <- function() {
    make_custom_filter("Smooth", function(x) {
      su <- try(stats::supsmu(seq_along(x), x), silent = TRUE)
      if (!inherits(su, "try-error")) return(su$y)
      k <- 5; filt <- stats::filter(x, rep(1/k, k), sides = 2)
      as.numeric(ifelse(is.na(filt), x, filt))
    })
  }
  custom_filter_hp <- function(lambda = 1600) {
    make_custom_filter(paste0("HP_", lambda), function(x) {
      if (requireNamespace("mFilter", quietly = TRUE)) {
        hp <- try(mFilter::hpfilter(x, freq = lambda, type = "lambda"), silent = TRUE)
        if (!inherits(hp, "try-error")) return(as.numeric(hp$trend))
      }
      x
    })
  }
  custom_filter_fft <- function(cutoff_ratio = 0.1) {
    make_custom_filter(paste0("FFT_lp", cutoff_ratio), function(x) {
      n <- length(x); if (n < 8) return(x)
      k <- max(2, floor(cutoff_ratio * n/2))
      fx <- fft(x)
      idx <- (k+2):(n-k)
      if (length(idx) > 0) fx[idx] <- 0+0i
      Re(fft(fx, inverse = TRUE)/n)
    })
  }
  custom_filter_seasonal_adj <- function(freq = NULL) {
    make_custom_filter(paste0("Seasonal_Adj", ifelse(is.null(freq), "", paste0("_f", freq))), function(x) {
      n <- length(x)
      if (is.null(freq)) {
        cand <- c(24, 7, 12)
        freq <- cand[which(cand < n)[1]]
        if (is.na(freq)) return(x)
      }
      if (freq <= 1 || freq*3 > n) return(x)
      tsx <- stats::ts(x, frequency = freq)
      decomp <- try(stats::stl(tsx, s.window = "periodic", robust = TRUE), silent = TRUE)
      if (inherits(decomp, "try-error")) return(x)
      adj <- decomp$time.series[, "trend"] + decomp$time.series[, "remainder"]
      as.numeric(adj)
    })
  }
  custom_filter_recursive <- function(alpha = 0.5) {
    make_custom_filter(paste0("Recursive_Filter_a", alpha), function(x) {
      if (!length(x)) return(x)
      y <- numeric(length(x)); y[1] <- x[1]
      for (i in 2:length(x)) y[i] <- alpha * y[i-1] + (1 - alpha) * x[i]
      y
    })
  }
  custom_filter_emd <- function() {
    make_custom_filter("EMD", function(x) {
      if (requireNamespace("EMD", quietly = TRUE)) {
        emd_res <- try(EMD::emd(x), silent = TRUE)
        if (!inherits(emd_res, "try-error") && !is.null(emd_res$residue)) return(as.numeric(emd_res$residue))
      }
      x
    })
  }
  custom_filter_remd <- function() {
    make_custom_filter("REMD", function(x) {
      if (requireNamespace("Rlibeemd", quietly = TRUE)) {
        em <- try(Rlibeemd::eemd(x), silent = TRUE)
        if (!inherits(em, "try-error")) {
          rec <- rowMeans(em)
          return(as.numeric(rec))
        }
      }
      x
    })
  }
  custom_filter_wavelet <- function() {
    make_custom_filter("Wavelet", function(x) {
      if (requireNamespace("wavelets", quietly = TRUE)) {
        wt <- try(wavelets::dwt(x, boundary = "reflection"), silent = TRUE)
        if (!inherits(wt, "try-error")) {
          for (i in seq_along(wt@W)) wt@W[[i]][] <- 0
          rec <- try(wavelets::idwt(wt), silent = TRUE)
          if (!inherits(rec, "try-error")) return(as.numeric(rec))
        }
      }
      x
    })
  }
}

# ================= Filtros customizados adicionais =====================
# Estrutura para permitir filtros além dos fornecidos pelo daltoolbox.
# Cada filtro custom precisa apenas de uma função fun(x) -> série filtrada.
make_custom_filter <- function(name, fun) structure(list(name = name, fun = fun), class = "custom_filter")

fit.custom_filter <- function(obj, serie) { obj }  # nada a ajustar
transform.custom_filter <- function(obj, serie) { 
  out <- try(obj$fun(serie), silent = TRUE)
  if (inherits(out, "try-error") || length(out) != length(serie)) return(serie) else return(out)
}

# Construtores de filtros customizados (fallback seguro caso pacotes ausentes)
custom_filter_identity <- function() make_custom_filter("None", function(x) x)
custom_filter_ema <- function(alpha = 0.2) {
  make_custom_filter(paste0("EMA", "_a", alpha), function(x) {
    if (!length(x)) return(x)
    y <- numeric(length(x)); y[1] <- x[1]
    for (i in 2:length(x)) y[i] <- alpha * x[i] + (1 - alpha) * y[i-1]
    y
  })
}
custom_filter_smooth <- function() {
  make_custom_filter("Smooth", function(x) {
    if (requireNamespace("stats", quietly = TRUE)) {
      # Tenta supsmu; fallback média móvel janela 5
      su <- try(stats::supsmu(seq_along(x), x), silent = TRUE)
      if (!inherits(su, "try-error")) return(su$y)
    }
    # fallback média móvel simples
    k <- 5; filt <- stats::filter(x, rep(1/k, k), sides = 2)
    as.numeric(ifelse(is.na(filt), x, filt))
  })
}
custom_filter_hp <- function(lambda = 1600) {
  make_custom_filter(paste0("HP_", lambda), function(x) {
    if (requireNamespace("mFilter", quietly = TRUE)) {
      hp <- try(mFilter::hpfilter(x, freq = lambda, type = "lambda"), silent = TRUE)
      if (!inherits(hp, "try-error")) return(as.numeric(hp$trend))
    }
    x  # fallback
  })
}
custom_filter_fft <- function(cutoff_ratio = 0.1) {
  make_custom_filter(paste0("FFT_lp", cutoff_ratio), function(x) {
    n <- length(x); if (n < 8) return(x)
    k <- max(2, floor(cutoff_ratio * n/2))
    fx <- fft(x)
    # Zera altas frequências (mantém k componentes de baixa frequência em cada lado)
    idx <- (k+2):(n-k)
    fx[idx] <- 0+0i
    Re(fft(fx, inverse = TRUE)/n)
  })
}
custom_filter_seasonal_adj <- function(freq = NULL) {
  make_custom_filter(paste0("Seasonal_Adj", ifelse(is.null(freq), "", paste0("_f", freq))), function(x) {
    n <- length(x)
    if (is.null(freq)) {
      # heurística simples: tenta 24, 7, 12
      cand <- c(24, 7, 12)
      freq <- cand[which(cand < n)[1]]
      if (is.na(freq)) return(x)
    }
    if (freq <= 1 || freq*3 > n) return(x)
    tsx <- stats::ts(x, frequency = freq)
    decomp <- try(stats::stl(tsx, s.window = "periodic", robust = TRUE), silent = TRUE)
    if (inherits(decomp, "try-error")) return(x)
    adj <- decomp$time.series[, "trend"] + decomp$time.series[, "remainder"]
    as.numeric(adj)
  })
}
custom_filter_recursive <- function(alpha = 0.5) {
  make_custom_filter(paste0("Recursive_Filter_a", alpha), function(x) {
    if (!length(x)) return(x)
    y <- numeric(length(x)); y[1] <- x[1]
    for (i in 2:length(x)) y[i] <- alpha * y[i-1] + (1 - alpha) * x[i]
    y
  })
}
# Placeholders para EMD / REMD / Wavelet com fallback seguro
custom_filter_emd <- function() {
  make_custom_filter("EMD", function(x) {
    if (requireNamespace("EMD", quietly = TRUE)) {
      emd_res <- try(EMD::emd(x), silent = TRUE)
      if (!inherits(emd_res, "try-error") && !is.null(emd_res$residue)) return(as.numeric(emd_res$residue))
    }
    x
  })
}
custom_filter_remd <- function() {
  make_custom_filter("REMD", function(x) {
    if (requireNamespace("Rlibeemd", quietly = TRUE)) {
      em <- try(Rlibeemd::eemd(x), silent = TRUE)
      if (!inherits(em, "try-error")) {
        # Reconstrução simplista: média das componentes
        rec <- rowMeans(em)
        return(as.numeric(rec))
      }
    }
    x
  })
}
custom_filter_wavelet <- function() {
  make_custom_filter("Wavelet", function(x) {
    if (requireNamespace("wavelets", quietly = TRUE)) {
      wt <- try(wavelets::dwt(x, boundary = "reflection"), silent = TRUE)
      if (!inherits(wt, "try-error")) {
        # Threshold simples: zera detalhes nível alto
        for (i in seq_along(wt@W)) wt@W[[i]][] <- 0
        rec <- try(wavelets::idwt(wt), silent = TRUE)
        if (!inherits(rec, "try-error")) return(as.numeric(rec))
      }
    }
    x
  })
}

get_extended_filters <- function(include_default = TRUE) {
  base <- if (include_default) get_default_filters() else list()
  # Adiciona todos os extras conforme TCC
  extras <- list(
    None = custom_filter_identity(),
    EMA = custom_filter_ema(),
    EMD = custom_filter_emd(),
    REMD = custom_filter_remd(),
    Wavelet = custom_filter_wavelet(),
    Smooth = custom_filter_smooth(),
    HP = custom_filter_hp(),
    FFT = custom_filter_fft(),
    Seasonal_Adj = custom_filter_seasonal_adj(),
    Recursive_Filter = custom_filter_recursive()
  )
  # Merge garantindo nomes únicos (prioriza explicitamente extras sobre base em conflito exceto filtros já existentes)
  dup <- intersect(names(base), names(extras))
  if (length(dup)) {
    # Mantém base para filtros originais (MA, Kalman, Lowess, Spline, Winsor)
    extras <- extras[setdiff(names(extras), dup)]
  }
  c(base, extras)
}

filter_category <- function(fname) {
  # Normaliza nome removendo sufixos parametrizados
  base <- fname
  base <- sub("_a[0-9.]+$", "", base)         # EMA_a0.2 / Recursive_Filter_a0.5
  base <- sub("_lp[0-9.]+$", "", base)        # FFT_lp0.1
  base <- sub("_f[0-9.]+$", "", base)         # Seasonal_Adj_f24
  base <- sub("_\\\u003d?[0-9]+$", "", base)  # HP_1600 (fallback)
  if (base %in% c("MA","EMA")) return("Media")
  if (base %in% c("EMD","REMD","Wavelet")) return("Decomposicao")
  if (base %in% c("Spline","Smooth","LOWESS")) return("Suavizacao")
  if (base %in% c("Kalman","HP")) return("Adaptativo")
  if (base %in% c("Winsor","Seasonal_Adj","FFT","Recursive_Filter")) return("Robusto")
  if (base %in% c("None")) return("Baseline")
  "Outro"
}

run_filters_yahoo <- function(benchmark = "A1", series_index = 1, value_col = NA,
                              filters = get_extended_filters(), plot = TRUE,
                              save_plots = FALSE, plots_dir = "figures_yahoo_filters",
                              return_models = FALSE,
                              truth_point_col = "orange", truth_point_pch = 21, # círculo vazado
                              detected_point_col = "green", detected_point_pch = 4) { # x
  prep <- prepare_yahoo_series(benchmark, series_index, value_col)
  serie <- na.omit(prep$serie)
  truth <- prep$truth
  label_col <- prep$label_col
  truth_source <- prep$truth_source
  if (save_plots) dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  metrics_list <- list(); results <- list()
  truth_idx <- if (!all(is.na(truth))) which(truth) else integer(0)
  for (fname in names(filters)) {
    filtro_obj <- filters[[fname]]
    if (plot) {
      op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
      plot(serie, type = "l", col = "#2c3e50", main = paste("Yahoo", benchmark, "Série", series_index, "- Filtro", fname),
           xlab = "Índice", ylab = prep$value_col)
      # Desenhar gabarito (ground truth) antes das detecções
      if (length(truth_idx)) points(truth_idx, serie[truth_idx], col = truth_point_col, bg = NA, pch = truth_point_pch)
    }
    filtro_fit <- try(fit(filtro_obj, serie), silent = TRUE)
    if (inherits(filtro_fit, "try-error")) {
      warning("Falha ao ajustar filtro ", fname, ": ", filtro_fit)
      filtered_series <- serie
      filtro_fit <- NULL
    } else {
      filtered_series <- try(transform(filtro_fit, serie), silent = TRUE)
      if (inherits(filtered_series, "try-error") || length(filtered_series) != length(serie)) filtered_series <- serie
    }
    if (plot) lines(filtered_series, col = "red", lty = 2)
    model <- arima_filter(filtro_fit)
    model <- fit(model, serie)
    detection <- detect(model, serie)
    pred_events <- detection$event
    detected_idx <- if (!is.null(pred_events) && length(pred_events)==length(serie)) which(pred_events) else integer(0)
    if (plot && length(detected_idx)) {
      points(detected_idx, serie[detected_idx], col = detected_point_col, pch = detected_point_pch)
      legend_entries <- c("Original", "Filtrada")
      legend_cols <- c("#2c3e50", "red")
      legend_lty <- c(1,2)
      legend_pch <- c(NA, NA)
      # Adicionar gabarito se existir
      if (length(truth_idx)) {
        legend_entries <- c(legend_entries, "Gabarito")
        legend_cols <- c(legend_cols, truth_point_col)
        legend_lty <- c(legend_lty, NA)
        legend_pch <- c(legend_pch, truth_point_pch)
      }
      # Adicionar detectado
      legend_entries <- c(legend_entries, "Detectado")
      legend_cols <- c(legend_cols, detected_point_col)
      legend_lty <- c(legend_lty, NA)
      legend_pch <- c(legend_pch, detected_point_pch)
      legend("topright", legend = legend_entries, col = legend_cols, lty = legend_lty, pch = legend_pch, bty = "n", cex = 0.8)
    } else if (plot) {
      # Sem detecções, ainda mostrar legenda incluindo gabarito se houver
      legend_entries <- c("Original", "Filtrada")
      legend_cols <- c("#2c3e50", "red")
      legend_lty <- c(1,2)
      legend_pch <- c(NA, NA)
      if (length(truth_idx)) {
        legend_entries <- c(legend_entries, "Gabarito")
        legend_cols <- c(legend_cols, truth_point_col)
        legend_lty <- c(legend_lty, NA)
        legend_pch <- c(legend_pch, truth_point_pch)
      }
      legend("topright", legend = legend_entries, col = legend_cols, lty = legend_lty, pch = legend_pch, bty = "n", cex = 0.8)
    }
    conf <- compute_confusion(pred_events, truth)
    conf <- c(TP=as.numeric(conf["TP"]), FP=as.numeric(conf["FP"]), FN=as.numeric(conf["FN"]), TN=as.numeric(conf["TN"]))
    scores <- compute_scores(conf["TP"], conf["FP"], conf["FN"], conf["TN"])
    fitted_vals <- fitted(model$model)
    errs <- compute_errors(serie, fitted_vals)
    detection_count <- if (!is.null(pred_events)) sum(pred_events) else NA
    detection_rate <- if (!is.null(pred_events)) detection_count / length(serie) else NA
    metrics_list[[fname]] <- data.frame(
      benchmark = benchmark,
      series = series_index,
      filtro = fname,
      categoria = filter_category(fname),
      truth_source = truth_source,
      TP = conf["TP"], FP = conf["FP"], FN = conf["FN"], TN = conf["TN"],
      precision = scores["precision"], recall = scores["recall"], f1 = scores["f1"], specificity = scores["specificity"],
      detection_count = detection_count, detection_rate = detection_rate,
      mae = errs["mae"], rmse = errs["rmse"], mape = errs["mape"],
      arima_p = model$p, arima_d = model$d, arima_q = model$q,
      stringsAsFactors = FALSE
    )
    results[[fname]] <- list(model = model, detection = detection, filtered_series = filtered_series)
  }
  metrics_df <- dplyr::bind_rows(metrics_list)
  rownames(metrics_df) <- NULL
  out <- list(metrics = metrics_df, results = if (return_models) results else NULL,
              meta = list(benchmark = benchmark, series_index = series_index, value_col = prep$value_col,
                          label_col = label_col, truth_source = truth_source))
  class(out) <- c("yahoo_filter_run", class(out))
  out
}

list_labeled_series <- function(benchmark = "A1") {
  load_yahoo_benchmarks(benchmark)
  obj_name <- yahoo_obj_name(benchmark)
  lst <- get(obj_name)
  tibble::tibble(
    series = seq_along(lst),
    has_event = sapply(lst, function(df) "event" %in% names(df) && is.logical(df$event) || (is.numeric(df$event) && all(na.omit(unique(df$event)) %in% c(0,1)) )),
    label_col = sapply(lst, find_label_column)
  ) %>% mutate(has_label = !is.na(label_col))
}

run_filters_yahoo_all_labeled <- function(benchmark = "A1", value_col = NA, filters = get_extended_filters(), plot = FALSE) {
  info <- list_labeled_series(benchmark)
  sel <- which(info$has_event | info$has_label)
  if (!length(sel)) stop("Nenhuma série com ground truth detectado em ", benchmark)
  rows <- list(); all_results <- list()
  for (s in sel) {
    res <- run_filters_yahoo(benchmark = benchmark, series_index = s, value_col = value_col, filters = filters, plot = plot)
    mdf <- res$metrics; mdf$series <- s
    rows[[length(rows)+1]] <- mdf
    all_results[[paste0("series_", s)]] <- res
  }
  agg <- dplyr::bind_rows(rows)
  list(aggregate_metrics = agg, per_series = all_results, labeled_info = info)
}

run_filters_yahoo_all <- function(benchmark = "A1", value_col = NA, filters = get_extended_filters(), plot = FALSE) {
  load_yahoo_benchmarks(benchmark)
  obj_name <- yahoo_obj_name(benchmark)
  if (!exists(obj_name)) stop("Benchmark não carregado: ", benchmark)
  lst <- get(obj_name)
  all_rows <- list(); all_results <- list()
  for (i in seq_along(lst)) {
    res <- try(run_filters_yahoo(benchmark = benchmark, series_index = i, value_col = value_col,
                                 filters = filters, plot = plot), silent = TRUE)
    if (inherits(res, "try-error")) {
      warning("Falha ao processar série ", i, ": ", res)
      next
    }
    mdf <- res$metrics; mdf$series <- i
    all_rows[[length(all_rows)+1]] <- mdf
    all_results[[paste0("series_", i)]] <- res
  }
  agg <- dplyr::bind_rows(all_rows)
  # Adicionar flag se há ground truth
  agg$has_truth <- !is.na(agg$truth_source) & agg$truth_source != ""
  list(aggregate_metrics = agg, per_series = all_results)
}

save_filter_metrics <- function(run_res, file = "metrics_yahoo_filters.csv") {
  stopifnot(inherits(run_res, "yahoo_filter_run"))
  utils::write.csv(run_res$metrics, file, row.names = FALSE)
  message("Métricas salvas em ", file)
  invisible(file)
}

save_aggregate_metrics <- function(agg_res, file = "metrics_yahoo_filters_all.csv") {
  utils::write.csv(agg_res$aggregate_metrics, file, row.names = FALSE)
  message("Métricas agregadas salvas em ", file)
  invisible(file)
}

count_events_true <- function(benchmark = "A1") {
  load_yahoo_benchmarks(benchmark)
  obj_name <- yahoo_obj_name(benchmark)
  if (!exists(obj_name)) stop("Benchmark não carregado: ", benchmark)
  lst <- get(obj_name)
  counts <- lapply(seq_along(lst), function(i) {
    df <- lst[[i]]
    if (!is.data.frame(df) || !("event" %in% names(df))) return(data.frame(series=i, event_true=NA_integer_))
    ev <- df$event
    # Normalizar para lógico
    ev_bin <- if (is.logical(ev)) ev else if (is.numeric(ev)) ev == 1 else rep(FALSE, length(ev))
    data.frame(series=i, event_true=sum(ev_bin, na.rm=TRUE))
  })
  res <- dplyr::bind_rows(counts)
  res$total_event_true <- sum(res$event_true, na.rm = TRUE)
  res
}

export_metrics_per_series <- function(aggregate_metrics, benchmark = "A1", dir = "metrics_por_serie") {
  # aggregate_metrics: data.frame retornado em all_a1$aggregate_metrics
  stopifnot(is.data.frame(aggregate_metrics))
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  series_ids <- sort(unique(aggregate_metrics$series))
  files <- c()
  for (s in series_ids) {
    df_s <- subset(aggregate_metrics, series == s)
    fname <- file.path(dir, sprintf("%s_series_%02d.csv", benchmark, s))
    utils::write.csv(df_s, fname, row.names = FALSE)
    files <- c(files, fname)
  }
  message("Exportadas ", length(files), " séries em ", dir)
  invisible(files)
}

print_metrics_all_series <- function(aggregate_metrics, series_limit = NULL) {
  stopifnot(is.data.frame(aggregate_metrics))
  series_ids <- sort(unique(aggregate_metrics$series))
  if (!is.null(series_limit)) series_ids <- head(series_ids, series_limit)
  for (s in series_ids) {
    cat("\n================ Série", s, "================\n")
    print(subset(aggregate_metrics, series == s))
  }
  invisible(TRUE)
}

# ==============================================================
# Execução multi-benchmark (A1..A4) e sumarizações
# ==============================================================

run_filters_yahoo_all_benchmarks <- function(benchmarks = c("A1","A2","A3","A4"),
                                             value_col = NA,
                                             filters = get_extended_filters(),
                                             plot = FALSE,
                                             out_dir = "metrics_benchmarks",
                                             save_per_benchmark = TRUE,
                                             combined_filename = "metrics_A1_A4_combinado.csv",
                                             summarize = TRUE) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  all_list <- list()
  for (b in benchmarks) {
    message("[Yahoo Filtros] Rodando benchmark ", b)
    res <- run_filters_yahoo_all(benchmark = b, value_col = value_col, filters = filters, plot = plot)
    agg <- res$aggregate_metrics
    if (save_per_benchmark) {
      fn <- file.path(out_dir, sprintf("metrics_%s.csv", b))
      utils::write.csv(agg, fn, row.names = FALSE)
      message("  -> CSV salvo: ", fn, " (", nrow(agg), " linhas)")
    }
    all_list[[b]] <- agg
  }
  combined <- dplyr::bind_rows(all_list)
  comb_path <- file.path(out_dir, combined_filename)
  utils::write.csv(combined, comb_path, row.names = FALSE)
  message("[Yahoo Filtros] CSV combinado salvo em ", comb_path, " (", nrow(combined), " linhas)")
  summaries <- list()
  if (summarize) {
    summaries$by_category <- summarize_metrics_by_category(combined)
    summaries$by_filter <- summarize_metrics_by_filter(combined)
    # Salvar sumarizações
    utils::write.csv(summaries$by_category, file.path(out_dir, "summary_por_categoria.csv"), row.names = FALSE)
    utils::write.csv(summaries$by_filter, file.path(out_dir, "summary_por_filtro.csv"), row.names = FALSE)
    message("[Yahoo Filtros] Sumarizações salvas (categoria e filtro)")
  }
  list(combined = combined, per_benchmark = all_list, summaries = summaries)
}

summarize_metrics_by_category <- function(aggregate_metrics) {
  req_cols <- c("benchmark","categoria","precision","recall","f1","specificity","mae","rmse","mape","detection_rate","series")
  missing <- setdiff(req_cols, names(aggregate_metrics))
  if (length(missing)) stop("Colunas ausentes para sumarização: ", paste(missing, collapse=","))
  aggregate_metrics %>%
    dplyr::group_by(benchmark, categoria) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_series = dplyr::n_distinct(series),
      precision = mean(precision, na.rm = TRUE),
      recall = mean(recall, na.rm = TRUE),
      f1 = mean(f1, na.rm = TRUE),
      specificity = mean(specificity, na.rm = TRUE),
      mae = mean(mae, na.rm = TRUE),
      rmse = mean(rmse, na.rm = TRUE),
      mape = mean(mape, na.rm = TRUE),
      detection_rate = mean(detection_rate, na.rm = TRUE),
      .groups = "drop"
    ) %>% dplyr::arrange(benchmark, categoria)
}

summarize_metrics_by_filter <- function(aggregate_metrics) {
  req_cols <- c("benchmark","filtro","precision","recall","f1","specificity","mae","rmse","mape","detection_rate","series")
  missing <- setdiff(req_cols, names(aggregate_metrics))
  if (length(missing)) stop("Colunas ausentes para sumarização: ", paste(missing, collapse=","))
  aggregate_metrics %>%
    dplyr::group_by(benchmark, filtro) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_series = dplyr::n_distinct(series),
      precision = mean(precision, na.rm = TRUE),
      recall = mean(recall, na.rm = TRUE),
      f1 = mean(f1, na.rm = TRUE),
      specificity = mean(specificity, na.rm = TRUE),
      mae = mean(mae, na.rm = TRUE),
      rmse = mean(rmse, na.rm = TRUE),
      mape = mean(mape, na.rm = TRUE),
      detection_rate = mean(detection_rate, na.rm = TRUE),
      .groups = "drop"
    ) %>% dplyr::arrange(benchmark, filtro)
}


AUTO_RUN_FILTERS_YAHOO <- TRUE
AUTO_RUN_SERIES_INDEX <- 1
AUTO_RUN_BENCHMARK <- "A1"
AUTO_SAVE_PLOTS <- FALSE
if (AUTO_RUN_FILTERS_YAHOO) {
  message("[Yahoo Filtros] Executando análise de filtros em ", AUTO_RUN_BENCHMARK, " série ", AUTO_RUN_SERIES_INDEX)
  run_res <- run_filters_yahoo(benchmark = AUTO_RUN_BENCHMARK, series_index = AUTO_RUN_SERIES_INDEX,
                               plot = TRUE, save_plots = AUTO_SAVE_PLOTS)
  print(run_res$metrics)
}

AUTO_RUN_ALL_A1 <- TRUE  # defina FALSE se não quiser rodar tudo automaticamente
AUTO_RUN_ALL_SAVE_CSV <- TRUE
AUTO_RUN_ALL_PLOT <- FALSE  # cuidado: TRUE gera 67* n_filtros plots
AUTO_RUN_ALL_OUTPUT_CSV <- "metrics_A1_todas_series.csv"

if (AUTO_RUN_ALL_A1) {
  message("[Yahoo Filtros] Rodando TODAS as séries do benchmark A1 (1..67)")
  all_a1 <- run_filters_yahoo_all("A1", plot = AUTO_RUN_ALL_PLOT)
  agg <- all_a1$aggregate_metrics
  if (AUTO_RUN_ALL_SAVE_CSV) {
    utils::write.csv(agg, AUTO_RUN_ALL_OUTPUT_CSV, row.names = FALSE)
    message("Arquivo agregado salvo em ", AUTO_RUN_ALL_OUTPUT_CSV)
  }
  # Resumo rápido
  message("Total de linhas (filtro x série): ", nrow(agg))
  message("Séries com ground truth: ", length(unique(agg$series[agg$has_truth])))
  message("Primeiras linhas:")
  print(head(agg))
}

# --------------------------------------------------------------
# Auto-run multi benchmarks (A1..A4) - desativado por padrão
# --------------------------------------------------------------
AUTO_RUN_ALL_BENCHMARKS <- FALSE
AUTO_RUN_ALL_BENCH_FILTERS <- TRUE
AUTO_RUN_ALL_BENCH_OUTDIR <- "metrics_benchmarks"
if (AUTO_RUN_ALL_BENCHMARKS) {
  multi_res <- run_filters_yahoo_all_benchmarks(
    benchmarks = c("A1","A2","A3","A4"),
    filters = if (AUTO_RUN_ALL_BENCH_FILTERS) get_extended_filters() else get_default_filters(),
    plot = FALSE,
    out_dir = AUTO_RUN_ALL_BENCH_OUTDIR,
    summarize = TRUE
  )
  message("[Yahoo Filtros] Multi-benchmark concluído. Linhas combinadas: ", nrow(multi_res$combined))
  message("Arquivos salvos em diretório: ", AUTO_RUN_ALL_BENCH_OUTDIR)
}

# Instruções rápidas:
# list_labeled_series("A1")
# run_res <- run_filters_yahoo("A1", 1)
# save_filter_metrics(run_res, "saida_metrics_A1_1.csv")
# agg <- run_filters_yahoo_all_labeled("A1")
# save_aggregate_metrics(agg, "saida_metrics_A1_all.csv")




