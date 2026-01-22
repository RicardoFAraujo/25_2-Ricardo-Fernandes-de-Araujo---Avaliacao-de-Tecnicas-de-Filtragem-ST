###############################################################################
# TCC - Análise Comparativa de Filtros com ARIMA (VERSÃO CORRIGIDA)
# Múltiplas Variáveis + Tabela Agregada por Dataset
###############################################################################

# Carregar bibliotecas
library(united)
library(daltoolbox) 
library(harbinger) 
library(tspredit)
library(ggplot2)
library(reshape2)
library(corrplot)

source("arima_filter.R")

# Carregar o dataset
data(gecco)

# Selecionar os dados (24 horas de dados)
gecco <- gecco$mult[16500:17939, ]  # 1440 amostras (60*24)

# MÚLTIPLAS VARIÁVEIS para análise (EXPLÍCITO, conforme solicitado)
vars <- c("tp", "cl", "ph", "redox", "leit", "trueb", "cl_2", "fm", "fm_2")
missing_vars <- setdiff(vars, names(gecco))
if (length(missing_vars) > 0) {
  stop("Variáveis faltando no gecco: ", paste(missing_vars, collapse = ", "))
}

series_list <- list(
  tp = as.numeric(gecco$tp),
  cl = as.numeric(gecco$cl),
  ph = as.numeric(gecco$ph),
  redox = as.numeric(gecco$redox),
  leit = as.numeric(gecco$leit),
  trueb = as.numeric(gecco$trueb),
  cl_2 = as.numeric(gecco$cl_2),
  fm = as.numeric(gecco$fm),
  fm_2 = as.numeric(gecco$fm_2)
)

# Filtros do tspredit
filters <- list(
  ts_fil_recursive(filter = 0.5), ts_fil_ema(), ts_fil_emd(), ts_fil_fft(), ts_fil_hp(),
  ts_fil_kalman(), ts_fil_lowess(), ts_fil_ma(), ts_fil_none(),
  ts_fil_qes(), ts_fil_remd(), ts_fil_seas_adj(), ts_fil_ses(gamma = FALSE),
  ts_fil_smooth(), ts_fil_spline(), ts_fil_wavelet(), ts_fil_winsor()
)

filter_names <- c(
  "Recursive_Filter","EMA", "EMD", "FFT", "HP", "Kalman", "Lowess", "MA", "None",
  "QES", "REMD", "Seasonal_Adj", "SES", "Smooth", "Spline", "Wavelet", "Winsor"
)

# Inicializar tabela de resultados DETALHADA (por variável)
detailed_results <- data.frame(
  Dataset = character(),
  Variable = character(),
  Filter = character(),
  TP = numeric(),
  FN = numeric(), 
  FP = numeric(),
  TN = numeric(),
  Precision = numeric(),
  Recall = numeric(),
  F1_Score = numeric(),
  Accuracy = numeric(),
  MAE = numeric(),
  RMSE = numeric(),
  stringsAsFactors = FALSE
)

cat("=== ANÁLISE CONSOLIDADA PARA MÚLTIPLAS VARIÁVEIS ===\n")
cat("Dataset: GECCO\n")
cat("Variáveis:", paste(names(series_list), collapse = ", "), "\n")
cat("Filtros:", length(filters), "\n\n")

# Loop principal por variável
for (variable_name in names(series_list)) {
  series <- series_list[[variable_name]]
  
  cat("\n", rep("=", 60), "\n")
  cat("PROCESSANDO VARIÁVEL:", toupper(variable_name), "\n")
  cat(rep("=", 60), "\n")
  
  # Remover NAs
  valid_idx <- which(!is.na(series))
  series <- series[valid_idx]
  aligned_event <- gecco$event[valid_idx]
  
  cat("Tamanho da série:", length(series), "observações\n")
  cat("Eventos reais:", sum(aligned_event), "\n\n")
  
  # Loop por filtro
  for (i in seq_along(filters)) {
    filter <- filters[[i]]
    filter_name <- filter_names[i]
    
    cat("Processando", variable_name, "+", filter_name, "...")
    
    # Aplicar filtro
    filter <- tryCatch(fit(filter, series), error = function(e) NULL)
    if (is.null(filter)) {
      cat(" ERRO no filtro\n")
      next
    }
    
    # ARIMA com filtro
    model <- arima_filter(filter)
    model <- tryCatch(fit(model, series), error = function(e) NULL)
    if (is.null(model)) {
      cat(" ERRO no ARIMA\n")
      next
    }
    
    # Detecção
    detection <- tryCatch(detect(model, series), error = function(e) NULL)
    if (is.null(detection)) {
      cat(" ERRO na detecção\n")
      next
    }
    
    # Calcular métricas (versão manual robusta)
    if (length(aligned_event) == length(detection$event)) {
      
      detected <- as.logical(detection$event)
      real <- as.logical(aligned_event)
      
      # Matriz de confusão manual
      TP <- sum(detected & real, na.rm = TRUE)
      FN <- sum(!detected & real, na.rm = TRUE)
      FP <- sum(detected & !real, na.rm = TRUE)
      TN <- sum(!detected & !real, na.rm = TRUE)
      
      # Métricas
      precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
      recall <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
      f1_score <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), 0)
      accuracy <- ifelse((TP + FN + FP + TN) > 0, (TP + TN) / (TP + FN + FP + TN), 0)
      
      # Erros de predição
      mae <- NA
      rmse <- NA
      tryCatch({
        if(!is.null(model$model) && length(fitted(model$model)) > 0) {
          n <- length(fitted(model$model))
          if(n <= length(series)) {
            series_adj <- tail(series, n)
            fitted_values <- fitted(model$model)
            if(length(series_adj) == length(fitted_values)) {
              mae <- mean(abs(series_adj - fitted_values), na.rm = TRUE)
              rmse <- sqrt(mean((series_adj - fitted_values)^2, na.rm = TRUE))
            }
          }
        }
      }, error = function(e) {})
      
      # Adicionar aos resultados detalhados
      detailed_results <- rbind(detailed_results, data.frame(
        Dataset = "GECCO",
        Variable = variable_name,
        Filter = filter_name,
        TP = TP, FN = FN, FP = FP, TN = TN,
        Precision = round(precision, 4),
        Recall = round(recall, 4),
        F1_Score = round(f1_score, 4),
        Accuracy = round(accuracy, 4),
        MAE = ifelse(is.na(mae), 0, round(mae, 4)),
        RMSE = ifelse(is.na(rmse), 0, round(rmse, 4)),
        stringsAsFactors = FALSE
      ))
      
      cat(" OK | F1:", round(f1_score, 3), "| Prec:", round(precision, 3), 
          "| Rec:", round(recall, 3), "\n")
      
    } else {
      cat(" ERRO: Tamanhos diferentes\n")
    }
  }
}

###############################################################################
# VERIFICAÇÃO DOS DADOS COLETADOS
###############################################################################

cat("\n", rep("=", 80), "\n")
cat("VERIFICAÇÃO DOS DADOS COLETADOS\n")
cat(rep("=", 80), "\n")

# Verificar se temos dados
if (nrow(detailed_results) == 0) {
  cat("ERRO: Nenhum resultado foi coletado!\n")
  cat("Verifique se os filtros estão funcionando corretamente.\n")
  stop("Análise interrompida - sem dados para processar.")
}

cat("Total de resultados coletados:", nrow(detailed_results), "\n")
cat("Variáveis processadas:", paste(unique(detailed_results$Variable), collapse = ", "), "\n")
cat("Filtros processados:", length(unique(detailed_results$Filter)), "\n\n")

# Salvar resultados detalhados
write.csv(detailed_results, "resultados_detalhados_por_variavel.csv", row.names = FALSE)
cat("Resultados detalhados salvos em: resultados_detalhados_por_variavel.csv\n\n")

# Mostrar amostra dos dados
cat("AMOSTRA DOS DADOS COLETADOS:\n")
print(head(detailed_results[, c("Variable", "Filter", "Precision", "Recall", "F1_Score", "Accuracy")]))

###############################################################################
# CONSOLIDAÇÃO DOS RESULTADOS POR DATASET (VERSÃO CORRIGIDA)
###############################################################################

cat("\n", rep("=", 80), "\n")
cat("CONSOLIDAÇÃO DOS RESULTADOS POR DATASET\n")
cat(rep("=", 80), "\n")

# Agregar métricas por filtro (média das variáveis)
cat("Calculando métricas agregadas por filtro...\n")

aggregated_results <- aggregate(
  cbind(Precision, Recall, F1_Score, Accuracy, MAE, RMSE) ~ Filter, 
  data = detailed_results, 
  FUN = function(x) round(mean(x, na.rm = TRUE), 4)
)

cat("Agregação concluída!\n")
cat("Filtros agregados:", nrow(aggregated_results), "\n\n")

# Verificar se a agregação funcionou
if (nrow(aggregated_results) == 0) {
  cat("ERRO: Falha na agregação dos dados!\n")
  stop("Análise interrompida - problema na agregação.")
}

print(head(aggregated_results))

###############################################################################
# TABELA CONSOLIDADA CONFORME SOLICITADO PELO PROFESSOR (VERSÃO CORRIGIDA)
###############################################################################

cat("\n", rep("=", 60), "\n")
cat("GERANDO TABELA CONSOLIDADA...\n")
cat(rep("=", 60), "\n")

# Criar tabela consolidada de forma mais robusta
create_consolidated_table_safe <- function(aggregated_data) {
  
  # Verificar se temos dados
  if (nrow(aggregated_data) == 0) {
    cat("ERRO: Dados agregados estão vazios!\n")
    return(NULL)
  }
  
  # Selecionar apenas as métricas principais
  metrics_cols <- c("Precision", "Recall", "Accuracy", "F1_Score")
  
  # Verificar se as colunas existem
  missing_cols <- setdiff(metrics_cols, names(aggregated_data))
  if (length(missing_cols) > 0) {
    cat("ERRO: Colunas faltando:", paste(missing_cols, collapse = ", "), "\n")
    return(NULL)
  }
  
  # Criar matriz de métricas
  metrics_matrix <- as.matrix(aggregated_data[, metrics_cols])
  rownames(metrics_matrix) <- aggregated_data$Filter
  
  # Transpor para ter métricas nas linhas
  consolidated_matrix <- t(metrics_matrix)
  
  # Converter para data.frame
  consolidated_df <- as.data.frame(consolidated_matrix)
  
  # Adicionar coluna de métricas como primeira coluna
  consolidated_df <- data.frame(
    Metrica = c("Precisão (PPV)", "Recall (TPR)", "Acurácia", "F1-Score"),
    consolidated_df,
    stringsAsFactors = FALSE
  )
  
  return(consolidated_df)
}

# Gerar tabela consolidada
tabela_consolidada_gecco <- create_consolidated_table_safe(aggregated_results)

# Verificar se a tabela foi criada
if (is.null(tabela_consolidada_gecco)) {
  cat("ERRO: Falha ao criar tabela consolidada!\n")
  stop("Análise interrompida.")
}

cat("=== TABELA CONSOLIDADA - DATASET GECCO ===\n")
cat("(Métricas calculadas por variável e agregadas por média aritmética)\n\n")

# Exibir tabela
print(tabela_consolidada_gecco)

# Salvar tabela consolidada
write.csv(tabela_consolidada_gecco, "tabela_consolidada_gecco.csv", row.names = FALSE)
cat("\nTabela consolidada salva em: tabela_consolidada_gecco.csv\n")

###############################################################################
# ANÁLISE ESTATÍSTICA DA TABELA CONSOLIDADA (VERSÃO CORRIGIDA)
###############################################################################

cat("\n", rep("=", 60), "\n")
cat("ANÁLISE ESTATÍSTICA DA TABELA CONSOLIDADA\n")
cat(rep("=", 60), "\n")

# Função segura para calcular estatísticas
calculate_row_stats <- function(row_data) {
  # Remover a primeira coluna (nome da métrica) e converter para numérico
  numeric_data <- as.numeric(row_data[-1])
  numeric_data <- numeric_data[!is.na(numeric_data)]
  
  if (length(numeric_data) == 0) {
    return(list(media = 0, desvio = 0, min_val = 0, max_val = 0, melhor_filtro = "N/A"))
  }
  
  # Encontrar o melhor filtro (maior valor)
  max_idx <- which.max(numeric_data)
  melhor_filtro <- names(row_data[-1])[max_idx]
  
  return(list(
    media = round(mean(numeric_data, na.rm = TRUE), 4),
    desvio = round(sd(numeric_data, na.rm = TRUE), 4),
    min_val = round(min(numeric_data, na.rm = TRUE), 4),
    max_val = round(max(numeric_data, na.rm = TRUE), 4),
    melhor_filtro = melhor_filtro
  ))
}

# Calcular estatísticas para cada métrica
cat("ESTATÍSTICAS POR MÉTRICA (Dataset GECCO):\n\n")

metrics_summary <- data.frame(
  Metrica = tabela_consolidada_gecco$Metrica,
  stringsAsFactors = FALSE
)

# Calcular estatísticas linha por linha
for (i in 1:nrow(tabela_consolidada_gecco)) {
  stats <- calculate_row_stats(tabela_consolidada_gecco[i, ])
  
  if (i == 1) {
    metrics_summary$Media <- stats$media
    metrics_summary$Desvio_Padrao <- stats$desvio
    metrics_summary$Min <- stats$min_val
    metrics_summary$Max <- stats$max_val
  } else {
    metrics_summary$Media[i] <- stats$media
    metrics_summary$Desvio_Padrao[i] <- stats$desvio
    metrics_summary$Min[i] <- stats$min_val
    metrics_summary$Max[i] <- stats$max_val
  }
}

print(metrics_summary)

# Melhores filtros por métrica
cat("\nMELHORES FILTROS POR MÉTRICA:\n")

best_filters <- data.frame(
  Metrica = tabela_consolidada_gecco$Metrica,
  stringsAsFactors = FALSE
)

for (i in 1:nrow(tabela_consolidada_gecco)) {
  stats <- calculate_row_stats(tabela_consolidada_gecco[i, ])
  
  if (i == 1) {
    best_filters$Melhor_Filtro <- stats$melhor_filtro
    best_filters$Valor <- stats$max_val
  } else {
    best_filters$Melhor_Filtro[i] <- stats$melhor_filtro
    best_filters$Valor[i] <- stats$max_val
  }
}

print(best_filters)

