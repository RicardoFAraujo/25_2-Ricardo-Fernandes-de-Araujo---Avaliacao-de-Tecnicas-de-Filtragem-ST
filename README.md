# Detecção de Anomalias em Séries Temporais (TCC)

**Título do TCC:** Avaliação de Técnicas de Filtragem para Detecção de Anomalias em Séries Temporais  
**Aluno:** Ricardo Fernandes de Araujo  
**Orientação:** Eduardo Soares Ogasawara, D.Sc. (Orientador) e Laís Baroni Villet, D.Sc. (Coorientadora)  
**Semestre de Defesa:** 2025.2

[PDF do TCC]


# TL;DR

Executar o pipeline do experimento (aplicar **filtros de pré-processamento** + rodar **detecção baseada em ARIMA** + **avaliar a detecção** vs. *ground truth* e exportar métricas em CSV) nos datasets Yahoo e GECCO:

```bash
Rscript run_yahoo_a1_a4_filters_arima_export_csv.R
Rscript run_gecco_filters_arima_export_csv.R
```


# Descrição Geral

Este repositório contém os scripts R usados no experimento do TCC para avaliar o impacto de técnicas de filtragem no pré-processamento de séries temporais para detecção supervisionada de anomalias com um detector baseado em ARIMA.

Conforme descrito no TCC, a arquitetura do experimento integra: **(i)** os datasets do repositório **UniTED** (via pacote R `united`), **(ii)** os filtros do framework **TSPredIT** (via pacote R `tspredit`) e **(iii)** a detecção baseada em ARIMA do framework **Harbinger** (via pacote R `harbinger`), usando uma versão **modificada** do detector `arima_filter` no arquivo `arima_filter.R` para receber o filtro como parâmetro.

Em alto nível, o fluxo é:

1) carregar os datasets (GECCO Industrial Challenge 2018 e Yahoo! Webscope S5 A1--A4);
2) aplicar filtros de pré-processamento;
3) ajustar ARIMA e executar detecção por eventos;
4) comparar com o ground truth e computar métricas (F1 como métrica central; além de precisão/recall/acurácia e erros de ajuste como MAE/RMSE);
5) exportar resultados em CSV.

Os resultados são sempre comparados com um baseline definido como **ARIMA sem filtragem**, e a melhoria é frequentemente expressa por ΔF1 = F1_filtro − F1_baseline.


# Frameworks e Integração

Este repositório integra quatro componentes principais:

* **UniTED / `united`**: fornece os datasets padronizados e rotulados usados no experimento (Yahoo! Webscope S5 e GECCO 2018).
* **TSPredIT / `tspredit`**: fornece os filtros de pré-processamento (ex.: Wavelet, EMD, Spline, HP, Kalman, LOWESS, Winsor, etc.).
* **Harbinger / `harbinger`**: fornece a infraestrutura do detector (distância L1, limiar gaussiano e interface de detecção), que é estendida no `arima_filter.R`.
* **`forecast`**: fornece o ajuste do modelo ARIMA (principalmente `auto.arima()` e `Arima()`), utilizado dentro do `arima_filter.R`.


# Funcionalidades

* Geração dos resultados do experimento em CSV (R)
* Execução no Yahoo (A1..A4) com múltiplos filtros e exportação por benchmark e combinado
* Execução no GECCO (recorte de 24h, variáveis explícitas) com múltiplos filtros e exportação completa por variável


# Filtros Avaliados

Conforme descrito no TCC, o experimento compara ARIMA sem filtragem (baseline) contra ARIMA acoplado a filtros de pré-processamento. Entre os filtros considerados estão métodos de média (MA/EMA), suavização exponencial (SES/QES), decomposição/multiescala (EMD/REMD/FFT/Wavelet), filtros adaptativos (Kalman/LOWESS/Recursive_Filter) e métodos robustos/estatísticos (HP/Spline/Winsor/Seasonal_Adj).

A lista exata (e possíveis variações/parametrizações) está explicitada nos próprios scripts `.R`.





# Como o arima_filter.R funciona (modificado)

O arquivo `arima_filter.R` implementa um detector `arima_filter(filtro)` que **estende o Harbinger** para aceitar um objeto de filtro (do TSPredIT) como parâmetro.

Em termos de pipeline:

* **Construção**: `arima_filter(filtro)` armazena `obj$filtro <- filtro` (ou `NULL` para o baseline sem filtragem).
* **Treino (`fit`)**:
    * remove `NA` da série;
    * se `obj$filtro` não é `NULL`, aplica `transform(obj$filtro, serie)`;
    * ajusta ARIMA com `forecast::auto.arima()` e guarda os parâmetros $(p,d,q)$ e a opção de drift para reutilizar na detecção.
* **Detecção (`detect`)**:
    * aplica o filtro apenas para **gerar a série de modelagem** (com salvaguardas: se o filtro falhar, inserir `NA` ou mudar o comprimento, o código reverte para a série original);
    * ajusta o ARIMA com os parâmetros estimados no `fit` (via `forecast::Arima()`, com fallback para `auto.arima()` se necessário);
    * calcula valores ajustados (`fitted`) e obtém resíduos como **série original − predição do modelo ajustado na série filtrada**;
    * aplica distância L1 e o detector de outliers do Harbinger para marcar anomalias.

Observação: o detector também ignora uma janela inicial de tamanho `sw_size` (relacionada a $(p,d,q)$) para evitar marcações espúrias no começo da série.

Esse desenho permite comparar, de forma controlada, **(i)** o baseline (ARIMA sem filtro) contra **(ii)** ARIMA acoplado a filtros, mantendo a etapa de detecção compatível com o Harbinger.


# Dependências

## R

Requer R e os pacotes:

* `united`
* `daltoolbox`
* `harbinger`
* `tspredit` (filtros do TSPredIT)
* `forecast` (modelo ARIMA)
* `dplyr` (para o script do Yahoo)
* `ggplot2`, `reshape2`, `corrplot` (usados no script do GECCO)

Observação: `arima_filter.R` é carregado via `source()` pelos scripts.

Instalação rápida (se aplicável no seu ambiente):

```r
install.packages(c(
    "united",
    "daltoolbox",
    "harbinger",
    "tspredit",
    "forecast",
    "dplyr",
    "ggplot2",
    "reshape2",
    "corrplot"
))
```


# Como Executar

## 1) Gerar CSVs do Yahoo (A1..A4)

Script principal: `run_yahoo_a1_a4_filters_arima_export_csv.R`.

```bash
Rscript run_yahoo_a1_a4_filters_arima_export_csv.R
```

Saídas típicas:

* `metrics_A1_todas_series.csv` (quando o modo auto-run de A1 está habilitado)
* `metrics_benchmarks/metrics_A1_A4_combinado.csv` e sumarizações em `metrics_benchmarks/` (se o auto-run multi-benchmark for habilitado)

## 2) Gerar CSVs do GECCO

Script principal: `run_gecco_filters_arima_export_csv.R`.

```bash
Rscript run_gecco_filters_arima_export_csv.R
```

Saídas:

* `resultados_detalhados_por_variavel.csv`
* `tabela_consolidada_gecco.csv`

Variáveis analisadas (explícitas): `tp`, `cl`, `ph`, `redox`, `leit`, `trueb`, `cl_2`, `fm`, `fm_2`.

## Observações

* Execute os comandos a partir da raiz do repositório (a mesma pasta onde estão os arquivos `.R`).
* Os scripts escrevem os CSVs no diretório de execução atual.
