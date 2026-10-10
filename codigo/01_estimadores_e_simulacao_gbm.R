# Estimadores e simulação GBM

library(tidyquant)
library(dplyr)
library(ggplot2)

#-------------------------------------------------------------------------------

calculo_vol <- function(ticker, from_date, to_date) { 
  # Recolha de dados
  dados <- tq_get(ticker, from = from_date, to = to_date) %>%
    select(date, open, high, low, close) %>%
    arrange(date) %>%
    mutate(
      u = log(high / open),
      d = log(low / open),
      c = log(close / open),
      li = log(high / low)) 
  
  dados <- dados %>%
    mutate(log_ret = log(close / lag(close)))
  
  dados <- na.omit(dados)  
  
  n_dias <- nrow(dados)
  
  # Método histórico
  sigma_historica <- sd(dados$log_ret) * sqrt(252)
  
  # Método - Parkinson
  const_park <- 1 / (4 * log(2)) 
  D_l <- (const_park / nrow(dados)) * sum(dados$li^2)
  sigma_parkinson <- sqrt(D_l) * sqrt(252)
  
  # Garman-Klass
  termos_gk <- 0.5 * (dados$u - dados$d)^2 - (2 * log(2) - 1) * dados$c^2
  var_media_gk <- mean(termos_gk)
  sigma_gk <- sqrt(var_media_gk) * sqrt(252)
  
  # Rogers-Satchell
  termos_rs <- dados$u * (dados$u - dados$c) + dados$d * (dados$d - dados$c)
  var_media_rs <- mean(termos_rs)
  sigma_rs <- sqrt(var_media_rs) * sqrt(252)
  
  print(
    ggplot(dados, aes(x = date, y = close)) +
      geom_line(color = "blue", size = 0.8) +
      labs(title = paste("Gráfico -", ticker),
           subtitle = paste("De", from_date, "a", to_date),
           x = "Data", y = "Preço") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
  )
  return(data.frame(
    n_observacoes = n_dias,
    ticker = ticker,
    vol_parkinson = sigma_parkinson,
    vol_gk = sigma_gk,
    vol_rs = sigma_rs,
    vol_historica = sigma_historica,
    ultimo_preco = tail(dados$close, 1) 
  ))
}

#-------------------------------------------------------

tsla_vol <- calculo_vol("TSLA", from_date = "2025-01-01", to_date = "2026-01-01"); tsla_vol 

meta_vol <- calculo_vol("META", from_date = "2025-01-01", to_date = "2026-01-01");meta_vol 
nvda_vol <- calculo_vol("NVDA", from_date = "2025-01-01", to_date = "2026-01-01");nvda_vol 
AMD_vol <- calculo_vol("AMD", from_date = "2025-01-01", to_date = "2026-01-01");AMD_vol 

Cola_vol <- calculo_vol("KO",  from_date = "2025-01-01", to_date = "2026-01-01"); Cola_vol 
McDonalds_vol <- calculo_vol("MCD",  from_date = "2025-01-01", to_date = "2026-01-01"); McDonalds_vol 

amzn_vol <- calculo_vol("AMZN",  from_date = "2025-01-01", to_date = "2026-01-01"); amzn_vol 
nike_vol <- calculo_vol("NKE",  from_date = "2025-01-01", to_date = "2026-01-01"); nike_vol 

BankAmerica_vol <- calculo_vol("BAC",  from_date = "2025-01-01", to_date = "2026-01-01"); BankAmerica_vol 

JNJ_vol <- calculo_vol("JNJ",  from_date = "2025-01-01", to_date = "2026-01-01"); JNJ_vol 

#-------------------------------------------------------

# Calculo do preço de opções de compra

# Volatilidades implícitas de mercado observadas a 02/01/2026.
# Fonte: AlphaQuery (https://www.alphaquery.com/)

#Fórmula de Black-Scholes - Call
Call <- function(S,K,TT,r,sigma){
  d1 <- (log(S/K)+(r+sigma^2/2)*TT)/(sigma*sqrt(TT))
  d2 <- d1-sigma*sqrt(TT)
  price <- S*pnorm(d1)-exp(-r*TT)*K*pnorm(d2)
  return(price)
}

precos_calls <- function(S, sigma, r = 0.01, T = 1) {
  K <- c(
    ITM = 0.97 * S,
    ATM = S,
    OTM = 1.03 * S
  )
  
  Call(S, K, T, r, sigma)
}

S <- tail(tsla_vol$ultimo_preco, 1);S 

# sigmas a usar
vol_parkinson_tsla <- tsla_vol$vol_parkinson;vol_parkinson_tsla 
vol_gk_tsla <- tsla_vol$vol_gk; vol_gk_tsla
vol_rs_tsla <- tsla_vol$vol_rs; vol_rs_tsla
vol_historica_tsla <- tsla_vol$vol_historica;vol_historica_tsla 

# Preço de Calls para TSLA - usando sigma_parkinson 
precos_calls(S, vol_parkinson_tsla)

#---------

# Preço de Calls para TSLA - usando sigma_historica 
precos_calls(S, vol_historica_tsla)

#----------

# Preço de Calls para TSLA - usando sigma_gk 
precos_calls(S, vol_gk_tsla)

#---------

# Preço de Calls para TSLA - usando sigma_rs 
precos_calls(S, vol_rs_tsla)

#-----------

# Preço de Calls para TSLA - usando sigma_implicita - 0.4796
vol_imp_tsla <- 0.4796;vol_imp_tsla
precos_calls(S, vol_imp_tsla)

#------------------

# Preço de Calls para JNJ: 

JNJ_vol

# sigmas a usar
vol_parkinson_JNJ <- JNJ_vol$vol_parkinson;vol_parkinson_JNJ 
vol_historica_JNJ <- JNJ_vol$vol_historica;vol_historica_JNJ 
vol_gk_JNJ <- JNJ_vol$vol_gk; vol_gk_JNJ
vol_rs_JNJ <- JNJ_vol$vol_rs; vol_rs_JNJ

S2 <- tail(JNJ_vol$ultimo_preco, 1);S2 

# Preço de Calls para JNJ - usando os vários sigmas
precos_calls(S2,vol_historica_JNJ)
precos_calls(S2,vol_parkinson_JNJ)
precos_calls(S2,vol_gk_JNJ)
precos_calls(S2,vol_rs_JNJ)

# Preço de Calls para JNJ - usando sigma_implicita - 0.2133
vol_imp_jnj <- 0.2133;vol_imp_jnj
precos_calls(S2,vol_imp_jnj)

#------------

meta_vol

# META
vol_historica_meta <- meta_vol$vol_historica;vol_historica_meta 
vol_parkinson_meta <- meta_vol$vol_parkinson;vol_parkinson_meta
vol_gk_meta <- meta_vol$vol_gk; vol_gk_meta
vol_rs_meta <- meta_vol$vol_rs; vol_rs_meta

S3 <- tail(meta_vol$ultimo_preco, 1);S3 

# Preço de Calls para META 
precos_calls(S3,vol_historica_meta)
precos_calls(S3,vol_parkinson_meta)
precos_calls(S3,vol_gk_meta)
precos_calls(S3,vol_rs_meta)

# Preço de Calls para META - usando sigma_implicita - 0.3759
vol_imp_meta <- 0.3759;vol_imp_meta
precos_calls(S3,vol_imp_meta)

#------------

# NVDA 

nvda_vol

vol_historica_nvda <- nvda_vol$vol_historica;vol_historica_nvda 
vol_parkinson_nvda <- nvda_vol$vol_parkinson;vol_parkinson_nvda 
vol_gk_nvda <- nvda_vol$vol_gk; vol_gk_nvda
vol_rs_nvda <- nvda_vol$vol_rs; vol_rs_nvda

S4 <- tail(nvda_vol$ultimo_preco, 1);S4 

# Preço de Calls para NVDA  
precos_calls(S4,vol_historica_nvda)
precos_calls(S4,vol_parkinson_nvda)
precos_calls(S4,vol_gk_nvda)
precos_calls(S4,vol_rs_nvda)

# Preço de Calls para NVDA - usando sigma_implicita - 0.3647
vol_imp_nvda <- 0.3647;vol_imp_nvda
precos_calls(S4,vol_imp_nvda)

#------------

# AMD

AMD_vol

vol_historica_AMD <- AMD_vol$vol_historica;vol_historica_AMD 
vol_parkinson_AMD <- AMD_vol$vol_parkinson;vol_parkinson_AMD 
vol_gk_AMD <- AMD_vol$vol_gk; vol_gk_AMD
vol_rs_AMD <- AMD_vol$vol_rs; vol_rs_AMD

S5 <- tail(AMD_vol$ultimo_preco, 1);S5 

# Preço de Calls para AMD 
precos_calls(S5,vol_historica_AMD)
precos_calls(S5,vol_parkinson_AMD)
precos_calls(S5,vol_gk_AMD)
precos_calls(S5,vol_rs_AMD)

# Preço de Calls para AMD - usando sigma_implicita - 0.4963
vol_imp_AMD <- 0.4963;vol_imp_AMD
precos_calls(S5,vol_imp_AMD)

#-------------

# KO

Cola_vol

vol_historica_Cola <- Cola_vol$vol_historica;vol_historica_Cola 
vol_parkinson_Cola <- Cola_vol$vol_parkinson;vol_parkinson_Cola 
vol_gk_Cola <- Cola_vol$vol_gk; vol_gk_Cola
vol_rs_Cola <- Cola_vol$vol_rs; vol_rs_Cola

S6 <- tail(Cola_vol$ultimo_preco, 1);S6 

# Preço de Calls para KO  
precos_calls(S6,vol_historica_Cola)
precos_calls(S6,vol_parkinson_Cola)
precos_calls(S6,vol_gk_Cola)
precos_calls(S6,vol_rs_Cola)

# Preço de Calls para KO - usando sigma_implicita - 0.1376
vol_imp_Cola <- 0.1376;vol_imp_Cola
precos_calls(S6,vol_imp_Cola)

#-------------

# MCD

McDonalds_vol

vol_historica_McD <- McDonalds_vol$vol_historica;vol_historica_McD 
vol_parkinson_McD <- McDonalds_vol$vol_parkinson;vol_parkinson_McD 
vol_gk_McD <- McDonalds_vol$vol_gk; vol_gk_McD
vol_rs_McD <- McDonalds_vol$vol_rs; vol_rs_McD

S7 <- tail(McDonalds_vol$ultimo_preco, 1);S7 

# Preço de Calls para MCD  
precos_calls(S7,vol_historica_McD)
precos_calls(S7,vol_parkinson_McD)
precos_calls(S7,vol_gk_McD)
precos_calls(S7,vol_rs_McD)

# Preço de Calls para MCD - usando sigma_implicita - 0.1618
vol_imp_McD <- 0.1618;vol_imp_McD
precos_calls(S7,vol_imp_McD)

#--------------------

# AMZN

amzn_vol

vol_historica_amzn <- amzn_vol$vol_historica;vol_historica_amzn 
vol_parkinson_amzn <- amzn_vol$vol_parkinson;vol_parkinson_amzn 
vol_gk_amzn <- amzn_vol$vol_gk; vol_gk_amzn
vol_rs_amzn <- amzn_vol$vol_rs; vol_rs_amzn

S8 <- tail(amzn_vol$ultimo_preco, 1);S8 

# Preço de Calls para AMZN 
precos_calls(S8,vol_historica_amzn)
precos_calls(S8,vol_parkinson_amzn)
precos_calls(S8,vol_gk_amzn)
precos_calls(S8,vol_rs_amzn)

# Preço de Calls para AMZN - usando sigma_implicita - 0.2953
vol_imp_amzn <- 0.2953;vol_imp_amzn
precos_calls(S8,vol_imp_amzn)


#---------------

# NIKE

nike_vol

vol_historica_nike <- nike_vol$vol_historica;vol_historica_nike 
vol_parkinson_nike <- nike_vol$vol_parkinson;vol_parkinson_nike 
vol_gk_nike <- nike_vol$vol_gk; vol_gk_nike
vol_rs_nike <- nike_vol$vol_rs; vol_rs_nike

S9 <- tail(nike_vol$ultimo_preco, 1);S9 

# Preço de Calls para NIKE  
precos_calls(S9,vol_historica_nike)
precos_calls(S9,vol_parkinson_nike)
precos_calls(S9,vol_gk_nike)
precos_calls(S9,vol_rs_nike)

# Preço de Calls para NIKE - usando sigma_implicita - 0.2895
vol_imp_nike <- 0.2895;vol_imp_nike
precos_calls(S9,vol_imp_nike)

#-------------

# BAC 

BankAmerica_vol

vol_historica_bac <- BankAmerica_vol$vol_historica;vol_historica_bac 
vol_parkinson_bac <- BankAmerica_vol$vol_parkinson;vol_parkinson_bac 
vol_gk_bac <- BankAmerica_vol$vol_gk; vol_gk_bac
vol_rs_bac <- BankAmerica_vol$vol_rs; vol_rs_bac

S10 <- tail(BankAmerica_vol$ultimo_preco, 1);S10 

# Preço de Calls para BAC - usando sigma_historica 
precos_calls(S10,vol_historica_bac)
precos_calls(S10,vol_parkinson_bac)
precos_calls(S10,vol_gk_bac)
precos_calls(S10,vol_rs_bac)

# Preço de Calls para BAC - usando sigma_implicita - 0.2501
vol_imp_bac <- 0.2501;vol_imp_bac
precos_calls(S10,vol_imp_bac)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# Dados simulados

# Simular trajetórias de preços - GBM

sim_trajetoria <- function(S0, mu, sigma, dias, passos_dia, delta_t) {
  n <- dias * passos_dia 
  
  retornos <- rnorm(n, mean = (mu - 0.5 * sigma^2) * delta_t, sd = sigma * sqrt(delta_t))
  log_precos <- c(log(S0), cumsum(retornos) + log(S0)) 
  precos <- exp(log_precos)
  
  precos_por_dia <- matrix(precos[-1], nrow = passos_dia)  
  fechos <- precos_por_dia[passos_dia, ] 
  aberturas <- c(S0, fechos[1:(dias-1)]) 
  maximos <- apply(precos_por_dia, 2, max) 
  minimos <- apply(precos_por_dia, 2, min)
  
  return(data.frame( 
    dia = 1:dias,
    open = aberturas,
    high = maximos,
    low = minimos,
    close = fechos
  ))
}

#----------

# Estimar volatilidades pelos diferentes métodos

estimar_volatilidades <- function(dados) { 
  
  # Parkinson
  li <- log(dados$high / dados$low)
  const_park <- 1 / (4 * log(2)) 
  D_l <- (const_park / nrow(dados)) * sum(li^2)
  sigma_parkinson <- sqrt(D_l) * sqrt(252) 
  
  # Histórico
  log_ret <- diff(log(dados$close)) 
  sigma_historica <- sd(log_ret) * sqrt(252) 
  
  # variáveis para GK e RS
  u <- log(dados$high / dados$open)
  d <- log(dados$low / dados$open)
  c <- log(dados$close / dados$open)
  
  # Garman-Klass
  termos_gk <- 0.5 * (u - d)^2 - (2 * log(2) - 1) * c^2
  var_media_gk <- mean(termos_gk)
  sigma_gk <- sqrt(var_media_gk) * sqrt(252)
  
  # Rogers-Satchell 
  termos_rs <- u * (u - c) + d * (d - c)
  var_media_rs <- mean(termos_rs)
  sigma_rs <- sqrt(var_media_rs) * sqrt(252)
  
  return(c(vol_parkinson = sigma_parkinson, 
           vol_historica = sigma_historica,
           vol_gk = sigma_gk,
           vol_rs = sigma_rs))
}

set.seed(110)
nsimul <- 1000 

S0 <- 100 
mu <- 0.5 
dias <- 252

# Cenários considerados na dissertação:
# sigma_real = 0.2, 0.4, 0.6, 0.8, 1.0
# passos_dia = 100 e 1000
#
# Alterar os dois parâmetros abaixo para reproduzir cada cenário.
sigma_real <- 1 
passos_dia <- 100 

delta_t <- 1/(dias*passos_dia)

vols_parkinson <- numeric(nsimul)
vols_historica <- numeric(nsimul)
vols_gk <- numeric(nsimul)       
vols_rs <- numeric(nsimul)

for (i in 1:nsimul) {
  traj <- sim_trajetoria(S0, mu, sigma_real, dias, passos_dia, delta_t) 
  vols <- estimar_volatilidades(traj) 
  vols_parkinson[i] <- vols[1]
  vols_historica[i] <- vols[2]
  vols_gk[i] <- vols[3]          
  vols_rs[i] <- vols[4]          
}
mean(vols_historica)
mean(vols_parkinson)  
mean(vols_gk)
mean(vols_rs)

sigma_real

# enviesamento
enviesamento_historica <- mean(vols_historica) - sigma_real;enviesamento_historica
enviesamento_parkinson <- mean(vols_parkinson) - sigma_real;enviesamento_parkinson
enviesamento_gk <- mean(vols_gk) - sigma_real;enviesamento_gk
enviesamento_rs <- mean(vols_rs) - sigma_real;enviesamento_rs

# desvios padrão
sd(vols_historica) 
sd(vols_parkinson) 
sd(vols_gk)
sd(vols_rs)

# variância
sd(vols_historica)^2
sd(vols_parkinson)^2 
sd(vols_gk)^2
sd(vols_rs)^2

#-------

# EQM - método histórico:

enviesamento_historica^2
EQM_historica <- enviesamento_historica^2 + var(vols_historica);EQM_historica 

#------

# EQM -  Parkinson:

enviesamento_parkinson^2
EQM_parkinson <- enviesamento_parkinson^2 + var(vols_parkinson);EQM_parkinson 

#-------

# EQM - Garman-Klass

enviesamento_gk^2
EQM_gk <- enviesamento_gk^2 + var(vols_gk);EQM_gk

#-------

# EQM - Rogers-Satchell

enviesamento_rs^2
EQM_rs <- enviesamento_rs^2 + var(vols_rs);EQM_rs

#--------

resumo_final <- data.frame(
  Metodo = c("Histórico", "Parkinson", "Garman-Klass", "Rogers-Satchell"),
  Media = c(mean(vols_historica), mean(vols_parkinson), mean(vols_gk), mean(vols_rs)),
  Enviesamento = c(enviesamento_historica, enviesamento_parkinson, enviesamento_gk, enviesamento_rs),
  Desvio_Padrao = c(sd(vols_historica), sd(vols_parkinson), sd(vols_gk), sd(vols_rs)),
  Variancia = c(var(vols_historica), var(vols_parkinson), var(vols_gk), var(vols_rs)),
  EQM = c(EQM_historica, EQM_parkinson, EQM_gk, EQM_rs)
)

print(resumo_final)

