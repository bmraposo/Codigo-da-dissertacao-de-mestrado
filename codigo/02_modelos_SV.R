# Modelos de Volatilidade Estocástica - estimação, simulação de trajetórias e 
# medição de risco 

library(tidyquant)
library(dplyr)
library(ggplot2)

calculo_var_diaria <- function(ticker, from_date, to_date) {
  
  # Recolha de dados
  dados <- tq_get(ticker, from = from_date, to = to_date) %>%
    select(date, open, high, low, close) %>%
    arrange(date) %>%
    mutate(
      u = log(high / open),
      d = log(low / open),
      c = log(close / open),
      li = log(high / low)
    ) %>%
    na.omit()
  
  # Parkinson
  dados <- dados %>%
    mutate(
      var_park = (1 / (4 * log(2))) * li^2)
  
  # Garman-Klass
  dados <- dados %>%
    mutate(var_gk = 0.5 * (u - d)^2 -(2*log(2)-1)*c^2)
  
  # Rogers-Satchell
  dados <- dados %>%
    mutate(var_rs = u * (u - c) + d * (d - c))
  
  
  print(
    ggplot(dados, aes(x = date, y = close)) +
      geom_line(color = "blue", size = 0.6) +
      labs(title = paste("Gráfico -", ticker),
           subtitle = paste("De", from_date, "a", to_date),
           x = "Data", y = "Preço") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
  )
  
  return(dados %>%
           select(date, open, high, low, close, var_park, var_gk, var_rs))
}

tsla_var_diaria <- calculo_var_diaria("TSLA", from_date = "2025-01-01", to_date = "2026-01-01")

# Estatisticas descritivas das series de variância diária - park, gk, rs
summary(tsla_var_diaria$var_park)
summary(tsla_var_diaria$var_gk)
summary(tsla_var_diaria$var_rs)
sd(tsla_var_diaria$var_park)
sd(tsla_var_diaria$var_gk)
sd(tsla_var_diaria$var_rs)

# Estatisticas descritivas das series de volatilidade diária - park, gk, rs
summary(sqrt(tsla_var_diaria$var_park))
summary(sqrt(tsla_var_diaria$var_gk))
summary(sqrt(tsla_var_diaria$var_rs))
sd(sqrt(tsla_var_diaria$var_park))
sd(sqrt(tsla_var_diaria$var_gk))
sd(sqrt(tsla_var_diaria$var_rs))

#---------------------------------
#---------------------------------

# Estimar os parâmetros dos modelos 

#-------------------------------------------------------------------------
# Hull-White
#-------------------------------------------------------------------------

estimacao_param_hw <- function(V, delta = 1) { 
  
  # calculo r_i
  r <- diff(log(V))
  
  n <- length(r)
  
  # média amostral dos retornos logarítmicos 
  r_bar <- mean(r)
  
  # estimador de sigma 
  soma_quadrados <- sum((r - r_bar)^2)
  sigma2_mle <- soma_quadrados / (n * delta)
  sigma_mle  <- sqrt(sigma2_mle)
  
  # estimador de alpha 
  alpha_mle <- (r_bar / delta) + 0.5 * sigma2_mle
  
  return(list(
    alpha = alpha_mle, 
    sigma   = sigma_mle  
  ))
}

estimacao_param_hw(tsla_var_diaria$var_park)
estimacao_param_hw(tsla_var_diaria$var_gk)
estimacao_param_hw(tsla_var_diaria$var_rs)

#-------------------------------------------------------------------------
# Heston 
#-------------------------------------------------------------------------

estimacao_param_heston <- function(V, delta = 1) { 
  
  # V_i e V_{i+1}
  V_i  <- V[-length(V)]
  V_ip <- V[-1]
  
  n <- length(V_i)
  
  # cálculos auxiliares
  sum_Vi        <- sum(V_i)
  sum_inv_Vi    <- sum(1 / V_i) 
  sum_Vip_invVi <- sum(V_ip / V_i) 
  
  # estimador de alpha
  num <- sum(V_ip) * sum_inv_Vi - n * sum_Vip_invVi
  den <- sum_Vi * sum_inv_Vi - n^2
  
  alpha_hat <- -(1 / delta) * log(num / den)
  
  # estimador de beta
  beta_hat <- (sum_Vip_invVi - n * exp(-alpha_hat * delta)) /
    ((1 - exp(-alpha_hat * delta)) * sum_inv_Vi)
  
  # estimador de sigma^2
  aux <- V_ip -
    V_i * exp(-alpha_hat * delta) -
    beta_hat * (1 - exp(-alpha_hat * delta))
  
  num_sigma <- 2 * alpha_hat * sum((aux^2) / V_i)
  
  den_sigma <- 2 * n * exp(-alpha_hat * delta) * (1 - exp(-alpha_hat * delta)) +
    (exp(-alpha_hat * delta) - 1)^2 * beta_hat * sum_inv_Vi
  
  sigma2_hat <- num_sigma / den_sigma
  sigma_hat  <- sqrt(sigma2_hat)
  
  return(list(
    alpha = alpha_hat,
    beta  = beta_hat,
    sigma = sigma_hat
  ))
}


estimacao_param_heston(tsla_var_diaria$var_park)
estimacao_param_heston(tsla_var_diaria$var_gk)
estimacao_param_heston(tsla_var_diaria$var_rs)

#-------------------------------------------------------------------------
# Ornstein-Uhlenbeck
#-------------------------------------------------------------------------

estimacao_param_OU <- function(V, delta = 1) {
  
  # V_i e V_{i+1}
  V_i  <- V[-length(V)] 
  V_ip <- V[-1] 
  
  n <- length(V_i)
  
  # somatórios
  Sx  <- sum(V_i)
  Sy  <- sum(V_ip)
  Sxx <- sum(V_i^2)
  Sxy <- sum(V_i * V_ip)
  
  # estimador de alpha
  alpha_hat <- -(1 / delta) * log((n * Sxy - Sy * Sx) / (n * Sxx - Sx^2))
  
  # estimador de beta
  beta_hat <- (Sy * Sxx - Sxy * Sx) /
    (n * (Sxx - Sxy) + Sx * (Sy - Sx))
  
  # estimador de sigma^2
  aux <- V_ip -
    V_i * exp(-alpha_hat * delta) -
    beta_hat * (1 - exp(-alpha_hat * delta))
  
  sigma2_hat <- ((2 * alpha_hat) /
                   (n * (1 - exp(-2 * alpha_hat * delta)))) * sum(aux^2)
  
  sigma_hat <- sqrt(sigma2_hat)
  
  return(list(
    alpha = alpha_hat,
    beta  = beta_hat,
    sigma = sigma_hat
  ))
}

estimacao_param_OU(tsla_var_diaria$var_park)
estimacao_param_OU(tsla_var_diaria$var_gk)
estimacao_param_OU(tsla_var_diaria$var_rs)

#-------------------------------------------------------------------------
# Ornstein-Uhlenbeck Geométrico
#-------------------------------------------------------------------------

estimacao_param_GOU <- function(V, delta = 1) {
  
  # transformação de variável
  X <- log(V)
  
  # X_i e X_{i+1} 
  X_i  <- X[-length(X)]  
  X_ip <- X[-1]  
  
  n <- length(X_i) 
  
  # somatórios
  Sx  <- sum(X_i)
  Sy  <- sum(X_ip)
  Sxx <- sum(X_i^2)
  Sxy <- sum(X_i * X_ip)
  
  # estimador de alpha 
  alpha_hat <- -(1 / delta) * log((n * Sxy - Sy * Sx) / (n * Sxx - Sx^2))
  
  # estimador de beta 
  beta_hat <- (Sy * Sxx - Sxy * Sx) / 
    (n * (Sxx - Sxy) + Sx * (Sy - Sx))
  
  # estimador de sigma^2 
  phi <- exp(-alpha_hat * delta)
  aux <- X_ip - X_i * phi - beta_hat * (1 - phi)
  
  sigma2_hat <- ((2 * alpha_hat) / 
                   (n * (1 - exp(-2 * alpha_hat * delta)))) * sum(aux^2)
  
  sigma_hat <- sqrt(sigma2_hat)
  
  return(list(
    alpha = alpha_hat, 
    beta  = beta_hat,  
    sigma = sigma_hat  
  ))
}

estimacao_param_GOU(tsla_var_diaria$var_park)
estimacao_param_GOU(tsla_var_diaria$var_gk)
estimacao_param_GOU(tsla_var_diaria$var_rs)

#---------------------------------------------------------
# Cálculo de \mu e \rho para a simulação de trajetórias de cada modelo
# Base: Série de variâncias diárias de Rogers-Satchell (var_rs)
#---------------------------------------------------------

log_retornos <- diff(log(tsla_var_diaria$close)) 

# calculo \mu 
mu_diario     <- mean(log_retornos, na.rm = TRUE) + 
  0.5 * mean(tsla_var_diaria$var_rs[-nrow(tsla_var_diaria)],
             na.rm = TRUE); mu_diario
mu_anualizado <- mu_diario * 252; mu_anualizado 

#--------

# Calculo rho 

choque_preco <- function(log_retornos, mu, V, delta = 1){
  Z_S <- (log_retornos - (mu - 0.5 * V) * delta) / sqrt(V * delta)
  return(Z_S)
}

V_preco <- tsla_var_diaria$var_rs[-nrow(tsla_var_diaria)] 

# choque do preço - igual pra todos os modelos 
Z_S <- choque_preco(log_retornos, mu_diario, V_preco,1)

#----------------------

# O choque da variância (Z_V) é diferente para cada modelo - *Apêndice*

#-----------------------------------------
# Modelo HW 
#-----------------------------------------

choque_var_hw <- function(V, alpha, sigma, delta = 1){
  delta_log_V <- diff(log(V))
  
  Z_V <- (delta_log_V - (alpha - 0.5 * sigma^2) * delta) / 
    (sigma * sqrt(delta))
  
  return(Z_V)
}

param_hw_rs <- estimacao_param_hw(tsla_var_diaria$var_rs);param_hw_rs

Z_V_hw <- choque_var_hw(tsla_var_diaria$var_rs, param_hw_rs$alpha,
                        param_hw_rs$sigma)
rho_hw <- cor(Z_S, Z_V_hw);rho_hw

#-----------------------------------------
# Modelo Heston 
#-----------------------------------------

choque_var_he <- function(V, alpha, beta, sigma, delta = 1){
  V_t  <- V[-length(V)] 
  
  V_pos <- pmax(V_t, 0)
  
  delta_V <- diff(V)
  
  Z_V <- (delta_V - alpha * (beta - V_pos) * delta) /
    (sigma * sqrt(V_pos * delta))
  return(Z_V)
}

param_he_rs <- estimacao_param_heston(tsla_var_diaria$var_rs);param_he_rs

Z_V_he <- choque_var_he(tsla_var_diaria$var_rs, param_he_rs$alpha, 
                        param_he_rs$beta, param_he_rs$sigma)

rho_he <- cor(Z_S, Z_V_he);rho_he  

#-----------------------------------------
# Modelo OU 
#-----------------------------------------

choque_var_ou <- function(V, alpha, beta, sigma, delta = 1){
  
  V_t <- V[-length(V)]      
  delta_V <- diff(V)        
  
  Z_V <- (delta_V - alpha * (beta - V_t) * delta) /
    (sigma * sqrt(delta))
  
  return(Z_V)
}

param_ou_rs <- estimacao_param_OU(tsla_var_diaria$var_rs);param_ou_rs

Z_V_ou <- choque_var_ou(tsla_var_diaria$var_rs, param_ou_rs$alpha, 
                        param_ou_rs$beta, param_ou_rs$sigma)

rho_ou <- cor(Z_S, Z_V_ou);rho_ou

#-----------------------------------------
# Modelo GOU 
#-----------------------------------------

choque_var_gou <- function(V, alpha, beta, sigma, delta = 1){
  
  X <- log(V)
  
  X_t <- X[-length(X)]      
  delta_X <- diff(X)        
  
  Z_X <- (delta_X - alpha * (beta - X_t) * delta) /
    (sigma * sqrt(delta))
  
  return(Z_X)
}

param_gou_rs <- estimacao_param_GOU(tsla_var_diaria$var_rs);param_gou_rs

# choque do logaritmo da variância
Z_X_gou <- choque_var_gou(tsla_var_diaria$var_rs, param_gou_rs$alpha, 
                          param_gou_rs$beta, param_gou_rs$sigma)

rho_gou <- cor(Z_S, Z_X_gou);rho_gou

# Resumindo
rho_hw 
rho_he 
rho_ou 
rho_gou 

#---------------------------------------------------------
# Valores iniciais (S0 e V0)
#---------------------------------------------------------

# O último preço e a última variância da amostra servem de ponto de partida
S0 <- tail(tsla_var_diaria$close, 1);S0
V0 <- tail(tsla_var_diaria$var_rs, 1);V0 

#----------------------------------------------------------
# Simulação de trajetórias
#----------------------------------------------------------
set.seed(110)
S0                  
N_sim   <- 10000     
N_steps <- 252 # no fim diminuir para 50 dias para testar a sensibilidade ao horizonte temporal em Hull-White     
delta_t   <- 1      
mu <- mu_diario  

#-------------------------------------------------------------------------
# Função auxiliar: geração de variáveis aleatórias correlacionadas
#-------------------------------------------------------------------------

# Gera matrizes de ruído para 10000 simulações x 252 dias
# Z1: Choque no Preço
# Z2: Choque na Variância (Correlacionado com Z1)

gerar_choques <- function(n_sim, n_steps, rho_hw, rho_he, rho_ou, rho_gou) {
  
  epsilon1 <- matrix(rnorm(n_sim * n_steps), nrow = n_sim, ncol = n_steps)
  epsilon2 <- matrix(rnorm(n_sim * n_steps), nrow = n_sim, ncol = n_steps)
  
  # Formulas de Z1, Z2_hw, Z2_he, Z2_ou, Z2_gou - transformação de Cholesky
  Z1 <- epsilon1
  Z2_hw  <- rho_hw  * epsilon1 + sqrt(1-rho_hw^2 ) * epsilon2
  Z2_he  <- rho_he  * epsilon1 + sqrt(1-rho_he^2 ) * epsilon2
  Z2_ou  <- rho_ou  * epsilon1 + sqrt(1-rho_ou^2 ) * epsilon2
  Z2_gou <- rho_gou * epsilon1 + sqrt(1-rho_gou^2) * epsilon2
  
  return(list(Z1 = Z1, Z2_hw = Z2_hw, Z2_he = Z2_he, Z2_ou = Z2_ou, 
              Z2_gou = Z2_gou))
}

choques <- gerar_choques(N_sim, N_steps, rho_hw, rho_he, rho_ou, rho_gou)
mat_Z1 <- choques$Z1
mat_Z2_hw  <- choques$Z2_hw
mat_Z2_he  <- choques$Z2_he
mat_Z2_ou  <- choques$Z2_ou
mat_Z2_gou <- choques$Z2_gou

#-------------------------------------------------------------------------
# Modelo Hull-White 
#-------------------------------------------------------------------------

param_hw_rs

alpha_hw <- param_hw_rs$alpha;alpha_hw # dividir por 252 para testar a sensibilidade aos parâmetros 
sigma_hw <- param_hw_rs$sigma;sigma_hw # dividir por sqrt(252) testar a sensibilidade aos parâmetros

# Matrizes para guardar resultados
S_hw <- matrix(NA, nrow=N_sim, ncol=N_steps+1); S_hw[,1] <- S0
V_hw <- matrix(NA, nrow=N_sim, ncol=N_steps+1); V_hw[,1] <- V0

for(t in 1:N_steps) {
  sqrt_V <- sqrt(V_hw[, t])
  
  # atualizar preço 
  S_hw[, t+1] <- S_hw[, t] * exp((mu - 0.5 * V_hw[, t]) * delta_t + sqrt_V * sqrt(delta_t) * mat_Z1[, t])
  
  # atualizar variância
  V_hw[, t+1] <- V_hw[, t] * exp((alpha_hw - 0.5 * sigma_hw^2) * delta_t + sigma_hw * sqrt(delta_t) * mat_Z2_hw[, t])
}

#-------------------------------------------------------------------------
# Modelo Heston 
#-------------------------------------------------------------------------

param_he_rs

alpha_he <- param_he_rs$alpha;alpha_he     
beta_he  <- param_he_rs$beta;beta_he      
sigma_he <- param_he_rs$sigma;sigma_he     

S_he <- matrix(NA, nrow=N_sim, ncol=N_steps+1); S_he[,1] <- S0
V_he <- matrix(NA, nrow=N_sim, ncol=N_steps+1); V_he[,1] <- V0

for(t in 1:N_steps) {
  
  # garantir que não há raiz de números negativos
  V_t <- pmax(V_he[, t], 0) 
  sqrt_V <- sqrt(V_t)
  
  # atualizar preço
  S_he[, t+1] <- S_he[, t] * exp((mu - 0.5 * V_t) * delta_t + sqrt_V * sqrt(delta_t) * mat_Z1[, t])
  
  # atualizar variância 
  V_next <- V_he[, t] + alpha_he * (beta_he - V_t) * delta_t + sigma_he * sqrt_V * sqrt(delta_t) * mat_Z2_he[, t]
  V_he[, t+1] <- pmax(V_next,0)
}

#-------------------------------------------------------------------------
# Modelo Ornstein-Uhlenbeck 
#-------------------------------------------------------------------------

param_ou_rs

alpha_ou <- param_ou_rs$alpha;alpha_ou 
beta_ou  <- param_ou_rs$beta;beta_ou     
sigma_ou <- param_ou_rs$sigma;sigma_ou   

S_ou <- matrix(NA, nrow=N_sim, ncol=N_steps+1); S_ou[,1] <- S0
V_ou <- matrix(NA, nrow=N_sim, ncol=N_steps+1); V_ou[,1] <- V0

for(t in 1:N_steps) {
  V_t <- V_ou[, t]
  sqrt_V <- sqrt(pmax(V_t, 0)) 
  
  # atualizar preço
  S_ou[, t+1] <- S_ou[, t] * exp((mu - 0.5 * pmax(V_t, 0)) * delta_t + sqrt_V * sqrt(delta_t) * mat_Z1[, t])
  
  # atualizar variância 
  V_ou[, t+1] <- V_t + alpha_ou * (beta_ou - V_t) * delta_t + sigma_ou * sqrt(delta_t) * mat_Z2_ou[, t]
}

#-------------------------------------------------------------------------
# Modelo Ornstein-Uhlenbeck geométrico 
#-------------------------------------------------------------------------

alpha_gou <- param_gou_rs$alpha;alpha_gou 
beta_gou  <- param_gou_rs$beta;beta_gou  
sigma_gou <- param_gou_rs$sigma;sigma_gou   

S_gou <- matrix(NA, nrow=N_sim, ncol=N_steps+1); S_gou[,1] <- S0
V_gou <- matrix(NA, nrow=N_sim, ncol=N_steps+1); V_gou[,1] <- V0

for(t in 1:N_steps) {
  V_t <- V_gou[, t]
  sqrt_V <- sqrt(V_t) 
  
  # Atualizar Preço 
  S_gou[, t+1] <- S_gou[, t] * exp((mu - 0.5 * V_t) * delta_t + sqrt_V * sqrt(delta_t) * mat_Z1[, t])
  
  # Atualizar Variância
  log_V_t <- log(V_t)
  log_V_next <- log_V_t + alpha_gou * (beta_gou - log_V_t) * delta_t + sigma_gou * sqrt(delta_t) * mat_Z2_gou[, t]
  
  V_gou[, t+1] <- exp(log_V_next)
}


#-------------------------------------------------------------------------
# Gráficos das trajetórias
#-------------------------------------------------------------------------

library(ggplot2)
library(reshape2)

plot_trajectories <- function(matriz, model_name, param, y_label = "Preço Simulado (USD)", ylim_range = NULL) {
  
  # Seleciona as primeiras 10 trajetórias
  df_plot <- as.data.frame(t(matriz[1:10, , drop = FALSE])) 
  df_plot$Time <- 0:(nrow(df_plot) - 1)
  df_long <- melt(df_plot, id="Time")
  
  p <- print(ggplot(df_long, aes(x=Time, y=value, group=variable)) +
               geom_line(alpha=0.5, color="blue") +
               geom_hline(yintercept = param, linetype="dashed", color="red") +
               labs(title = paste("Simulação -", model_name), 
                    y = y_label, 
                    x = "Dias (t)") +
               theme_minimal())
  if (!is.null(ylim_range)) { 
    p <- p + coord_cartesian(ylim = ylim_range)
  }
  
  print(p)
}

min_ou <- min(V_ou[1:10, ]);min_ou
max_heston <- max(V_he[1:10, ]);max_heston
lims <- c(min_ou * 0.9, max_heston * 1.1);lims

plot_trajectories(S_hw, "Hull-White", S0, "Preço (USD)")
plot_trajectories(S_he, "Heston", S0,"Preço (USD)")
plot_trajectories(S_ou, "Ornstein-Uhlenbeck", S0, "Preço (USD)")
plot_trajectories(S_gou, "Geométrico Ornstein-Uhlenbeck (GOU)", S0, "Preço (USD)")

plot_trajectories(V_hw, "Hull-White", V0, "Variância Diária")
plot_trajectories(V_he, "Heston", V0, "Variância Diária", lims) 
plot_trajectories(V_ou, "Ornstein-Uhlenbeck", V0, "Variância Diária",lims)
plot_trajectories(V_gou, "Geométrico Ornstein-Uhlenbeck (GOU)", V0, "Variância Diária")

# Estatísticas descritivas das simulações
S_hw_final <- S_hw[, ncol(S_hw)]
summary(S_hw_final)
V_hw_final <- V_hw[, ncol(V_hw)]
summary(V_hw_final)

S_he_final <- S_he[, ncol(S_he)]
summary(S_he_final)
V_he_final <- V_he[, ncol(V_he)]
summary(V_he_final)

S_ou_final <- S_ou[, ncol(S_ou)]
summary(S_ou_final)
V_ou_final <- V_ou[, ncol(V_ou)]
summary(V_ou_final)

S_gou_final <- S_gou[, ncol(S_gou)]
summary(S_gou_final)
V_gou_final <- V_gou[, ncol(V_gou)]
summary(V_gou_final)

S_esperado <- S0*exp(mu*252);S_esperado # próximo da media dos preços simulados

# Contar quantas trajetórias colapsaram para menos de 1 dólar
trajetorias_mortas <- sum(S_hw_final < 1, na.rm = TRUE)
percentagem_morte <- (trajetorias_mortas / length(S_hw_final)) * 100

cat("Número de trajetórias que colapsaram:", trajetorias_mortas, "\n")
cat("Percentagem de colapso:", round(percentagem_morte, 2), "%\n")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# VaR e Expected Shortfall

calcular_risco <- function(S_matrix, nome_modelo) {
  
  precos_finais <- S_matrix[, ncol(S_matrix)]
  retornos <- (precos_finais / S0) - 1
  
  var_95 <- quantile(retornos, 0.05) 
  es_95 <- mean(retornos[retornos <= var_95])
  
  var_99 <- quantile(retornos, 0.01)
  es_99 <- mean(retornos[retornos <= var_99])
  
  return(data.frame(
    Modelo = nome_modelo,
    VaR_95 = round(var_95 * -100, 2),
    ES_95 = round(es_95 * -100, 2),
    VaR_99 = round(var_99 * -100, 2),
    ES_99 = round(es_99 * -100, 2)
  ))
}

tabela_risco <- rbind(
  calcular_risco(S_he, "Heston"),
  calcular_risco(S_ou, "Ornstein-Uhlenbeck"),
  #calcular_risco(S_hw, "Hull-White"),
  calcular_risco(S_gou, "Ornstein-Uhlenbeck Geom.")
)
print(tabela_risco)

#----------- 

# VaR e ES teórico 

calcular_risco_teorico <- function(mu, beta, gamma, T_dias = 252) {
  
  sigma_bs <- sqrt(beta) 
  media_teorica <- (mu - 0.5 * sigma_bs^2) * T_dias
  desvio_teorico <- sigma_bs * sqrt(T_dias)
  
  quantil_log <- qnorm(gamma, mean = media_teorica, sd = desvio_teorico)
  var_teorico <- -(exp(quantil_log) - 1) * 100
  
  z_gamma <- qnorm(gamma)
  fator_exp <- exp(media_teorica + 0.5 * desvio_teorico^2)
  es_teorico <- (1 - (fator_exp / gamma) * pnorm(z_gamma - desvio_teorico)) * 100
  
  resultados <- data.frame(
    Nivel_confianca = paste0((1- gamma) * 100, "%"),
    VaR_Teorico = round(var_teorico, 2),
    ES_Teorico = round(es_teorico, 2)
  )
  
  return(resultados)
}

risco_95_he <- calcular_risco_teorico(mu, beta_he, 0.05);risco_95_he
risco_99_he <- calcular_risco_teorico(mu, beta_he, 0.01);risco_99_he
risco_95_ou <- calcular_risco_teorico(mu, beta_ou, 0.05);risco_95_ou
risco_99_ou <- calcular_risco_teorico(mu, beta_ou, 0.01);risco_99_ou

# Para GOU
beta_gou_var <- exp(beta_gou + (sigma_gou^2/(4*alpha_gou)));beta_gou_var

risco_95_gou <- calcular_risco_teorico(mu, beta_gou_var, 0.05);risco_95_gou
risco_99_gou <- calcular_risco_teorico(mu, beta_gou_var, 0.01);risco_99_gou
