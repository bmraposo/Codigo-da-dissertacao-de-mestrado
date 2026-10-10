# Calibração Heston

#------------------------------------------------------
# FÓRMULA SEMIFECHADA DE HESTON
#------------------------------------------------------

# função caracteristica (f_j) 
heston_char_func <- function(phi, S0, K, r, T, V0, alpha, beta, sigma, rho, j) {
  
  # u_j e b_j 
  if (j == 1) {
    u <- 0.5    
    b <- alpha - rho * sigma
  } else {
    u <- -0.5
    b <- alpha
  }
  
  i <- 0 + 1i
  
  # d_j
  d <- sqrt((rho * sigma * phi * i - b)^2 - (sigma^2) * (2 * u * phi * i - phi^2))
  
  # g_j
  g <- (b - rho * sigma * phi * i + d) / (b - rho * sigma * phi * i - d)
  
  # D_j
  D_j <- ((b - rho * sigma * phi * i + d) / sigma^2) * ((1 - exp(d * T)) / (1 - g * exp(d * T)))
  
  # C_j
  C_j <- r * phi * i * T + (alpha * beta / sigma^2) * (
    (b - rho * sigma * phi * i + d) * T - 2 * log((1 - g * exp(d * T)) / (1 - g))
  )
  
  x <- log(S0)
  
  # Função Característica f_j
  f_j <- exp(C_j + D_j * V0 + i * phi * x)
  
  return(f_j)
}

#--------

# calculo do integrando de P_j
heston_integrand <- function(phi, S0, K, r, T, V0, alpha, beta, sigma, rho, j) {
  i <- 0 + 1i
  
  # Chamada da função característica
  f_j <- heston_char_func(phi, S0, K, r, T, V0, alpha, beta, sigma, rho, j)
  
  # Expressão de dentro do integral
  expr <- exp(-i * phi * log(K)) * f_j / (i * phi)
  
  return(Re(expr))
}

#--------

# função principal de apreçamento (C_0)
heston_call_price <- function(S0, K, r, T, V0, alpha, beta, sigma, rho) {
  
  limit <- 64 # para evitar overflow 
  
  # P_1 
  int_1 <- integrate(heston_integrand, lower = 1e-6, upper = limit, 
                     S0=S0, K=K, r=r, T=T, V0=V0, alpha=alpha, beta=beta, sigma=sigma, rho=rho
                     , j=1, subdivisions = 2000, stop.on.error = FALSE)
  P1 <- 0.5 + (1 / pi) * int_1$value
  
  # P_2 
  int_2 <- integrate(heston_integrand, lower = 1e-6, upper = limit, 
                     S0=S0, K=K, r=r, T=T, V0=V0, alpha=alpha, beta=beta, sigma=sigma, rho=rho
                     , j=2, subdivisions = 2000, stop.on.error = FALSE)
  P2 <- 0.5 + (1 / pi) * int_2$value
  
  C0 <- S0 * P1 - K * exp(-r * T) * P2
  
  return(C0)
}

#-------------------------------------------------------------------------------

# Calibração - Heston

# ------------------------------------------------------------------------------
# Preparação da tabela de preços de opções 
# ------------------------------------------------------------------------------

library(readxl)  

options_data <- read_excel(file.choose(), sheet = "Data_R")

names(options_data)[names(options_data) == "Time_to_Maturity_T"] <- "T"
names(options_data)[names(options_data) == "Risk_Free_Rate_r"] <- "r"

colunas_para_eliminar <- c("Last_Trade_Date", "Last_Price", "Change", "%_Change", 
                           "Volume", "Open_Interest", "Implied_Volatility", "Maturity_Date")

options_data <- options_data[, !(names(options_data) %in% colunas_para_eliminar)]
options_data$Option_Type <- as.factor(options_data$Option_Type)

View(options_data)

str(options_data)

S0 <- 402.90;S0 # preço de fecho do dia 7/7/2026 da Tesla

#------------

# Formula semi-fechada para tratar calls e puts
heston_price_universal <- function(Option_Type, S0, K, r, T, V0, alpha, beta, sigma, rho) {
  
  # Calcula sempre o preço da Call primeiro 
  C0 <- heston_call_price(S0, K, r, T, V0, alpha, beta, sigma, rho)
  
  if (Option_Type == "Put") {
    # Aplica a Paridade Put-Call
    P0 <- C0 - S0 + K * exp(-r * T)
    return(P0)
  } else {
    return(C0)
  }
}

# parametros iniciais para calibração
rho_he # do ficheiro 02_modelos_SV
params_init <- c(V0 = 0.25, alpha = 2, beta = 0.25, 
                 sigma = 0.5, rho = rho_he);params_init 

#-----------

heston_sse <- function(params, data, S0) {
  
  V0    <- params[1]
  alpha <- params[2]
  beta  <- params[3]
  sigma <- params[4]
  rho   <- params[5]
  
  model_prices <- mapply(FUN = heston_price_universal,
                         Option_Type = as.character(data$Option_Type),
                         K = data$Strike,
                         r = data$r,
                         T = data$T,
                         MoreArgs = list(S0 = S0, V0 = V0, alpha = alpha, 
                                         beta = beta, sigma = sigma, rho = rho))
  # Cálculo do SSE
  sse <- sum((data$Market_Price - model_prices)^2)
  
  return(sse)
}

lower_bounds <- c(1e-4, 1e-4, 1e-4, 1e-4, -0.99)
upper_bounds <- c(10.0, 300.0, 20.0, 20.0, 0.99)

#-------------

# Execução da calibração

fit_heston <- optim(
  par = params_init,                
  fn = heston_sse,                  
  data = options_data,              
  S0 = S0,                          
  method = "L-BFGS-B",              
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(trace = 1,          
                 maxit = 1000)      
)

# Resultados
print(fit_heston$par)
sse_final <- fit_heston$value;sse_final
N <- nrow(options_data);N
rmse_final <- sqrt(sse_final/N);rmse_final

cat("SSE:", fit_heston$value, "\n")
cat("RMSE: $", round(rmse_final, 4), "\n")

#-----------------

# Comparação: preços calibrados vs preços de mercado

library(ggplot2)

fit_heston$par

V0_cal    <- fit_heston$par[1]
alpha_cal <- fit_heston$par[2]
beta_cal  <- fit_heston$par[3]
sigma_cal <- fit_heston$par[4]
rho_cal   <- fit_heston$par[5]

options_data$Model_Price <- mapply(
  FUN = heston_price_universal,
  Option_Type = as.character(options_data$Option_Type),
  K = options_data$Strike,
  r = options_data$r,
  T = options_data$T,
  MoreArgs = list(S0 = S0, V0 = V0_cal, alpha = alpha_cal, 
                  beta = beta_cal, sigma = sigma_cal, rho = rho_cal)
)

fit_heston$value

SST <- sum(
  (options_data$Market_Price -
     mean(options_data$Market_Price))^2
)

R2_calibracao <- 1 - fit_heston$value / SST;R2_calibracao

plot_fit <- ggplot(options_data, aes(x = Market_Price, y = Model_Price, color = Option_Type)) +
  geom_point(alpha = 0.7, size = 2.5) +
  
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Preço de Mercado (USD)",
    y = "Preço Estimado pelo Modelo (USD)",
    color = "Tipo de Opção"
  ) +
  scale_color_manual(values = c("Call" = "#1b9e77", "Put" = "#d95f02")) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12, color = "darkgray"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(plot_fit)

#-----------------------

# Avaliar erro: Model_Price - Market_Price

options_data$Residual <- options_data$Model_Price - options_data$Market_Price

plot_residuals <- ggplot(options_data, aes(x = Strike, y = Residual, color = Option_Type)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Preço de Exercício (USD)",
    y = "Resíduo (USD)",
    color = "Tipo de Opção"
  ) +
  scale_color_manual(values = c("Call" = "#1b9e77", "Put" = "#d95f02")) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12, color = "darkgray"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(plot_residuals)

#-----------------

# Spread medio da amostra 

spreads_absolutos <- options_data$Ask - options_data$Bid
mean_spread_abs <- mean(spreads_absolutos, na.rm = TRUE)
cat("Spread Médio: $", round(mean_spread_abs, 2), "\n")

#----------

# Preços Calls, Heston calibrado

# parâmetros 
S0 <- 402.90 
strikes <- c(0.97 * S0, S0, 1.03 * S0);strikes
r <- 0.0406  
T_mat <- 1

V0_cal    
alpha_cal 
beta_cal  
sigma_cal 
rho_cal

p_he_cal_ITM <- heston_call_price(S0, strikes[1], r, T_mat, V0_cal, alpha_cal, beta_cal, sigma_cal, rho_cal);p_he_cal_ITM
p_he_cal_ATM <- heston_call_price(S0, strikes[2], r, T_mat, V0_cal, alpha_cal, beta_cal, sigma_cal, rho_cal);p_he_cal_ATM
p_he_cal_OTM <- heston_call_price(S0, strikes[3], r, T_mat, V0_cal, alpha_cal, beta_cal, sigma_cal, rho_cal);p_he_cal_OTM

#--------------------

# Comparar apreçamento de Heston e Black-Scholes

library(tidyquant)
library(dplyr)
library(ggplot2)

calculo_vol <- function(ticker, from_date, to_date) { 

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
  
  sigma_historica <- sd(dados$log_ret) * sqrt(252)
  
  const_park <- 1 / (4 * log(2)) 
  D_l <- (const_park / nrow(dados)) * sum(dados$li^2)
  sigma_parkinson <- sqrt(D_l) * sqrt(252)
  
  termos_gk <- 0.5 * (dados$u - dados$d)^2 - (2 * log(2) - 1) * dados$c^2
  var_media_gk <- mean(termos_gk)
  sigma_gk <- sqrt(var_media_gk) * sqrt(252)
  
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

tsla_vol_heston <- calculo_vol("TSLA", from_date = "2025-07-07", to_date = "2026-07-08"); tsla_vol_heston 

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

precos_calls(S0,tsla_vol_heston$vol_parkinson,r = 0.0406) 
precos_calls(S0,tsla_vol_heston$vol_gk,r = 0.0406) 
precos_calls(S0,tsla_vol_heston$vol_rs,r = 0.0406) 

p_he_cal_ITM
p_he_cal_ATM
p_he_cal_OTM