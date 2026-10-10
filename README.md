# Codigo da dissertação de mestrado
Este repositório reúne o código em R desenvolvido no âmbito da dissertação de mestrado *Análise comparativa de estimadores de volatilidade: impacto no apreçamento de opções, na modelação estocástica e na medição de risco*.

## Organização dos ficheiros
A pasta `codigo` contém três ficheiros R:
- **01_estimadores_e_simulacao_gbm.R** — Comparação dos estimadores de volatilidade Histórico, Parkinson, Garman–Klass e Rogers–Satchell, aplicação ao modelo de Black–Scholes e simulações de movimento Browniano geométrico.
- **02_modelos_SV.R** — Estimação de parâmetros dos modelos de volatilidade estocástica Hull–White, Heston, Ornstein–Uhlenbeck e Ornstein–Uhlenbeck geométrico, simulação de trajetórias e análise de risco.
- **03_calibracao_heston.R** — Calibração do modelo de Heston a preços de opções da Tesla e comparação com o modelo de Black–Scholes.

A pasta `dados` contém o ficheiro `tesla_options.xlsx`, com os dados de opções e as taxas de juro utilizadas na calibração, incluindo os cálculos de interpolação.
