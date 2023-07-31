# DIMENSÕES DO TROCADOR DE CALOR
import math

d_h = 2E-3  # DIÂMETRO DO CANAL QUENTE
t_h = 4E-3  # DISTÂNCIA ENTRE OS CANAIS QUENTES
d_c = 2E-3  # DIÂMETRO DO CANAL FRIO
t_c = 4E-3  # DISTÂNCIA ENTRE OS CANAIS FRIOS

w = 0.5E-3  # DISTÂNCIA ENTRE AS CAMADAS
n_h = 12  # NÚMERO DE CANAIS POR CAMADA QUENTE
N_h = 7  # NÚMERO DE CAMADAS QUENTE
n_c = 12  # NÚMERO DE CANAIS POR CAMADA FRIO
N_c = 7  # NÚMERO DE CAMADAS FRIO

L_h = t_h*n_h  # COMPRIMENTO DO CANAL QUENTE
L_c = t_c*n_c  # COMPRIMENTO DO CANAL FRIO
A = (N_h+N_c)*d_h + (N_h + N_c + 1)*w  # ALTURA

Ast_h = math.pi * (d_h/2)**2  # ÁREA DE S.T DO CANAL
Ast_c = math.pi * (d_c/2)**2  # ÁREA DE S.T DO CANAL

Am_h = d_h * math.pi * L_h  # ÁREA DE TROCA DE CALOR DO CANAL
Am_c = d_c * math.pi * L_c  # ÁREA DE TROCA DE CALOR DO CANAL

Atot_h = Am_h * n_h * N_h  # ÁREA DE TROCA DE CALOR TOTAL CANAL QUENTE
Atot_c = Am_c * n_c * N_c  # ÁREA DE TROCA DE CALOR TOTAL CANAL FRIO

V = A * L_h * L_c  # VOLUME DO NÚCLEO
Vch_h = Ast_h * L_h * n_h * N_h  # VOLUME DE TODOS OS CANAIS QUENTES
Vch_c = Ast_c * L_c * n_c * N_c  # VOLUME DE TODOS OS CANAIS FRIOS

Poro_h = Vch_h/V  # POROSIDADE CANAL QUENTE
Poro_c = Vch_c/V  # POROSIDADE CANAL FRIO
Afree_h = Poro_h * A * L_h  # ÁREA LIVRE DE ESCOAMENTO QUENTE
Afree_c = Poro_c * A * L_c  # ÁREA LIVRE DE ESCOAMENTO FRIO

beta = (Atot_h + Atot_c) / V  # GRAU DE COMPACTAÇÃO
k = 0.23  # CONDUTIVIDADE TÉRMICA DO MATERIAL DO TROCADOR
