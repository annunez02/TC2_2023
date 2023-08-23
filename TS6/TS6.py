import sympy as sp
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as ss

from pytc2.dibujar import dibujar_Pi, dibujar_Tee, dibujar_lattice
from pytc2.general import print_latex, print_subtitle, a_equal_b_latex_s
from pytc2.general import print_latex

from pytc2.cuadripolos import calc_MAI_impedance_ij, calc_MAI_vtransf_ij_mn, calc_MAI_ztransf_ij_mn

#----------------------Definicion de las variables---------------------
S = sp.symbols('S', complex=True)
R = 1
L1 = 1.5
C2 = 4/3
L3= 0.5

Z_C2, Z_L1, Z_L3 = sp.symbols('Z_C2 Z_L1 Z_L3', real=True)
R1 = sp.symbols('R1', real=True, positive=True)

#---------------------------Definicion de matrices--------------------------

Ta_sym = sp.Matrix([  
                    [ Z_L1 * 1/Z_C2 + 1,     Z_L1],
                    [ 1/Z_C2,                   1]
                ])


Tb_sym = sp.Matrix([  
                    [ Z_L3 * 1/R1 + 1,     Z_L3],
                    [ 1/R1,                   1]
                ])

#---------------------------Definicion de matrices--------------------------

Ta = sp.Matrix([  
                    [ S * L1 * (S * C2) + 1,    S * L1],
                    [ (S * C2),                   1]
                ])


Tb = sp.Matrix([  
                    [ S * L3 * 1/R + 1,     S * L3],
                    [ 1/R,                   1]
                ])

#-------------------------Calculo de la matriz resultante-------------------------

Ttotal_sym = Ta_sym * Tb_sym
Ttotal_sym = sp.simplify(Ttotal_sym)

Ttotal = Ta * Tb
Ttotal = sp.simplify(Ttotal)
# con_detalles = False
con_detalles = True

#-----------------------------Visualizacion----------------------------
print_subtitle('Ta')
print_latex(a_equal_b_latex_s('T_a', Ta))

print_subtitle('Tb')
print_latex(a_equal_b_latex_s('T_b', Tb))

print_subtitle('Ttotal_sym')
print_latex(a_equal_b_latex_s('T_{tot}', Ttotal_sym))

print_subtitle('Ttotal')
print_latex(a_equal_b_latex_s('T_{tot}', Ttotal))



#################################### MAI #########################################

Ymai = sp.Matrix([  
                    [ 1/(S * Z_L1),           0,              0,             -1/(S * Z_L1)],
                    [ 0,            1/(S * Z_L3) + 1/R1,    -1/R1,          -1/(S * Z_L3)],
                    [ 0,                    -1/R1,      1/R1 + S * Z_C2,     -S * Z_C2],
                    [ -1/(S * Z_L1),      -1/(S * Z_L3),     -S * Z_C2,      1/(S * Z_L1) + 1/(S * Z_L3) + S * Z_C2 ]
                 ])

# con_detalles = False
con_detalles = True

# Calculo la transf a partir de la MAI

print('Transferencia de tensión:')
Vmai = calc_MAI_vtransf_ij_mn(Ymai, 1, 2, 0, 2, verbose=con_detalles)
Vmai = sp.simplify(Vmai)
print_latex( r'T^{{ {:d}{:d} }}_{{ {:d}{:d} }} = '.format(1, 2, 0, 2) +  sp.latex(Vmai))

# Calculo la Z en el puerto de entrada a partir de la MAI
Zmai = calc_MAI_impedance_ij(Ymai, 0, 2, verbose=con_detalles)
print_latex( r'Z_{{ {:d}{:d} }} = '.format(0,2) +  sp.latex(Zmai) )

