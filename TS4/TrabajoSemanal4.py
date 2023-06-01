
"""
@author: Ana Nuñez
Simulacion Tarea semanal 4
"""

# Librerías externas NumPy, SciPy y Matplotlib
import scipy.signal as ss
import matplotlib.pyplot as plt
import numpy as np

from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot, analyze_sys

# Definicion de los valores para la simulacion de la transferencia normalizada
xi = 0.51
Q = 1
w0 = 1/(xi**(1/3))

all_sys = []

# Definicion de la función transferencia normalizada pasa bajos

num_nor = [1]
den_nor = [1, 2, 2, 1]

H_norLP = ss.TransferFunction(num_nor, den_nor)

# Definicion de la función transferencia pasa bajos 

[numLP, denLP] = ss.lp2lp(num_nor, den_nor, w0)
H_LP = ss.TransferFunction(numLP, denLP)

# Definicion de la función transferencia pasa altos

[numHP, denHP] = ss.lp2hp(numLP, denLP)
H_HP = ss.TransferFunction(numHP, denHP)


all_sys.append(H_norLP)
all_sys.append(H_LP)
all_sys.append(H_HP)
# Visualizacion

plt.close('all')

analyze_sys(all_sys, ['H_norLP', 'H_LP', 'H_HP'])
   

  
#analyze_sys(H_HP, sys_name='pasa altos tercer orden Q={:d}'.format(Q))  
  
    