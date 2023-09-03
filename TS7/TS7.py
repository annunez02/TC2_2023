# Librerías externas NumPy, SciPy y Matplotlib
import scipy.signal as ss
import matplotlib.pyplot as plt
import numpy as np

from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot, analyze_sys

# Definicion de los valores para la simulacion de la transferencia normalizada
Q = 1/np.sqrt(2)
w0 = 0.01

all_sys = []

# Definicion de la función transferencia normalizada pasa bajos

num_nor = [1, 0, 1.25]
den_nor = [1, 1/Q, 1]

H_norLP = ss.TransferFunction(num_nor, den_nor)

# Definicion de la función transferencia pasa bajos 

[numLP, denLP] = ss.butter(2, 1, btype='low', analog=False, output='ba', fs=100)
H_LP = ss.TransferFunction(numLP, denLP)


all_sys.append(H_norLP)
all_sys.append(H_LP)
# Visualizacion

plt.close('all')

analyze_sys(all_sys, ['H_norLP', 'H_LP'])
   

  
#analyze_sys(H_HP, sys_name='pasa altos tercer orden Q={:d}'.format(Q))  