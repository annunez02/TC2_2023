
"""
@author: Ana Nuñez
Simulacion Tarea semanal 3
"""

# Librerías externas NumPy, SciPy y Matplotlib
from scipy.signal import TransferFunction
import matplotlib.pyplot as plt
import numpy as np

from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot

# Definicion de los valores para la simulacion de la transferencia normalizada
xi = 0.51
Q = 1
w0 = 1/(xi**(1/3))

# Definicion de la transferencia normalizada
T_nor = TransferFunction( [w0**3], [1, 2*w0, 2*(w0**2), w0**3] )


# Visualizacion

plt.close('all')
   
bodePlot(T_nor, fig_id = 1)
pzmap(T_nor, fig_id = 2) #S plane pole/zero plot
  
    