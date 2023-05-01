
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


# Definicion de los valores para la simulacion de la transferencia 
R = 848.83
C = 100e-9
L = 72.05e-3


# Definicion de la transferencia normalizada
T_nor = TransferFunction( [w0**3], [1, 2*w0, 2*(w0**2), w0**3] )


# Definicion de la transferencia 
num_RC = [1/(R*C)]
den_RC = [1, 1/(R*C)]

num_RLC = [1/(L*C)]
den_RLC = [1, R/L, 1/(L*C)]

num_T = np.convolve(num_RC, num_RLC)
den_T = np.convolve(den_RC, den_RLC)

T = TransferFunction(num_T ,den_T)    

# Visualizacion

plt.close('all')
   
bodePlot(T_nor, fig_id = 1)
pzmap(T_nor, fig_id = 2) #S plane pole/zero plot

# Nota: al graficar la transferencia desnormalizada me queda mal la escala de magnitud y no logra apreciarse

bodePlot(T, fig_id = 3)    
pzmap(T, fig_id = 4) #S plane pole/zero plot
    
    
    