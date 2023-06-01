# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 18:38:34 2023

@author: Ana Nuñez

Simulacion de transferencia normalizada en frecuencia para distintos valores de R2/R1
"""

# Librerías externas NumPy, SciPy y Matplotlib
from scipy.signal import TransferFunction
import matplotlib.pyplot as plt
import numpy as np

from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot


w0 = 1

plt.close('all')

R1 = 1
R2 = [ 0.5, 1, 2]

for qq in range(len(R2)):
    
    my_tf = TransferFunction( [1, -(R2[qq]/R1)*w0], [1, w0] )
    
    bodePlot(my_tf, fig_id=1, filter_description = 'R2/R1={:3.3f}'.format(R2[qq]/R1) )
    
    pzmap(my_tf, fig_id=2, filter_description = 'R2/R1={:3.3f}'.format(R2[qq]/R1)) #S plane pole/zero plot
    
    
    
    