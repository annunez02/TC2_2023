# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 18:38:34 2023

@author: Ana Nuñez

Simulacion para un valor determinado de R2/R1, R3 y C
"""

# Librerías externas NumPy, SciPy y Matplotlib
from scipy.signal import TransferFunction
import matplotlib.pyplot as plt
import numpy as np

from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot

R1 = 1
R2 = 1
R3 = 1000
C = 1e-6

w0 = 1/(R3*C)

plt.close('all')
    
my_tf = TransferFunction( [1, -(R2/R1)*w0], [1, w0] )
    
bodePlot(my_tf, fig_id=1)
    
pzmap(my_tf, fig_id=2) #S plane pole/zero plot
    
    
    