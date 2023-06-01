from scipy.signal import TransferFunction
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as ss

from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot

#----------------------Definicion de las variables---------------------

xi = 1
BW = 9/20
k_10 = 3.162

#---------------------------Definicion de H_LP-------------------------

num = [1]
den = [1, 2, 2, 1 ]

Hnor = TransferFunction(num, den)

#-------------------------Transformacion LP a BP-----------------------

num, den = ss.lp2bp(num, den, xi , BW)
z, p, k = ss.tf2zpk(num, den)

#----Multiplico por la constante para lograr una ganancia de 10 dB-----

H = TransferFunction(k_10 * num, den)

#-----------------------------Visualizacion----------------------------

pzmap(Hnor, fig_id = 1) #S plane pole/zero plot
pzmap(H, fig_id = 2) #S plane pole/zero plot
bodePlot(H, fig_id = 3)
GroupDelay(H, fig_id = 4)
print("num = ", num, "\n",  "den = ", den)
print("z = ", z, "\n", "p = ", p, "\n", "k =", k)
  
