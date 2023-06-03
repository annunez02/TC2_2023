from scipy.signal import TransferFunction
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as ss

from pytc2.sistemas_lineales import analyze_sys, pretty_print_bicuad_omegayq, tf2sos_analog, pretty_print_SOS

from pytc2.general import print_subtitle

#----------------------Definicion de las variables---------------------

xi = 0.35
n = 3
rp = 0.5

#---------------------------Definicion de H-------------------------

[z, p, k] = ss.cheb1ap(n, rp)

[num, den] = ss.zpk2tf(z, p, k)

this_sos = tf2sos_analog(num, den)

Hnor = TransferFunction(num, den)

#-----------------------------Visualizacion----------------------------

# pzmap(Hnor, fig_id = 1) #S plane pole/zero plot
# bodePlot(Hnor, fig_id = 2)
# GroupDelay(Hnor, fig_id = 3)
pretty_print_SOS(this_sos, mode='omegayq')
print("z = ", z, "\n", "p = ", p, "\n", "k =", k)
analyze_sys(Hnor, "Cheby")
  
