import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as ss

from pytc2.sistemas_lineales import analyze_sys, tf2sos_analog, pretty_print_SOS

#----------------------Definicion de las variables---------------------

xi = 0.35
n = 3
rp = 0.5

#---------------------------Definicion de H-------------------------

[z, p, k] = ss.cheb1ap(n, rp)

[num, den] = ss.zpk2tf(z, p, k)

sos_lp = tf2sos_analog(num, den)
sos_lp[sos_lp < 1e-6] = 0.0

#-----------------------------Visualizacion----------------------------

pretty_print_SOS(sos_lp, mode='omegayq')
analyze_sys(sos_lp, "Cheby")
  
