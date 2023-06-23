from scipy.signal import TransferFunction
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as ss

from pytc2.sistemas_lineales import analyze_sys, pretty_print_bicuad_omegayq, tf2sos_analog, pretty_print_SOS


#----------------------Definicion de las variables---------------------
#---- Para el filtro pasa-altos ----
w0 = 1
wze = 1/3

#---- Para el filtro pasa-bajos ----
W0 = 1/w0
Wze = 1/wze

#---- Z P K LP ----
[r1, r2] = np.roots([1, 1, 1])

zLP = [3j, -3j]
pLP = [-1, r1, r2]
kLP = 1/9

#---------------------------Definicion de HLP--------------------------

[numLP, denLP] = ss.zpk2tf(zLP, pLP, kLP)

LP_sos = tf2sos_analog(numLP, denLP)

#-------------------------Transformación LP-HP-------------------------

[numHP, denHP] = ss.lp2hp(numLP, denLP, w0)

[zHP, pHP, kHP] = ss.tf2zpk(numHP, denHP)

HP_sos = tf2sos_analog(numHP, denHP)

#-----------------------------Visualizacion----------------------------

# -------- LP --------
print("HLP(s): ")
LP_sos[LP_sos < 1e-6] = 0.0
pretty_print_SOS(LP_sos, mode='omegayq')
#print("zLP = ", zLP, "\n", "pLP = ", pLP, "\n", "kLP =", kLP)
analyze_sys(LP_sos, "LP prototype")

#print("\n")

# -------- HP --------
print("HHP(s): ")
HP_sos[HP_sos < 1e-6] = 0.0
pretty_print_SOS(HP_sos, mode='omegayq')
analyze_sys(HP_sos, "HP transformation", same_figs=False)