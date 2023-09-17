# Librerías externas NumPy, SciPy y Matplotlib
import scipy.signal as ss
import matplotlib.pyplot as plt
import numpy as np

from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot, analyze_sys, pretty_print_SOS

# Definicion de parámetros
fs1 = 10000
fs2 = 10
n = 2

# Obtencion de la transferencia butter normalizada
zz, pp, kk = ss.buttap(N=n)

num, den = ss.zpk2tf(zz, pp, kk)
TF = ss.TransferFunction(num, den)

# Aplicacion de la transformada bilineal para la obtencion del filtro digital

z, p, k = ss.bilinear_zpk(zz, pp, kk, fs1)

numz , denz = ss.zpk2tf(z, p, k)
TFZ = ss.TransferFunction(numz, denz)
#w, HTFZ = ss.freqz(b=numz, a=denz, fs=fs2)

# Repito para fs = 10

z_, p_, k_ = ss.bilinear_zpk(zz, pp, kk, fs2)

numz_ , denz_ = ss.zpk2tf(z_, p_, k_)
TFZ_ = ss.TransferFunction(numz_, denz_)
#w_, HTFZ_ = ss.freqz(b=numz_, a=denz_, fs=fs2)

# Visualizo
analyze_sys(TFZ, digital = True, sys_name='TFZ fs = 1K')
analyze_sys(TFZ_, digital = True, sys_name='TFZ fs = 10')