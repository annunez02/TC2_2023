import scipy.signal as sig
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.io as sio
from pytc2.sistemas_lineales import plot_plantilla

fs = 40.0   # KHz

wp1 = 2.0   # KHz
ws1 = 4.0   # KHz
ws2 = 6.0   # KHz
wp2 = 8.0   # KHz

ripple = 1.0        # dB
atenuacion = 20.0   # dB
nyq_frec = fs/2     # KHz

frecs = np.array([    0.0,      wp1,             ws1,            ws2,       wp2,      nyq_frec]) / nyq_frec
gains = np.array([-ripple,  -ripple,     -atenuacion,    -atenuacion,   -ripple,      -ripple ])
gains = 10**(gains/20)

num = sig.firls(5001, bands=frecs, desired=gains, weight=None, fs=fs)
den = 1.0

w  = np.append(np.logspace(-1, 0.8, 300), np.logspace(0.9, 1.6, 300) )
w  = np.append(w, np.linspace(110, nyq_frec, 100, endpoint=True) ) / nyq_frec * np.pi


_, TF_FIR = sig.freqz(num, den, w)

# Renormalizo el eje de frecuencia

w = w / np.pi * nyq_frec

phase = np.angle(TF_FIR)
wg, gd = sig.group_delay((num, den), w = 2000, whole = False, fs = fs)

# ----- Modulo -----

plt.figure()

plt.plot(w, 20 * np.log10(abs(TF_FIR)), label='FIR-Win {:d}'.format(num.shape[0]))

plt.title('Filtro FIR')
plt.xlabel('Frecuencia [Hz]')
plt.ylabel('Módulo [dB]')
plt.grid()
plt.axis([0, 20, -30, 20 ])

axes_hdl = plt.gca()
axes_hdl.legend()

plot_plantilla(filter_type = 'bandstop', fpass = [wp1, wp2], ripple = ripple , fstop = [ws1, ws2], attenuation = atenuacion, fs = fs)
