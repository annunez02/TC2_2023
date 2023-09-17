# Inicialización e importación de módulos

# Módulos para Jupyter
import warnings
warnings.filterwarnings('ignore')

# Módulos importantantes
import scipy.signal as sig
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.io as sio
from pytc2.sistemas_lineales import plot_plantilla

# -----------------------------------------ECG---------------------------------------------

fig_sz_x = 10
fig_sz_y = 7
fig_dpi = 100 # dpi

fig_font_size = 16

mpl.rcParams['figure.figsize'] = (fig_sz_x,fig_sz_y)
plt.rcParams.update({'font.size':fig_font_size})

###
## Señal de ECG registrada a 1 kHz, con contaminación de diversos orígenes.
###

# para listar las variables que hay en el archivo
#io.whosmat('ecg.mat')
mat_struct = sio.loadmat('ecg.mat')

ecg_one_lead = mat_struct['ecg_lead']
ecg_one_lead = ecg_one_lead.flatten()
cant_muestras = len(ecg_one_lead)


# --------------------------------Definicion de parámetros---------------------------------

fs = 1000 # Hz
nyq_frec = fs / 2

# Plantilla

# filter design
ripple = 0.1 # dB
atenuacion = 40 # dB

ws1 = 1.0 #Hz
wp1 = 3.0 #Hz
wp2 = 25.0 #Hz
ws2 = 35.0 #Hz

# A partir de aca empiezo a crear

# -----------------------------------------IIR---------------------------------------------

system = sig.iirdesign(wp = [wp1, wp2], ws = [ws1, ws2], gpass = ripple, gstop = atenuacion, analog = False, ftype = 'butter', output = 'sos', fs = 1000)

print(f"El filtro tiene {system.shape[0]} SOS")

w, TF_IIR = sig.sosfreqz(system, worN = 2000, fs = fs)

phase = np.angle(TF_IIR)

num_IIR, den_IIR = sig.sos2tf(system)
w, gd = sig.group_delay((num_IIR, den_IIR), w = 2000, whole = False, fs = fs)

# ----- Modulo -----

plt.figure(1)

plt.plot(w, 20 * np.log10(abs(TF_IIR)), label = 'IIR {:d} SOS'.format(system.shape[0]))

plt.title('Filtro IIR')
plt.xlabel('Frecuencia [Hz]')
plt.ylabel('Módulo [dB]')
plt.grid()
plt.axis([0, 100, -100, 5 ])

plot_plantilla(filter_type = 'bandpass', fpass = [wp1, wp2], ripple = ripple , fstop = [ws1, ws2], attenuation = atenuacion, fs = fs)

# ----- Fase -----

plt.figure(2)

plt.plot(w, phase)
plt.title('Fase del filtro IIR')
plt.xlabel('Frecuencia [Hz]')
plt.ylabel('Fase [rad]')
plt.grid()
plt.axis([0, 100, -2*np.pi, 2*np.pi ])

# ----- Group delay -----

plt.figure(3)

plt.plot(w, gd)
plt.title('Retardo de grupo del filtro IIR')
plt.ylabel('Group delay [samples]')
plt.xlabel('Frequency [rad/sample]')
plt.grid()

# -----------------------------------------FIR---------------------------------------------

#-----Definicion de parámetros-----

frecs = np.array([        0.0,          ws1,         wp1,        wp2,           ws2,         nyq_frec   ]) / nyq_frec
gains = np.array([-atenuacion,  -atenuacion,     -ripple,    -ripple,   -atenuacion,        -atenuacion ])
gains = 10**(gains/20)

#-----Obtengo la transferencia-----

cant_coef = 5501

num_win = sig.firwin2(cant_coef, frecs, gains , window='hamming')
den = 1.0

# Muestreo el filtro donde me interesa verlo según la plantilla

w  = np.append(np.logspace(-1, 0.8, 250), np.logspace(0.9, 1.6, 250) )
w  = np.append(w, np.linspace(110, nyq_frec, 100, endpoint=True) ) / nyq_frec * np.pi

_, TF_FIR = sig.freqz(num_win, den, w)

# Renormalizo el eje de frecuencia

w = w / np.pi * nyq_frec

phase = np.angle(TF_FIR)
wg, gd = sig.group_delay((num_win, den), w = 2000, whole = False, fs = fs)

# ----- Modulo -----

plt.figure(4)

plt.plot(w, 20 * np.log10(abs(TF_FIR)), label='FIR-Win {:d}'.format(num_win.shape[0]))

plt.title('Filtro FIR')
plt.xlabel('Frecuencia [Hz]')
plt.ylabel('Módulo [dB]')
plt.grid()
plt.axis([0, 100, -60, 5 ])

axes_hdl = plt.gca()
axes_hdl.legend()

plot_plantilla(filter_type = 'bandpass', fpass = [wp1, wp2], ripple = ripple , fstop = [ws1, ws2], attenuation = atenuacion, fs = fs)


# ----- Fase -----

plt.figure(5)

plt.plot(w / np.pi * nyq_frec, phase)
plt.title('Fase del filtro FIR')
plt.xlabel('Frecuencia [Hz]')
plt.ylabel('Fase [rad]')
plt.grid()
plt.axis([0, nyq_frec, -2*np.pi, 2*np.pi ])

# ----- Group delay -----

plt.figure(3)

plt.plot(wg, gd)
plt.title('Retardo de grupo del filtro FIR')
plt.ylabel('Group delay [samples]')
plt.xlabel('Frequency [rad/sample]')
plt.grid()
plt.axis([0, nyq_frec, 0, 5000])


 