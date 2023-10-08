import scipy.signal as sig
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.io as sio
import pandas   as pd
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

num = sig.firls(201, bands=frecs, desired=gains, weight=None, fs=fs)
print(f"num: {num}")
df = pd.DataFrame({'num': num})
df.to_csv("FIR_LSQR.csv", index=False, encoding='utf-8-sig')
den = 1.0

w, TF_FIR = sig.freqz(num, den, fs=fs)

# Renormalizo el eje de frecuencia
w_hz = w / (2 * np.pi)

# Plot the magnitude response in dB
plt.plot(w_hz, 20 * np.log10(np.abs(TF_FIR)))
plt.grid(True)
plt.ylim(-50, 10)
plt.xlabel('Normalized Frequency (Hz)')
plt.ylabel('Magnitude (dB)')
plt.show()