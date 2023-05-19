from scipy.signal import TransferFunction
import matplotlib.pyplot as plt
import numpy as np
import scipy

from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot
xi = 0.51
q1 = 1/1.84
q2 = 1/0.76

num = [1, 0, 0 ]
den = [1, 2.61, 3.41, 2.61, 1 ]


Hnor = TransferFunction(num, den)
num, den = scipy.signal.lp2lp(num, den, xi**(-1/4))

H = TransferFunction(num, den)

pzmap(Hnor, fig_id = 1) #S plane pole/zero plot
pzmap(H, fig_id = 2) #S plane pole/zero plot
print("H = ", H)
