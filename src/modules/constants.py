# coding: UTF-8
# This is created 2020/04/17 by Y. Shinohara
# This is lastly modified 2020/04/20 by Y. Shinohara
import numpy as np
#Mathematicl constants
pi  = np.pi
tpi = 2.0*np.pi
fpi = 4.0*np.pi
zI  = 1.0j
#Physical constants
#https://ja.wikipedia.org/wiki/%E5%8E%9F%E5%AD%90%E5%8D%98%E4%BD%8D%E7%B3%BB
aB = 0.0529177210903 #nanometer
Hartree = 27.211386245988 #eV
Atomtime = 0.024188843265857 #fs
Atomfield = Hartree/aB #V/nm
Atomvolume = aB**3 #nm^3
#https://ja.wikipedia.org/wiki/%E5%BE%AE%E7%B4%B0%E6%A7%8B%E9%80%A0%E5%AE%9A%E6%95%B0
sol = 137.035999084 #speed of light
ch = 1241.5 #eV * nm
chbar = 197.3 # eV * nm
halfepsc = 3.509e16 # W/cm^2 \frac{1}{2}*\epsilon_0 * c
Atomfluence = halfepsc*Atomtime*1.0e-15 # J/cm^2 ,W/cm^2 * fs = femto J/cm^2
