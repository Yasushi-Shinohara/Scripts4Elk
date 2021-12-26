# coding: UTF-8
# This is created 2021/12/24 by Y. Shinohara
# This is lastly modified 2021/xx/yy by Y. Shinohara
import sys
from modules.constants import pi, Hartree
import numpy as np

class GenerateDoS:
    def __init__(self):
        self.DoS = None #Density of States
        self.omega = None
        self.NoS = None #Number of States
        self.occDoS = None #Occupied Density of States
        self.occNoS = None #Occupied Number of States

    @classmethod
    def generate(self, ED, NDoS = 2000, ewidth = 0.01, plot_option = True):
        emin = np.amin(ED.eigval) - 0.2*(np.amax(ED.eigval) - np.amin(ED.eigval))
        emax = np.amax(ED.eigval) + 0.2*(np.amax(ED.eigval) - np.amin(ED.eigval))
        self.omega = np.linspace(emin, emax, NDoS)
        print(emin, emax, self.omega)
        self.DoS = np.zeros(NDoS, dtype='float64')
        self.NoS = np.zeros(NDoS, dtype='float64')
        self.occDoS = np.zeros(NDoS, dtype='float64')
        self.occNoS = np.zeros(NDoS, dtype='float64')
        for ik in range(ED.Nk):
            for ib in range(ED.Nb):
                self.DoS[:] = self.DoS[:] + np.exp(-(self.omega[:]-ED.eigval[ib,ik])**2/ewidth**2)
                self.occDoS[:] = self.occDoS[:] + np.exp(-(self.omega[:]-ED.eigval[ib,ik])**2/ewidth**2)*ED.occ[ib,ik]
        self.DoS = self.DoS/np.sqrt(pi)/ewidth/ED.Nk*ED.spin_degeneracy
        self.occDoS = self.occDoS/np.sqrt(pi)/ewidth/ED.Nk
        domega = self.omega[1] - self.omega[0]
        for ie in range(NDoS-1):
            self.NoS[ie+1] = self.NoS[ie] + 0.5*(self.DoS[ie+1] + self.DoS[ie])
            self.occNoS[ie+1] = self.occNoS[ie] + 0.5*(self.occDoS[ie+1] + self.occDoS[ie])
        self.NoS = domega*self.NoS
        self.occNoS = domega*self.occNoS
        print('# Number of energy grid: NDoS =', NDoS)
        print('# Energy width for DoS: ewidth =', ewidth, '[a.u.] =', ewidth*Hartree, '[eV]')
        print('# Energy minimum and maximum for DoS: emin, emax =', emin, emax, '[a.u.] =', emin*Hartree, emax*Hartree, '[eV]')
        print('# Number of integrad DoS:', self.NoS[NDoS-1])
        print(self.NoS[NDoS-1], ED.Nb)
        print(self.occNoS[NDoS-1])
        if (plot_option):
            import matplotlib.pyplot as plt
            plt.figure()
            plt.title('Density of States')
            plt.fill_between(self.omega, 0.0*self.NoS, self.NoS/np.amax(self.NoS)*np.amax(self.DoS)*0.8, facecolor = 'k', alpha=0.25, label='NoS(normalized)')
            plt.plot(self.omega, self.DoS, label='DoS')
            plt.plot(self.omega, self.occDoS, label='occDoS')
            plt.grid()
            plt.legend()
            plt.show()
        return self.omega, self.DoS, self.NoS, self.occDoS, self.occNoS
