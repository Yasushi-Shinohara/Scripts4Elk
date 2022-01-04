# coding: UTF-8
# This is created 2021/12/27 by Y. Shinohara
# This is lastly modified 2021/xx/yy by Y. Shinohara
import sys
from modules.constants import pi, fpi, zI, Hartree
import numpy as np

class GenerateSigmaEpsilon:
    def __init__(self):
        self.omega = None
        self.sigma = None #Optical Conductivity
        self.epsilon = None #Dielectric function

    @classmethod
    def generate(self, ED, Nene = 2000, ewidth = 0.01, plot_option = True):
        emax = 1.2*(np.amax(ED.eigval) - np.amin(ED.eigval))
        emin = -1.0*emax
        self.omega = np.linspace(emin, emax, Nene)
        print(emin, emax, self.omega)
        self.sigma = np.zeros([Nene, 3, 3], dtype='complex128')
        Nevery = int(ED.Nk/10.0)
        print('# Following is progress of GenerateSigmaEpsilon.generate function. ')
        for ik in range(ED.Nk):
            if (ik%Nevery ==0):
                print('# '+str(np.round(ik/ED.Nk, decimals = 2)*100.0)+' % is done.')
            for ib in range(ED.Nb):
                for jb in range(ED.Nb):
                    if (ib != jb):
                        ene_denominator = self.omega[:] - (ED.eigval[jb,ik] - ED.eigval[ib,ik]) + zI*ewidth
                        ene_denominator = ene_denominator*(ED.eigval[jb,ik] - ED.eigval[ib,ik])
                        docc = ED.occ[ib,ik] - ED.occ[jb,ik]
                        for ixyz in range(3):
                            for jxyz in range(3):
                                moment = ED.pmat[ixyz,ib,jb,ik]*ED.pmat[jxyz,jb,ib,ik]
                                self.sigma[:,ixyz,jxyz] = self.sigma[:,ixyz,jxyz] + docc*moment/ene_denominator[:]
        self.sigma = self.sigma*zI/ED.vcell/ED.Nk
        #Constructing epsilon from sigma
        self.epsilon = np.zeros([Nene, 3, 3], dtype='complex128')
        for ixyz in range(3):
            for jxyz in range(3):
                self.epsilon[:,ixyz,jxyz] = self.sigma[:,ixyz,jxyz]*fpi*zI/self.omega
        for ixyz in range(3):
            self.epsilon[:,ixyz,ixyz] = np.ones(Nene, dtype='complex128') + self.epsilon[:,ixyz,ixyz]
        print('# Number of energy grid: Nene =', Nene)
        print('# Energy width for sigma-epsilon: ewidth =', ewidth, '[a.u.] =', ewidth*Hartree, '[eV]')
        print('# Energy minimum and maximum for DoS: emin, emax =', emin, emax, '[a.u.] =', emin*Hartree, emax*Hartree, '[eV]')
#        print('# Number of integrad DoS:', self.NoS[NDoS-1])
        if (plot_option):
            import matplotlib.pyplot as plt
            plt.figure()
            plt.title('Optical conductivity')
            plt.plot(self.omega, np.real(self.sigma[:,0,0]), label='Real part')
            plt.plot(self.omega, np.imag(self.sigma[:,0,0]), label='Imaginary part')
            plt.grid()
            plt.legend()
            plt.show()
            #
            plt.figure()
            plt.title('Dielectrc function')
            plt.plot(self.omega, np.real(self.epsilon[:,0,0]), label='Real part')
            plt.plot(self.omega, np.imag(self.epsilon[:,0,0]), label='Imaginary part')
            plt.grid()
            plt.legend()
            plt.show()
        return self.omega, self.sigma, self.epsilon
