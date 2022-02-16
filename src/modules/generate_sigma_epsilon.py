# coding: UTF-8
# This is created 2021/12/27 by Y. Shinohara
# This is lastly modified 2021/xx/yy by Y. Shinohara
import sys
import os
from modules.constants import pi, fpi, zI, Hartree, Atomvolume
import numpy as np
import time
import ctypes as ct

class GenerateSigmaEpsilon:
    def __init__(self):
        self.Nene = None
        self.omega = None
        self.sigma = None #Optical Conductivity
        self.epsilon = None #Dielectric function
        self.epsilon_inv = None #Inverse dielectric function
        self.ewidth = None
        self.sum_epsilon = None
        self.sum_epsilon_inv = None
        self.omega_plasma = None #The plasma frequency
#
        self.ref_ewidth = None
        self.ref_Nw     = None
        self.ref_Nk     = None
        self.ref_Nb     = None
        self.FL = None

    @classmethod
    def generate(self, ED, Fortlib_option, Nene = 2000, ewidth = 0.001, plot_option = True, algorithm_option = 'ver1'):
#    def generate(self, ED, Fortlib_option, Nene = 2000, ewidth = 0.001, plot_option = True, algorithm_option = 'org'):
        self.Nene = Nene
        self.ewidth = ewidth
        emax = 1.2*(np.amax(ED.eigval) - np.amin(ED.eigval))
        emin = -1.0*emax
        self.omega = np.linspace(emin, emax, self.Nene)
        print('# emin, emax =',emin, emax)
        self.sigma = np.zeros([self.Nene, 3, 3], dtype='complex128')
        #
        print('Fortlib_option', Fortlib_option)
        if (not Fortlib_option):
            print('Fortlib_option', Fortlib_option)
            print('hoge')
            algorithm_option = algorithm_option + '_fort'
            self.Prep4Fortlib(self, ED)
        print('# ',algorithm_option, 'is chosen as algorithm_option.')
        if (algorithm_option == 'org'):
            self._compute_sigma_original(self, ED)
        elif (algorithm_option == 'ver1'):
            self._compute_sigma_ver1(self, ED)
        elif (algorithm_option == 'ver2'):
            self._compute_sigma_ver2(self, ED)
        if (algorithm_option == 'org_fort'):
            self._compute_sigma_original_Fortran(self, ED)
        else:
            print('# ERROR: No avilable option for algorithm_option with ', algorithm_option)
            sys.exit()
        #
        self.sigma = self.sigma*zI/ED.vcell/ED.Nk
        #Constructing epsilon from sigma
        self.epsilon = np.zeros([self.Nene, 3, 3], dtype='complex128')
        for ixyz in range(3):
            for jxyz in range(3):
                self.epsilon[:,ixyz,jxyz] = self.sigma[:,ixyz,jxyz]*fpi*zI/self.omega
        for ixyz in range(3):
            self.epsilon[:,ixyz,ixyz] = np.ones(self.Nene, dtype='complex128') + self.epsilon[:,ixyz,ixyz]
        self.epsilon_inv = 1.0/self.epsilon  
        print('# Number of energy grid: Nene =', self.Nene)
        print('# Energy width for sigma-epsilon: ewidth =', self.ewidth, '[a.u.] =', self.ewidth*Hartree, '[eV]')
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
        return self.omega, self.sigma, self.epsilon, self.epsilon_inv

    @classmethod
    def check_sum(self, ED, plot_option = True, SI_CGS = 'CGS'):
        omegamin = np.sqrt(np.amin(self.omega**2))
        Nzero = np.argmin((self.omega - omegamin)**2)
        omega = self.omega[Nzero:]
        epsilon = self.epsilon[Nzero:,:,:]
        epsilon_inv = self.epsilon_inv[Nzero:,:,:]
        domega = omega[1] - omega[0]
        self.sum_epsilon = np.zeros([3,3],dtype='float64')
        self.sum_epsilon_inv = np.zeros([3,3],dtype='float64')
        for ixyz in range(3):
            for jxyz in range(3):
                self.sum_epsilon[ixyz,jxyz] = np.sum(omega[:]*np.imag(epsilon[:,ixyz,jxyz]))*domega
                self.sum_epsilon_inv[ixyz,jxyz] = -np.sum(omega[:]*np.imag(epsilon_inv[:,ixyz,jxyz]))*domega
        if (SI_CGS == 'SI'):
            #SI value for the plasma freq. calculaiton
            epsilon0 = 8.8541878128e-12  #F/m, The vacuum permitivity
            me = 9.1093837015e-31        #kg, The electron mass
            ee = 1.602176634e-19         #C, The elementary charge
            hbar = 6.582119569e-16       #eV s, The Dirac constant
            h = 4.135667696e-15          #eV s, The Planck constatnt
            eledensity = np.sum(ED.occ)/ED.Nk/(ED.vcell*Atomvolume)*1.0e27 #/m^3
            self.omega_plasma = np.sqrt((eledensity*ee**2)/(me*epsilon0)) #/s
            self.omega_plasma = self.omega_plasma*hbar/Hartree #Hartree
        elif (SI_CGS == 'CGS'):
            #With the CGS, plasma_freq = sqrt(4.0*pi*Nele/Volume)
            eledensity = np.sum(ED.occ)/ED.Nk/ED.vcell #Atomic unit
            self.omega_plasma = np.sqrt(fpi*eledensity) #Atomic unit
        else :
            print('# ERROR: No avilable option for SI_CGS with ', SI_CGS)
            sys.exit()
        print('# ', SI_CGS, 'way to evaluate plasma frequency is chosen.')
        print('# The conductivity sum rule: ',self.sum_epsilon[0,0])
        print('# The longtudinal f-sum rule',self.sum_epsilon_inv[0,0])
        print('# The sum should be ',0.5*pi*self.omega_plasma**2)
        return self.sum_epsilon, self.sum_epsilon_inv, self.omega_plasma

    def _compute_sigma_original(self, ED):
    #Trivial implementation based on the formula without any trick
        ts = time.time()
        Nevery = int(ED.Nk/10.0)
        print('# Following is progress of GenerateSigmaEpsilon.generate function. ')
        for ik in range(ED.Nk):
            if (ik%Nevery ==0):
                te = time.time()
                print('# '+str(np.round(ik/ED.Nk*100.0, decimals = 2))+' % is done.')
                print('# ', te-ts, 'sec')
            for ib in range(ED.Nb):
                for jb in range(ED.Nb):
                    if (ib != jb):
                        ene_denominator = self.omega[:] - (ED.eigval[jb,ik] - ED.eigval[ib,ik]) + zI*self.ewidth
                        ene_denominator = ene_denominator*(ED.eigval[jb,ik] - ED.eigval[ib,ik])
                        docc = ED.occ[ib,ik] - ED.occ[jb,ik]
                        for ixyz in range(3):
                            for jxyz in range(3):
                                moment = ED.pmat[ixyz,ib,jb,ik]*ED.pmat[jxyz,jb,ib,ik]
                                self.sigma[:,ixyz,jxyz] = self.sigma[:,ixyz,jxyz] + docc*moment/ene_denominator[:]
#
    def _compute_sigma_ver1(self, ED):
    #Reduction of orbital sum by thinking of index exchange
        ts = time.time()
        Nevery = int(ED.Nk/10.0)
        print('# Following is progress of GenerateSigmaEpsilon.generate function. ')
        for ik in range(ED.Nk):
            if (ik%Nevery ==0):
                te = time.time()
                print('# '+str(np.round(ik/ED.Nk*100.0, decimals = 2))+' % is done.')
                print('# ', te-ts, 'sec')
            for ib in range(ED.Nb):
                for jb in range(ib):
                    weight1 = (ED.occ[ib,ik] - ED.occ[jb,ik]) \
                        /(self.omega[:] - (ED.eigval[jb,ik] - ED.eigval[ib,ik]) + zI*self.ewidth) \
                        /(ED.eigval[jb,ik] - ED.eigval[ib,ik])
                    weight2 = (ED.occ[jb,ik] - ED.occ[ib,ik]) \
                        /(self.omega[:] - (ED.eigval[ib,ik] - ED.eigval[jb,ik]) + zI*self.ewidth) \
                        /(ED.eigval[ib,ik] - ED.eigval[jb,ik])
                    for ixyz in range(3):
                        for jxyz in range(3):
                            moment1 = ED.pmat[ixyz,ib,jb,ik]*ED.pmat[jxyz,jb,ib,ik]
                            moment2 = ED.pmat[ixyz,jb,ib,ik]*ED.pmat[jxyz,ib,jb,ik]
                            self.sigma[:,ixyz,jxyz] = self.sigma[:,ixyz,jxyz] + moment1*weight1 + moment2*weight2
#
    def _compute_sigma_ver2(self, ED):
    #Using vector-operation as much as possible, but slow
        ts = time.time()
        Nevery = int(self.Nene/10.0)
        print('# Following is progress of GenerateSigmaEpsilon.generate function. ')
        deigvalbk = np.zeros([ED.Nb, ED.Nb, ED.Nk], dtype='float64')
        doccbk = np.zeros([ED.Nb, ED.Nb, ED.Nk], dtype='float64')
        weightbk0 = np.zeros([ED.Nb, ED.Nb, ED.Nk], dtype='complex128')
        for ik in range(ED.Nk):
            for ib in range(ED.Nb):
                for jb in range(ED.Nb):
                    deigvalbk[ib,jb,ik] = ED.eigval[ib,ik] - ED.eigval[jb,ik]
                    doccbk[ib,jb,ik] = ED.occ[ib,ik] - ED.occ[jb,ik]
                    if (ib != jb):
                        weightbk0[ib,jb,ik] = -doccbk[ib,jb,ik]/deigvalbk[ib,jb,ik]
        te = time.time()
        print('# preparation of deigvalbk arrays is done.')
        print('# ', te-ts, 'sec')
        for i in range(self.Nene):
            if (i%Nevery == 0):
                te = time.time()
                print('# '+str(np.round(i/self.Nene*100.0, decimals = 2))+' % is done.')
                print('# ', te-ts, 'sec')
            for ixyz in range(3):
                for jxyz in range(3):
                    A = ED.pmat[ixyz,:,:,:]*np.transpose(ED.pmat[jxyz,:,:,:], (1,0,2))*weightbk0[:,:,:]/(self.omega[i] + deigvalbk[:,:,:] + zI*self.ewidth)
                    self.sigma[i,ixyz,jxyz] = np.sum(A)

#
    def Prep4Fortlib(self, ED):
        self.ref_ewidth = ct.byref(ct.c_double(self.ewidth))
        self.ref_Nene   = ct.byref(ct.c_int32(self.Nene))
        self.ref_Nk     = ct.byref(ct.c_int32(ED.Nk))
        self.ref_Nb     = ct.byref(ct.c_int32(ED.Nb))
        dir_name = os.path.dirname(os.path.abspath(__file__)).strip('modules')
        print('# Fortlib.so: ',dir_name+"Fortlib.so")
        self.FL = np.ctypeslib.load_library(dir_name+"Fortlib.so",".")
        self.FL.compute_sigma_.argtypes = [
            np.ctypeslib.ndpointer(dtype='float64'),     #omega
            np.ctypeslib.ndpointer(dtype='float64'),     #eigval
            np.ctypeslib.ndpointer(dtype='float64'),     #occ
            np.ctypeslib.ndpointer(dtype='complex128'),  #pmat
            ct.POINTER(ct.c_double),                     #ewidth
            ct.POINTER(ct.c_int32),                      #Nene
            ct.POINTER(ct.c_int32),                      #Nk
            ct.POINTER(ct.c_int32),                      #Nb
            np.ctypeslib.ndpointer(dtype='complex128'),] #sigma
        self.FL.compute_sigma_.restype = ct.c_void_p
#
    def _compute_sigma_original_Fortran(self, ED):
    #Trivial implementation based on the formula without any trick, with Fortran implementation
        self.FL.compute_sigma_(self.omega, ED.eigval, ED.occ, ED.pmat, self.ref_ewidth, self.ref_Nene, self.ref_Nk, self.ref_Nb, self.sigma)

