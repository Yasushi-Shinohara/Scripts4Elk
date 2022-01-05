# coding: UTF-8
# This is created 2022/01/05 by Y. Shinohara
# This is lastly modified 2021/xx/yy by Y. Shinohara
import sys
from modules.constants import pi, Hartree
import numpy as np

class ChangeOccupation:
    def __init__(self):
        self.system_type = None #'ST' or 'DT': single- or double-temperature
        self.nVT = None
        self.nCB = None
        self.EVT = None         #Energy of valence top
        self.ECB = None         #Energy of conduction bottom
        #
        self.T = None           #Temperature for whole system within 'ST'
        self.beta = None
        self.mu
        #
        self.Tp = None          #Temperature for "particle" in conduction band
        self.Th = None          #Temperature for hole in valence band
        self.betap = None
        self.betah = None
        self.Np = None          #Number of particle per unit cell
        self.Nh = None          #Number of hole per unit cell
        self.mup = None         #Chemical potential for particle
        self.muh = None         #Chemical potential for hole

    @classmethod
    def get_EVTECB(self, ED):
        self.nVT = int(ED.Ne/ED.spin_degeneracy)
        self.nCB = int(ED.Ne/ED.spin_degeneracy) + 1
        self.EVT = np.amax(ED.eigval[:self.nVT,:])
        self.ECB = np.amin(ED.eigval[self.nCB-1:,:])
        print('# EVT =', self.EVT)
        print('# ECB =', self.ECB)
        
    @classmethod
    def single_temperature(self, ED, T, plot_option = True, Ne_epsilon = 1.0e-8, Nitermax = 200):
        exp_maxvalue = 1.0e2 #Maximum value for exponential argument
        self.T = 1.0*T
        self.beta = 1.0/self.T
        mumin = np.amin(ED.eigval)
        mumax = np.amax(ED.eigval)
        #
        self.mu = 0.5*(mumin + mumax)
        print('# Initial chemical potential: ',self.mu)
        temp = self.beta*(ED.eigval[:,:] - self.mu)
        temp = np.clip(temp, None, exp_maxvalue)
        occ = np.exp(temp) + 1.0
        occ = np.float(ED.spin_degeneracy)/occ
        Ne = np.sum(occ)/ED.Nk
        print('# Initial number of electron',Ne)
        print('# Iteration to find proper chemical potentials starts...')
        print('# ncounter, chemical potential, Ne, expected Ne ')
        ncounter = 0
        while (np.abs(Ne - ED.Ne) > Ne_epsilon):
            ncounter += 1
            if (Ne < ED.Ne):
                mumin = self.mu  
            else:
                mumax = self.mu  
            self.mu = 0.5*(mumin + mumax)
            temp = np.clip(self.beta*(ED.eigval[:,:] - self.mu), None, exp_maxvalue)
            occ = np.float(ED.spin_degeneracy)/(np.exp(temp) + 1.0)
            Ne = np.sum(occ)/ED.Nk
            print('# ',ncounter, self.mu, Ne, ED.Ne)
            if (ncounter == Nitermax):
                print('# The iteration is NOT converged in ',ncounter,' iterations.')
                break
        else:
            print('# The iteration is converged properly.')
            print('# ====================================')
        Nv = np.sum(occ[:self.nVT,:])/ED.Nk
        Nc = np.sum(occ[self.nCB-1:,:])/ED.Nk
        print('# The fundamental gap', (self.ECB - self.EVT), '[a.u.] =',(self.ECB - self.EVT)*Hartree, '[eV] .')
        print('# The temperature', self.T, '[a.u.] =',self.T*Hartree, '[eV] .')
        print('# The chimical potential', self.mu, '[a.u.] =',self.mu*Hartree, '[eV] .')
        print('# Number of thermally excited electron in the conduction band:', Nc)
        ED.occ = occ

#
    @classmethod
    def double_temperature(self, ED, Np_in, Tp, Th, plot_option = True, Ne_epsilon = 1.0e-8, Nitermax = 200):
        exp_maxvalue = 1.0e2 #Maximum value for exponential argument
        self.Np = 1*Np_in
        self.Nh = self.Np
        self.Tp = 1.0*Tp
        self.betap = 1.0/self.Tp
        self.Th = 1.0*Th
        self.betah = 1.0/self.Th
        #Particle
        mumin = np.amin(ED.eigval)
        mumax = np.amax(ED.eigval)
        self.mup = 0.5*(mumin + mumax)
        print('# Initial chemical potential for particle: ',self.mup)
        temp = self.betap*(ED.eigval[:,:] - self.mup)
        temp = np.clip(temp, None, exp_maxvalue)
        occp = np.exp(temp) + 1.0
        occp = np.float(ED.spin_degeneracy)/occp
        Np = np.sum(occp[self.nCB-1:,:])/ED.Nk
        print('# Initial number of particle',Np)
        print('# Iteration to find proper chemical potentials starts...')
        print('# ncounter, chemical potential, Np, expected Np ')
        ncounter = 0
        while (np.abs(Np - self.Np) > Ne_epsilon):
            ncounter += 1
            if (Np < self.Np):
                mumin = self.mup
            else:
                mumax = self.mup 
            self.mup = 0.5*(mumin + mumax)
            temp = np.clip(self.betap*(ED.eigval[:,:] - self.mup), None, exp_maxvalue)
            occp = np.float(ED.spin_degeneracy)/(np.exp(temp) + 1.0)
            Np = np.sum(occp[self.nCB-1:,:])/ED.Nk
            print('# ',ncounter, self.mup, Np, self.Np)
            if (ncounter == Nitermax):
                print('# The iteration is NOT converged in ',ncounter,' iterations.')
                break
        else:
            print('# The iteration is converged properly.')
            print('# ====================================')
        Np = np.sum(occp[self.nCB-1:,:])/ED.Nk
        #Hole
        mumin = np.amin(ED.eigval)
        mumax = np.amax(ED.eigval)
        self.muh = 0.5*(mumin + mumax)
        print('# Initial chemical potential for particle: ',self.muh)
        temp = self.betap*(ED.eigval[:,:] - self.muh)
        temp = np.clip(temp, None, exp_maxvalue)
        occh = np.exp(temp) + 1.0
        occh = np.float(ED.spin_degeneracy)/occh
        Nh = ED.Ne - np.sum(occh[:self.nVT,:])/ED.Nk
        print('# Initial number of hole',Nh)
        print('# Iteration to find proper chemical potentials starts...')
        print('# ncounter, chemical potential, Nh, expected Nh ')
        ncounter = 0
        while (np.abs(Nh - self.Nh) > Ne_epsilon):
            ncounter += 1
            if (Nh > self.Nh):
                mumin = self.muh
            else:
                mumax = self.muh 
            self.muh = 0.5*(mumin + mumax)
            temp = np.clip(self.betah*(ED.eigval[:,:] - self.muh), None, exp_maxvalue)
            occh = np.float(ED.spin_degeneracy)/(np.exp(temp) + 1.0)
            Nh = ED.Ne - np.sum(occh[:self.nVT,:])/ED.Nk
            print('# ',ncounter, self.muh, Nh, self.Nh)
            if (ncounter == Nitermax):
                print('# The iteration is NOT converged in ',ncounter,' iterations.')
                break
        else:
            print('# The iteration is converged properly.')
            print('# ====================================')
        Nh = ED.Ne - np.sum(occh[:self.nVT,:])/ED.Nk
        print('# The fundamental gap', (self.ECB - self.EVT), '[a.u.] =',(self.ECB - self.EVT)*Hartree, '[eV] .')
        print('# The particle temperature', self.Tp, '[a.u.] =',self.Tp*Hartree, '[eV] .')
        print('# The particle chimical potential', self.mup, '[a.u.] =',self.mup*Hartree, '[eV] .')
        print('# Number of particle :', Np)
        print('# The hole temperature', self.Th, '[a.u.] =',self.Th*Hartree, '[eV] .')
        print('# The hole chimical potential', self.muh, '[a.u.] =',self.muh*Hartree, '[eV] .')
        print('# Number of hole :', Nh)
        ED.occ[self.nCB-1:,:] = occp[self.nCB-1:,:]
        ED.occ[:self.nVT,:] = occh[:self.nVT,:]

    
