# coding: UTF-8
# This is created 2021/12/24 by Y. Shinohara
# This is lastly modified 2021/xx/yy by Y. Shinohara
import sys
from modules.constants import tpi, Atomtime, Hartree, Atomfield, aB
import numpy as np

class ElkInfo:
    def __init__(self, elk_ver):
        self.version = elk_ver

#
class ElkFiles:
    def __init__(self):
        self.elkin = 'elk.in'
        self.BANDLINES = 'BANDLINES.OUT'
        self.BAND = 'BAND.OUT'
        self.DTOTENERGY = 'DTOTENERGY.OUT'
        self.EFERMI = 'EFERMI.OUT'
        self.EIGVAL = 'EIGVAL.OUT'
        self.EPSILON_11 = 'EPSILON_11.OUT'
        self.FERMIDOS = 'FERMIDOS.OUT'
        self.GAP = 'GAP.OUT'
        self.GEOMETRY = 'GEOMETRY.OUT'
        self.IDOS = 'IDOS.OUT'
        self.INFO = 'INFO.OUT'
        self.KPOINTS = 'KPOINTS.OUT'
        self.LATTICE = 'LATTICE.OUT'
        self.PMAT = 'PMAT.OUT'
        self.SIGMA_11 = 'SIGMA.OUT'
        self.STATE = 'STATE.OUT'
        self.TDOS = 'TDOS.OUT'
        self.TOTENERGY = 'TOTENERGY.OUT'
#
class ElkData:
    def __init__(self):
        self.Nk = None
        self.vkl = None
        self.a1 = None
        self.a2 = None
        self.a3 = None
        self.b1 = None
        self.b2 = None
        self.b3 = None
        self.vcell = None
        self.Nb = None
        self.eigval = None
        self.vkc = None
        self.pmat = None
        self.occ = None
        self.Ne = None
        self.spin_treatment = None
        self.spin_degeneracy = None
#
    def _read_KPOINTS(self, EFs, dir_path):
        file_KPOINTS = dir_path + '/' + EFs.KPOINTS
        f = open(file_KPOINTS,'r')
        lines = f.readlines()
        f.close()
        #
        temp = lines[0].split()
        self.Nk = int(temp[0])
        print ('Nk =',self.Nk)
        self.vkl = np.zeros([3,self.Nk])
        for ik in range(1,len(lines)):
            temp = lines[ik].split()
            for i in range(3):
                self.vkl[i,ik-1] = float(temp[i+1])

    def _read_LATTICE(self, EFs, dir_path):
        file_LATTICE = dir_path + '/' + EFs.LATTICE
        f = open(file_LATTICE,'r')
        lines = f.readlines()
        f.close()
        for i in range(len(lines)):
            temp = lines[i].split(':')
            if (str(temp[0]) == 'vector a1 ') :
                self.a1 = temp[1].split()
                print ('a1 = ',self.a1)
            elif (str(temp[0]) == 'vector a2 ') :
                self.a2 = temp[1].split()
                print ('a2 = ',self.a2)
            elif (str(temp[0]) == 'vector a3 ') :
                self.a3 = temp[1].split()
                print ('a3 = ',self.a3)
            elif (str(temp[0]) == 'vector b1 ') :
                self.b1 = temp[1].split()
                print ('b1 = ',self.b1)
            elif (str(temp[0]) == 'vector b2 ') :
                self.b2 = temp[1].split()
                print ('b2 = ',self.b2)
            elif (str(temp[0]) == 'vector b3 ') :
                self.b3 = temp[1].split()
                print ('b3 = ',self.b3)
#
        temp1 = self.a1
        temp2 = self.a2
        temp3 = self.a3
        self.a1 = np.zeros(3)
        self.a2 = np.zeros(3)
        self.a3 = np.zeros(3)
        for ixyz in range(3):
            self.a1[ixyz] = float(temp1[ixyz])
            self.a2[ixyz] = float(temp2[ixyz])
            self.a3[ixyz] = float(temp3[ixyz])
        self.vcell = np.dot(np.cross(self.a1,self.a2),self.a3)
        print ('volume of cell in real-space, vcell = ',self.vcell)
        temp1 = self.b1
        temp2 = self.b2
        temp3 = self.b3
        self.b1 = np.zeros(3)
        self.b2 = np.zeros(3)
        self.b3 = np.zeros(3)
        for ixyz in range(3):
            self.b1[ixyz] = float(temp1[ixyz])
            self.b2[ixyz] = float(temp2[ixyz])
            self.b3[ixyz] = float(temp3[ixyz])
#
    def _read_INFO(self, EFs, dir_path):
        file_INFO = dir_path + '/' + EFs.INFO
        f = open(file_INFO,'r')
        lines = f.readlines()
        f.close()
        for i in range(len(lines)):
            temp = lines[i].split(':')
            if (str(temp[0]) == 'k-point grid ') :
                ngridk = temp[1].split()
                print ('# ngridk = ',ngridk)
            elif (str(temp[0]) == 'k-point offset ') :
                self.vkloff = temp[1].split()
                print ('# vkloff = ',self.vkloff)
            elif (str(temp[0]) == 'Total valence charge    ') :
                self.Ne = float(temp[1])
                print ('# Ne = ',self.Ne)
            elif (str(temp[0]) == 'Spin treatment ') :
                self.spin_treatment = str(lines[i+1]).strip()
                if (self.spin_treatment.strip() == 'spin-unpolarised'):
                    self.spin_degeneracy = 2
                print('# Spin treatment: ', self.spin_treatment)
                print('# Spin degeneracy: ', self.spin_degeneracy)
        Nk1 = int(ngridk[0])
        Nk2 = int(ngridk[1])
        Nk3 = int(ngridk[2])
        Nkt = Nk1*Nk2*Nk3
        if (self.Nk != Nkt):
            print ('# Error stop: Inconsistency between Nkt and Nk')
            sys.exit()
        
        self.vkc = np.zeros([3,self.Nk])
        for ik in range(self.Nk):
            self.vkc[:,ik] = self.b1*self.vkl[0,ik] + self.b2*self.vkl[1,ik] + self.b3*self.vkl[2,ik]
        kx = self.vkc[0,:]
        ky = self.vkc[1,:]
        kz = self.vkc[2,:]
#
    def _read_EIGVAL(self, EFs, dir_path):
        file_EIGVAL = dir_path + '/' + EFs.EIGVAL
        f = open(file_EIGVAL,'r')
        lines = f.readlines()
        f.close()
        temp = lines[0].split(':')
        Nkt = int(temp[0])
        if (self.Nk != Nkt):
            print ('# Error stop: Inconsistency between Nkt and Nk')
            sys.exit()
        temp = lines[1].split(':')
        self.Nb = int(temp[0])
        print ('Nb =', self.Nb)
        self.eigval = np.zeros([self.Nb,self.Nk])
        for ik in range(self.Nk):
            for ib in range(self.Nb):
                i = ik*(self.Nb+4) + ib+5
                temp = lines[i].split()
                self.eigval[ib,ik] = float(temp[1])
#
    def get_eigval(self, EFs, dir_path):
        self._read_KPOINTS(EFs, dir_path)
        self._read_LATTICE(EFs, dir_path)
        self._read_INFO(EFs, dir_path)
        self._read_EIGVAL(EFs, dir_path)
#
    def _read_PMAT(self, EFs, dir_path, elk_ver):
        file_PMAT = dir_path + '/' + EFs.PMAT
#Sector for read a binary file                                              
#See: http://ig.hateblo.jp/entry/2014/05/30/225607                          
        head = ("head","<i") #In the elk code, head is not substituted due to use access="direct" option
        tail = ("tail","<i") #In the elk code, tail is not substituted due to use access="direct" option
        chartemp = '<'+str(self.Nb**2*3)+'c16'
        datatype = np.dtype([("vec","<3d"), ("N","<i"), ("array",chartemp)]) #c16 means 128-bit complex
#datatype = np.dtype([("vec","<3d"), ("N","<i"), ("array","<1323c16")]) #c16 means 128-bit complex
# "<3d" , "<i" mean three double precision and integer                      
        f = open(file_PMAT,'r')
        chunk = np.fromfile(f, dtype=datatype, count=self.Nk) #NK times read of file "f" with given dtype
        f.close()
        vkl_d = np.transpose(chunk[:]["vec"])
        NB_d =  np.transpose(chunk[:]["N"])
        self.pmat = np.zeros([3,self.Nb,self.Nb,self.Nk],dtype=np.complex128) #np.comple128 means 128-bit complex
#In the elk, p matrix-element is one for \psi_k = e^ikr * u_k rather than lattice periodic part u_k.
        for ik in range(self.Nk):
            self.pmat[:,:,:,ik] = chunk[ik]["array"].reshape((3,self.Nb,self.Nb),order='F')
#Order in the array, pmat, is different among elk-2.x.y and elk-3.x.y.      
        if (elk_ver >= 3):
            for ik in range(self.Nk):
                temp = chunk[ik]["array"].reshape((self.Nb,self.Nb,3),order='F')
                self.pmat[0,:,:,ik] = temp[:,:,0]
                self.pmat[1,:,:,ik] = temp[:,:,1]
                self.pmat[2,:,:,ik] = temp[:,:,2]
#
    def get_eigval_pmat(self, EFs, dir_path,elk_ver):
        self._read_KPOINTS(EFs, dir_path)
        self._read_LATTICE(EFs, dir_path)
        self._read_INFO(EFs, dir_path)
        self._read_EIGVAL(EFs, dir_path)
        self._read_PMAT(EFs, dir_path,elk_ver)
