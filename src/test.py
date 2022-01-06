#!/usr/bin/python
# coding: UTF-8
# This is created 2021/12/24 by Y. Shinohara
import sys
from modules.constants import *
from modules.parameters import ElkFiles, ElkData

ED = ElkData()
EFs = ElkFiles()
#
argv = sys.argv
argc = len(argv)
#First standard input is the directory path we want
if (argc == 1 or argc ==2):
    print('# Error: No direcotry path.')
    sys.exit()
elif (argc == 3):
    dir_path = argv[1]
    print('# The direcotry path is "'+dir_path+'".')
    elk_ver = int(argv[2])
    print('# Elk version is "'+str(elk_ver)+'".')

#ED._read_KPOINTS(EFs, dir_path)
#ED._read_LATTICE(EFs, dir_path)
#ED._read_INFO(EFs, dir_path)
#ED._read_EIGVAL(EFs, dir_path)
#ED.get_eigval(EFs, dir_path)

#ED._read_PMAT(EFs, dir_path, elk_ver)
ED.get_eigval_pmat(EFs, dir_path, elk_ver)

from modules.generate_DoS import GenerateDoS

GDoS = GenerateDoS()
omega, DoS, NoS, occDoS, occNoS = GenerateDoS.generate(ED)

from modules.change_occ import ChangeOccupation

ChangeOccupation.get_EVTECB(ED)
T = 5.0e-2
#ChangeOccupation.single_temperature(ED, T)
Np = 0.4121255861054197
Np = 2.0*Np
Tp = 1.0*T/10.0
Th = 1.0*T/10.0
ChangeOccupation.double_temperature(ED, Np, Tp, Th)

omega, DoS, NoS, occDoS, occNoS = GenerateDoS.generate(ED)

#
from modules.generate_sigma_epsilon import GenerateSigmaEpsilon

omega, sigma, epsilon, epsilon_inv = GenerateSigmaEpsilon.generate(ED, ewidth = 0.005)
#omega_org, sigma_org, epsilon_org, epsilon_inv_org = GenerateSigmaEpsilon.generate_org(ED, ewidth = 0.005)
sum_epsilon, sum_epsilon_inv, omega_plasma = GenerateSigmaEpsilon.check_sum(ED)

