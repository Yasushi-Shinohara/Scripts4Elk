#!/usr/bin/python
# coding: UTF-8
# This is created 2022/01/05 by Y. Shinohara
import sys
from modules.constants import *
from modules.parameters import ElkFiles, ElkData

ED = ElkData()
EFs = ElkFiles()
#
argv = sys.argv
argc = len(argv)
#First standard input is the directory path we want
if (argc == 1 or argc ==2 or argc ==3):
    print('# Error: No direcotry path.')
    sys.exit()
elif (argc == 4):
    dir_path = argv[1]
    print('# The direcotry path is "'+dir_path+'".')
    elk_ver = int(argv[2])
    print('# Elk version is "'+str(elk_ver)+'".')
    Fortlib_option = argv[3]
    print('# Fortlib_option is "'+str(Fortlib_option)+'".')

ED.get_eigval_pmat(EFs, dir_path, elk_ver)

from modules.generate_sigma_epsilon import GenerateSigmaEpsilon
omega, sigma, epsilon, epsilon_inv = GenerateSigmaEpsilon.generate(ED, Fortlib_option, ewidth = 0.002)
sum_epsilon, sum_epsilon_inv, omega_plasma = GenerateSigmaEpsilon.check_sum(ED)

