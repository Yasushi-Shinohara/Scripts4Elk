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
ED.get_eigval(EFs, dir_path)
#ED._read_PMAT(EFs, dir_path, elk_ver)
ED.get_eigval_pmat(EFs, dir_path, elk_ver)