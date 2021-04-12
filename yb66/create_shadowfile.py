import sys
import os
import numpy

from orangecontrib.xoppy.util.xoppy_util import locations, XoppyPhysics
from xoppy_xraylib_util2 import bragg_calc2
from crystal_shadow import crystal_shadow
from dabax_util import Crystal_GetCrystalsList

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Calculation structure factor')
    parser.add_argument('-n','--name',dest='descriptor', default=['YB66'],type=str, nargs=1, help='Crystal name')
    parser.add_argument('-m','--m', metavar='H K L', default=[0,0,6],type=int, nargs=1, help='Miller indic [H, K, L]')
    parser.add_argument('-e','--e', dest='EngRange', default=[2006,2194,0.5],type=float, nargs=3, help='[emin,emax,estep]')
    parser.add_argument('-s','--SHADOWFILE', dest='SHADOW_NAME', default=[""],type=str, nargs=3, help='SHADOW filename')
    
    args = parser.parse_args()
    descriptor = args.descriptor[0]
    HMILLER = args.m[0]
    KMILLER = args.m[1]
    LMILLER = args.m[2]
    ENERGY = args.EngRange[0]
    ENERGY_END = args.EngRange[1]
    estep = args.EngRange[2]
    NPOINTS = int((ENERGY_END-ENERGY)/estep + 1)
    SHADOW_NAME = args.SHADOW_NAME[0]
    energy = numpy.linspace(ENERGY,ENERGY_END,NPOINTS)

    print("Using crystal descriptor: ",descriptor)
    bragg_dictionary = bragg_calc2(descriptor=descriptor,hh=HMILLER,kk=KMILLER,ll=LMILLER,temper=1.0,
                                            emin=ENERGY,emax=ENERGY_END,estep=estep,fileout=None)   #50eV, replaced with estep

    if SHADOW_NAME=='':
        SHADOW_NAME = descriptor + '_' + str(HMILLER) + str(KMILLER) + str(LMILLER) + '_sha.dat'
    crystal_shadow(SHADOW_NAME,bragg_dictionary,energy)    
