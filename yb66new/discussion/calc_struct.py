import sys
import numpy
from xoppy_xraylib_util2 import crystal_fh2, bragg_calc2
import math



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Calculation structure factor')
    parser.add_argument('-n','--name',dest='descriptor', default=['YB66'],type=str, nargs=1, help='Crystal name')
    parser.add_argument('-m','--m', metavar='H K L', default=[4,0,0],type=int, nargs='+', help='Miller indic [H, K, L]')
    parser.add_argument('-e','--e', dest='EngRange', metavar='emin,emax,estep', default=[8040,8050,1],type=float, nargs='+', help='[emin,emax,estep]')
    
    args = parser.parse_args()
    descriptor = args.descriptor[0].strip()
    HMILLER = args.m[0]
    KMILLER = args.m[1]
    LMILLER = args.m[2]
    ENERGY = args.EngRange[0]
    ENERGY_END = args.EngRange[1]
    estep = args.EngRange[2]
    NPOINTS = int((ENERGY_END-ENERGY)/estep + 1)
    
    print("Using crystal descriptor: ",descriptor)
    bragg_dictionary = bragg_calc2(descriptor=descriptor,hh=HMILLER,kk=KMILLER,ll=LMILLER,temper=1.0,
                                            emin=ENERGY,emax=ENERGY_END,estep=estep,fileout=None)   #50eV, replaced with estep
    energy = numpy.linspace(ENERGY,ENERGY_END,NPOINTS)

    for i,ienergy in enumerate(energy):
        dic2 = crystal_fh2(bragg_dictionary,ienergy)
        print("Energy=%g eV FH=(%g,%g)"%(ienergy,dic2["STRUCT"].real,dic2["STRUCT"].imag))

