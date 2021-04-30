import sys
import os
import numpy
from crystal_util import bragg_calc2, crystal_fh2
from run_diff_pat import run_diff_pat
from srxraylib.plot.gol import plot
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Calculation structure factor')
    parser.add_argument('-n','--name',dest='descriptor', default=['YB66'],type=str, nargs=1, help='Crystal name')
    parser.add_argument('-m','--m', metavar='H K L', default=[4,0,0],type=int, nargs=3, help='Miller indic [H, K, L]')
    parser.add_argument('-e','--e', dest='Energy', default=[8040],type=float, nargs=1, help='Energy')
    parser.add_argument('-s','--s', dest='ScanRange', default=[-100,100,200],type=float, nargs=3, help='[SCANFROM,SCANTO,SCANPOINTS]')
    parser.add_argument('-t','--t', dest='Thickness', default=[0.7],type=float, nargs=1, help='Crystal thicknes [cm]')
    parser.add_argument('-a','--a', dest='Anisotropy', default=[0],type=int, nargs=1, help='Anisotropy (0-3)')
    parser.add_argument('-g','--g', dest='GEOMETRY', default=[0],type=int, nargs=1, help='0 BRAGG: diffr beam, 1 LAUE: diffr beam, 2 BRAGG: transm beam, 3 LAUE: transm beam')
    parser.add_argument('--MOSAIC', dest='MOSAIC', default=[0],type=int, nargs=1, help='MOSAIC (0-1)')
    parser.add_argument('--MOSAIC_FWHM', dest='MOSAIC_FWHM', default=[0.1],type=float, nargs=1, help='MOSAIC_FWHM (Deg)')
    parser.add_argument('--ASYMMETRY_ANGLE', dest='ASYMMETRY_ANGLE', default=[0.0],type=float, nargs=1, help='ASYMMETRY_ANGLE')
    parser.add_argument('--SCAN', dest='SCAN', default=[1],type=int, nargs=1, help='0 Theta (absolute), 1 Th - Th Bragg (corrected), 2 Th - Th Bragg, Energy [eV],3 y (Zachariasen)')
    parser.add_argument('--RSAG', dest='RSAG', default=[125.0],type=float, nargs=1, help='Scan: 0 Th absolute 1 Th-Th Bragg (cor) 2 Th-Th Bragg')
    parser.add_argument('--RMER', dest='RMER', default=[1290.0],type=float, nargs=1, help='Scan: 0 Th absolute 1 Th-Th Bragg (cor) 2 Th-Th Bragg')
    parser.add_argument('--POISSON', dest='POISSON', default=[0.22],type=float, nargs=1, help='Scan: 0 Th absolute 1 Th-Th Bragg (cor) 2 Th-Th Bragg')
    parser.add_argument('--CUT',dest='CUT', default=["2 -1 -1;  1 1 1 ; 0 0 0"],type=str, nargs=1, help='Cut')
    parser.add_argument('--FILECOMPLIANCE',dest='FILECOMPLIANCE', default=["mycompliance.dat"],type=str, nargs=1, help='FILECOMPLIANCE')

    args = parser.parse_args()
    descriptor = args.descriptor[0]
    UNIT = 1
    TEMPER = 1.0
    THICKNESS = args.Thickness[0]
    MILLER_INDEX_H = args.m[0]
    MILLER_INDEX_K = args.m[1]
    MILLER_INDEX_L = args.m[2]
    SCANFROM = args.ScanRange[0]
    SCANTO = args.ScanRange[1]
    SCANPOINTS = int(args.ScanRange[2])
    estep = (SCANTO-SCANFROM)/SCANPOINTS
    ENERGY = args.Energy[0]
    ANISOTROPY = args.Anisotropy[0]
    GEOMETRY = args.GEOMETRY[0]
    MOSAIC = args.MOSAIC[0]
    MOSAIC_FWHM = args.MOSAIC_FWHM[0]
    ASYMMETRY_ANGLE = args.ASYMMETRY_ANGLE[0]
    SCAN = args.SCAN[0]
    RSAG = args.RSAG[0]
    RMER = args.RMER[0]
    POISSON = args.POISSON[0]
    CUT = args.CUT[0]
    FILECOMPLIANCE = args.FILECOMPLIANCE[0]



    if SCAN == 3: # energy scan
        emin = SCANFROM - 1
        emax = SCANTO + 1
    else:
        emin = ENERGY - 100.0
        emax = ENERGY + 100.0


    print("Using crystal descriptor: ",descriptor)

    bragg_dictionary = bragg_calc2(descriptor=descriptor,
                                            hh=MILLER_INDEX_H,kk=MILLER_INDEX_K,ll=MILLER_INDEX_L,
                                            temper=TEMPER,
                                            emin=emin,emax=emax,estep=estep,fileout="xcrystal.bra")

    run_diff_pat(
        MOSAIC = 0,
        GEOMETRY = 0,
        SCAN = 2,
        UNIT = 1,
        SCANFROM = -100.0,
        SCANTO = 100.0,
        SCANPOINTS = 200,
        ENERGY = 8040.0,
        ASYMMETRY_ANGLE = 0.0,
        THICKNESS = 0.7,
        MOSAIC_FWHM = 0.1,
        RSAG = 125.0,
        RMER = 1290.0,
        ANISOTROPY = 0,
        POISSON = 0.22,
        CUT = "2 -1 -1 ; 1 1 1 ; 0 0 0",
        FILECOMPLIANCE = "mycompliance.dat")

    a = numpy.loadtxt("diff_pat.dat",skiprows=5)

    #
    # plot
    #

    plot(a[:, 0], a[:, -1])

