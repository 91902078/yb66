import numpy
from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_calc
from xoppy_dabax_util import bragg_calc2
from run_diff_pat import run_diff_pat
from srxraylib.plot.gol import plot
import os

if __name__ == "__main__":
    descriptor = 'Muscovite'

    #
    # old code
    #

    os.system("rm -f xcrystal.bra xoppy.inp")
    dic1a = bragg_calc(descriptor=descriptor,hh=1,kk=1,ll=1,temper=1.0,emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra")

    print("KEYS: ",dic1a.keys())
    print(dic1a)
    os.system("cp xcrystal.bra xcrystal.bra.old")

    run_diff_pat(
        MOSAIC = 0,
        GEOMETRY = 0,
        SCAN = 2,
        UNIT = 1,
        SCANFROM = 25.0,
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

    a1 = numpy.loadtxt("diff_pat.dat",skiprows=5)


    #
    # New code
    #

    dic2a = bragg_calc2(descriptor=descriptor,hh=1,kk=1,ll=1,temper=1.0,emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra")
    print("KEYS: ",dic2a.keys())
    print(dic2a)
    os.system("cp xcrystal.bra xcrystal.bra.new")


    # compare
    try:
        for key in dic1a.keys():
            print(key,dic1a[key],dic2a[key])
    except:
        pass


    run_diff_pat(
        MOSAIC = 0,
        GEOMETRY = 0,
        SCAN = 2,
        UNIT = 1,
        SCANFROM = 25.0,
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

    a2 = numpy.loadtxt("diff_pat.dat",skiprows=5)

    #
    # comparison
    #

    plot(a1[:, 0], a1[:, -1],
         a2[:, 0], a2[:, -1],
         linestyle=[None,''],
         marker=[None,"+"],
         legend=["XOPPY bragg_calc code","NEW bragg_calc2 code"])



