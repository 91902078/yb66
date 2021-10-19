import numpy
from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_calc
from xoppy_dabax_util import bragg_calc2
from run_diff_pat import run_diff_pat
from srxraylib.plot.gol import plot
import os

if __name__ == "__main__":

    #
    # old code
    #

    os.system("rm -f xcrystal.bra xoppy.inp")
    dic1a = bragg_calc(descriptor="Si",hh=1,kk=1,ll=1,temper=1.0,
                       emin=5000.0,emax=12000.0,estep=5.0,fileout="xcrystal.bra")
    print("KEYS: ",dic1a.keys())
    print(dic1a)


    run_diff_pat(
        GEOMETRY = 0,
        SCAN = 2,
        UNIT = 1,
        SCANFROM = -100.0,
        SCANTO = 100.0,
        SCANPOINTS = 200,
        ENERGY = 10000.0,
        ASYMMETRY_ANGLE = 0.0,
        THICKNESS = 0.7,)

    a1 = numpy.loadtxt("diff_pat.dat",skiprows=5)

    #
    # New code
    #

    dic2a = bragg_calc2(descriptor="Si",hh=1,kk=1,ll=1,temper=1.0,
                       emin=5000.0,emax=12000.0,estep=5.0,fileout="xcrystal.bra")
    print("KEYS: ",dic2a.keys())
    print(dic2a)




    run_diff_pat(
        GEOMETRY = 0,
        SCAN = 2,
        UNIT = 1,
        SCANFROM = -100.0,
        SCANTO = 100.0,
        SCANPOINTS = 200,
        ENERGY = 10000.0,
        ASYMMETRY_ANGLE = 0.0,
        THICKNESS = 0.7,)

    a2 = numpy.loadtxt("diff_pat.dat",skiprows=5)

    #
    # comparison
    #

    print("\n\n\n\n\n\n")
    for key in dic1a.keys():
        print(key,dic1a[key],dic2a[key])


    plot(a1[:, 0], a1[:, -1],
         a2[:, 0], a2[:, -1],
         linestyle=[None,''],
         marker=[None,"+"],
         legend=["XOPPY bragg_calc code","NEW bragg_calc2 code"])



