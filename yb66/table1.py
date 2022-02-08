import numpy
import os

from srxraylib.plot.gol import plot, set_qt
# from xoppylib.decorators.xraylib_decorated import XraylibDecorated
from xoppylib.decorators.dabax_decorated import DabaxDecorated
import scipy.constants as codata

from crystalpy.diffraction.DiffractionSetupDabax import DiffractionSetupDabax


def run(descriptor, hkl, material_constants_library):
    import os
    import xraylib

    from xoppylib.crystals.tools import bragg_calc, bragg_calc2
    from xoppylib.crystals.tools import run_diff_pat



    os.system("rm -f xcrystal.bra xoppy.inp")

    dic1a = bragg_calc2(descriptor=descriptor,hh=int(hkl[0]),kk=int(hkl[1]),ll=int(hkl[2]),temper=1.0,
                       emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra",
                       material_constants_library=material_constants_library)

    print("KEYS: ",dic1a.keys())


    run_diff_pat(
        dic1a,
        preprocessor_file="xcrystal.bra",
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

    aa = numpy.loadtxt("diff_pat.dat",skiprows=5)
    return aa, dic1a



if __name__ == "__main__":

    set_qt()

    # xrl = XraylibDecorated()
    dx0 = DabaxDecorated(file_Crystals="Crystals.dat")
    dx1 = DabaxDecorated(file_Crystals="Crystals_xrayserver.dat")

    do_plot = 1


#
# table 3
#

    list0 = dx0.Crystal_GetCrystalsList()
    list1 = dx1.Crystal_GetCrystalsList()

    mylist =        ["KCl", "Diamond", "InSb", "ADP", "AlphaQuartz", "PET", "Beryl", "Muscovite", "TlAP", "RbAP", "KAP"]
    myreflections = ["200", "002",     "111",  "200", "110",         "002", "110",   "002",       "001",  "001",  "001"]
    darwin_ener   = [2.0,    2.0,       1.7,   1.73,   1.5,          1.48,  0.81,    0.65,        0.5,    0.5,    0.5  ]
    dspacing=[]
    emin = []
    darwin_width = []
    peak_reflectivity = []
    file_found = []

    i = 0
    i0found = 0
    i1found = 0
    a0 = None
    a1 = None

    n = len(mylist)

    for i in range(n):
        if mylist[i] in list0:
            print(">>> %s found in Crystals.dat" % (mylist[i]))
            c0 = dx0.Crystal_GetCrystal(mylist[i])
            ds = dx0.Crystal_dSpacing(c0,int(myreflections[i][0]),int(myreflections[i][1]),int(myreflections[i][2]))
            ee = codata.h * codata.c / codata.e / (2 * ds *1e-10)
            dspacing.append(ds)
            emin.append(ee)
            file_found.append(0)
            # a0,dic1a = run(mylist[i], myreflections[i], dx0)

            cpy = DiffractionSetupDabax(geometry_type=0, crystal_name=mylist[i], thickness=1e-5,
                 miller_h=int(myreflections[i][0]), miller_k=int(myreflections[i][1]), miller_l=int(myreflections[i][2]),
                 asymmetry_angle=0.0,
                 azimuthal_angle=0.0,)

            darwin_width.append( cpy.darwinHalfwidthS(1e3*darwin_ener[i]) )

        elif mylist[i] in list1:
            print(">>> %s found in Crystals_xrayserver.dat" % (mylist[i]))
            c0 = dx1.Crystal_GetCrystal(mylist[i])
            ds = dx1.Crystal_dSpacing(c0,int(myreflections[i][0]),int(myreflections[i][1]),int(myreflections[i][2]))
            ee = codata.h * codata.c / codata.e / (2 * ds *1e-10)
            dspacing.append(ds)
            emin.append(ee)
            file_found.append(1)
            # a0,dic1a = run(mylist[i], myreflections[i], dx0)
            darwin_width.append(0.0)
        else:
            print(">>> %s NOT found" % (mylist[i]))
            dspacing.append(-1)
            emin.append(0.0)
            file_found.append(-1)
            darwin_width.append(0.0)


                # if mylist[i] in list0:
    #     print(">>> %s found in Crystals_xrayserver.dat" % (mylist[i]))
    #     a1,dic1a = run(mylist[i], myreflections[i], dx1)

    #
    #
    # plot(a0[:, 0], a0[:, -1],
    #      a1[:, 0], a1[:, -1],
    #      linestyle=[None,'--'],
    #      marker=[None,"+"],
    #      title=mylist[i],
    #      legend=["using DABAX Crystals.dat","using DABAX Crystals_xrayserver.dat"], show=1,
    #      xtitle=r"$\theta-\theta_B$ [$\mu$rad]", ytitle="reflectivity")
    #
    #
    # print(a0.shape)

    f = open("table1.txt", 'w')
    f.write("         NAME (hkl)     dspacing   found in  emin  DarwinWidth\n")
    for i in range(len(mylist)):
        f.write("%15s (%s)   %10.5g   %3d  %g  %g\n" % (
            mylist[i],
            myreflections[i],
            2 * dspacing[i],
            file_found[i],
            emin[i],
            2e6*darwin_width[i],
        ))




