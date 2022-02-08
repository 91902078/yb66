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
    # dx0 = DabaxDecorated(file_Crystals="Crystals.dat")
    dx1 = DabaxDecorated(file_Crystals="Crystals_xrayserver.dat")

    print(dx1.info())

    do_plot = 1


#
# table 3
#

    # list0 = dx0.Crystal_GetCrystalsList()
    list1 = dx1.Crystal_GetCrystalsList()

    print(list1)

    mylist = list1
    mylib  = dx1


    n = len(mylist)
    darwin_ener   = [0.0] * n
    dspacing=[0.0] * n
    emin = [0.0] * n
    darwin_width = [0.0] * n
    peak_reflectivity = [0.0] * n
    H = [0] * n
    K = [0] * n
    L = [0] * n




    for i in range(n):
        c0 = mylib.Crystal_GetCrystal(mylist[i])
        for h in [2,1,0]:
            for k in [1,0]:
                for l in [1,0]:
                    ds = mylib.Crystal_dSpacing(c0, h, k, l)
                    ee = codata.h * codata.c / codata.e / (2 * ds * 1e-10)

                    cpy = DiffractionSetupDabax(geometry_type=0, crystal_name=mylist[i], thickness=1e-5,
                         miller_h=h,miller_k=k, miller_l=l,
                         asymmetry_angle=0.0,
                         azimuthal_angle=0.0,
                         dabax=mylib)

                    FH = cpy.FH(ee * 1.1)

                    # print(mylist[i], h, k, l, ds, numpy.abs(FH) )

                    if ( (h+k+l) > 0):
                        if numpy.abs(FH) > 1:
                            print(">>>>>> FOUND!", h,k,l)
                            H[i] = h
                            K[i] = k
                            L[i] = l
                            emin[i] = ee
                            dspacing[i] = ds
                            darwin_width[i] = cpy.darwinHalfwidthS(ee * 1.1)


    f = open("table1_explore.txt", 'w')
    f.write("         NAME (hkl)     2*dspacing   emin  DarwinWidth\n")
    for i in range(n):
        f.write("%15s (%d%d%d)   %10.5g   %g  %g\n" % (
            mylist[i],
            H[i],K[i],L[i],
            2 * dspacing[i],
            emin[i],
            2e6*darwin_width[i],
        ))




