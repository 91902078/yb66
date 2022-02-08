import numpy
import os
import matplotlib.pylab as plt
from srxraylib.plot.gol import plot, set_qt

import xraylib
from dabax.dabax_xraylib import DabaxXraylib
# from xoppylib.decorators.xraylib_decorated import XraylibDecorated
# from xoppylib.decorators.dabax_decorated import DabaxDecorated

from xoppylib.crystals.tools import bragg_calc, bragg_calc2
from xoppylib.crystals.tools import run_diff_pat



if __name__ == "__main__":

    set_qt()
    do_plot = 1


    #
    # muscovite profile
    #


    #
    # old code
    #
    os.system("rm -f xcrystal.bra xoppy.inp")

    dic1a = bragg_calc2(descriptor='Muscovite',hh=1,kk=1,ll=1,temper=1.0,
                       emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra",
                       material_constants_library=xraylib)

    print("KEYS: ",dic1a.keys())
    print(dic1a)
    os.system("cp xcrystal.bra xcrystal.bra.old")

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

    a1 = numpy.loadtxt("diff_pat.dat",skiprows=5)



    #
    # New code
    #
    dic2a = bragg_calc2(descriptor='Muscovite', hh=1, kk=1, ll=1, temper=1.0, emin=7900.0, emax=8100.0, estep=5.0,
                        fileout="xcrystal.bra",
                        material_constants_library=DabaxXraylib(), )
    print("KEYS: ", dic2a.keys())
    print(dic2a)
    os.system("cp xcrystal.bra xcrystal.bra.new")

    run_diff_pat(
        dic2a,
        preprocessor_file="xcrystal.bra",
        MOSAIC=0,
        GEOMETRY=0,
        SCAN=2,
        UNIT=1,
        SCANFROM=25.0,
        SCANTO=100.0,
        SCANPOINTS=200,
        ENERGY=8040.0,
        ASYMMETRY_ANGLE=0.0,
        THICKNESS=0.7,
        MOSAIC_FWHM=0.1,
        RSAG=125.0,
        RMER=1290.0,
        ANISOTROPY=0,
        POISSON=0.22,
        CUT="2 -1 -1 ; 1 1 1 ; 0 0 0",
        FILECOMPLIANCE="mycompliance.dat")

    a2 = numpy.loadtxt("diff_pat.dat", skiprows=5)


    #
    # DABAX Crystals_xrayserver.dat
    #
    dic3a = bragg_calc2(descriptor="Mica",kk=1,ll=1,temper=1.0,emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra",
                        material_constants_library=DabaxXraylib(file_Crystals="stepanov/Crystals_xrayserver.dat"),)
    print("KEYS: ",dic2a.keys())
    print(dic3a)
    os.system("cp xcrystal.bra xcrystal.bra.3a")



    run_diff_pat(
        dic3a,
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

    a3 = numpy.loadtxt("diff_pat.dat",skiprows=5)





    #
    # comparison
    #

    plot(a1[:, 0], a1[:, -1],
         a2[:, 0], a2[:, -1],
         a3[:, 0], a3[:, -1],
         linestyle=[None,'--',''],
         marker=[None,"+",'o'],
         title="Muscovite",
         legend=["using xraylib","using DABAX Crystals.dat","using DABAX Crystals_xrayserver.dat"], show=0,
         xtitle=r"$\theta-\theta_B$ [$\mu$rad]", ytitle="reflectivity")

    plt.savefig("muscovite.png")
    print("File muscovite.png written to disk.")
    plt.show()

    #
    # crystals in Crystals_xrayserver.dat
    #
    from xoppylib.decorators.dabax_decorated import DabaxDecorated
    dx = DabaxDecorated(file_Crystals="stepanov/Crystals_xrayserver.dat")
    print(dx.Crystal_GetCrystalsList())