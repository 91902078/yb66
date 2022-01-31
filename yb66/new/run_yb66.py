import numpy
import os
import matplotlib.pylab as plt
from srxraylib.plot.gol import plot
from run_diff_pat import run_diff_pat

import xraylib
from xoppylib.decorators.xraylib_decorated import XraylibDecorated
from xoppylib.decorators.dabax_decorated import DabaxDecorated

from xoppylib.crystals.tools import bragg_calc, bragg_calc2

from check_FH import check_structure_factor

if __name__ == "__main__":

    xrl = XraylibDecorated()
    dx = DabaxDecorated()

    do_plot = 1


    #
    # YB66
    #

    if False:
        descriptor = 'YB66'
        SCANFROM = 0  # in microradiants
        SCANTO = 100  # in microradiants
        TEMPER = 1.0
        ENERGY = 8040.0
        SCANPOINTS = 200

        print("Using crystal descriptor: ", descriptor)
        bragg_dictionary = bragg_calc2(descriptor=descriptor,
                                       hh=4, kk=0, ll=0,
                                       temper=1.0,
                                       emin=ENERGY - 100.0, emax=ENERGY + 100.0,
                                       estep=(SCANTO - SCANFROM) / SCANPOINTS, fileout="xcrystal.bra",
                                       material_constants_library=dx)

        run_diff_pat(
            MOSAIC=0,
            GEOMETRY=0,
            SCAN=2,
            UNIT=1,
            SCANFROM=SCANFROM,
            SCANTO=SCANTO,
            SCANPOINTS=SCANPOINTS,
            ENERGY=ENERGY,
            ASYMMETRY_ANGLE=0.0,
            THICKNESS=0.7,
            MOSAIC_FWHM=0.1,
            RSAG=125.0,
            RMER=1290.0,
            ANISOTROPY=0,
            POISSON=0.22,
            CUT="2 -1 -1 ; 1 1 1 ; 0 0 0",
            FILECOMPLIANCE="mycompliance.dat")

        a = numpy.loadtxt("diff_pat.dat", skiprows=5)

        #
        # plot
        #
        plot(a[:, 0], a[:, -1], show=0,
             xtitle=r"$\theta-\theta_B$ [$\mu$rad]", ytitle="reflectivity", title="YB66 400 8040eV")

        plt.savefig("yb66_profile.png")
        print("File yb66_profile.png written to disk.")
        plt.show()



    #
    # method 2 (new)
    #
    if True:
        #
        # script to make the calculations (created by XOPPY:crystal)
        #

        import numpy
        from xoppylib.crystals.tools import run_diff_pat as run_diff_pat_new
        from xoppylib.crystals.tools import bragg_calc2
        import xraylib
        from dabax.dabax_xraylib import DabaxXraylib

        descriptor = 'YB66'
        SCANFROM = 0  # in microradiants
        SCANTO = 100  # in microradiants
        TEMPER = 1.0
        ENERGY = 8040.0
        SCANPOINTS = 200

        dx = DabaxXraylib()

        #
        # compute
        #

        print("Using crystal descriptor: ", descriptor)
        bragg_dictionary = bragg_calc2(descriptor=descriptor,
                                       hh=4, kk=0, ll=0,
                                       temper=1.0,
                                       emin=ENERGY - 100.0, emax=ENERGY + 100.0,
                                       estep=(SCANTO - SCANFROM) / SCANPOINTS, fileout="xcrystal.bra",
                                       material_constants_library=dx)


        run_diff_pat_new(
            bragg_dictionary,
            preprocessor_file="xcrystal.bra",
            MOSAIC=0,
            GEOMETRY=0,
            SCAN=2,
            UNIT=1,
            SCANFROM=SCANFROM,
            SCANTO=SCANTO,
            SCANPOINTS=SCANPOINTS,
            ENERGY=ENERGY,
            ASYMMETRY_ANGLE=0.0,
            THICKNESS=0.7,
            MOSAIC_FWHM=0.1,
            RSAG=125.0,
            RMER=1290.0,
            ANISOTROPY=0,
            POISSON=0.22,
            CUT="2 -1 -1 ; 1 1 1 ; 0 0 0",
            FILECOMPLIANCE="mycompliance.dat",
        )

        # import os
        # command = "..\diff_pat.exe < xoppy.inp"
        # # print("Running command '%s' in directory: %s " % (command, locations.home_bin_run()))
        # print("\n--------------------------------------------------------\n")
        # os.system(command)
        # print("\n--------------------------------------------------------\n")


        #
        # example plot
        #
        from srxraylib.plot.gol import plot

        data = numpy.loadtxt("diff_pat.dat", skiprows=5)
        plot(data[:, 0], data[:, -1])

        #
        # end script
        #
