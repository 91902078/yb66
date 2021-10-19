


import numpy

from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_calc as bragg_calc_old

from xoppy_xraylib_util import bragg_calc

from xoppy_dabax_util import bragg_calc2

import os

from run_diff_pat import run_diff_pat

if __name__ == "__main__":


    for descriptor in ["Muscovite", "Si"]:
        #
        # old code muscovite
        #
        dic1a = bragg_calc_old(descriptor=descriptor,hh=1,kk=1,ll=1,temper=1.0,emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra")
        os.system("cp xcrystal.bra xcrystal.bra.old")

        run_diff_pat(
            MOSAIC=0,
            GEOMETRY=0,
            SCAN=2,
            UNIT=1,
            SCANFROM=-100.0,
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

        a1 = numpy.loadtxt("diff_pat.dat", skiprows=5)


        #
        # New code muscovite
        #

        dic2a = bragg_calc(descriptor=descriptor,hh=1,kk=1,ll=1,temper=1.0,emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra")
        print("KEYS: ",dic2a.keys())
        os.system("cp xcrystal.bra xcrystal.bra.new")

        run_diff_pat(
            MOSAIC=0,
            GEOMETRY=0,
            SCAN=2,
            UNIT=1,
            SCANFROM=-100.0,
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
        # Xiaojiang code
        #


        dic3a = bragg_calc2(descriptor=descriptor, hh=1, kk=1, ll=1, temper=1.0, emin=7900.0, emax=8100.0, estep=5.0,
                            fileout="xcrystal.bra")
        print("KEYS: ", dic3a.keys())
        os.system("cp xcrystal.bra xcrystal.bra.xj")

        run_diff_pat(
            MOSAIC=0,
            GEOMETRY=0,
            SCAN=2,
            UNIT=1,
            SCANFROM=-100.0,
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

        a3 = numpy.loadtxt("diff_pat.dat", skiprows=5)





        for key in dic1a.keys():
            if key != "info":
                print(">>>", key,dic1a[key],dic2a[key],dic3a[key])
                # tmp = numpy.abs(dic1b[key] - dic2b[key])
                # if tmp.size == 1:
                #     assert (tmp < 1e-6)
                # else:
                #     assert(tmp.sum() < 1e-6)

        #
        # comparison
        #

        print(a1.shape, a2.shape, a3.shape)
        from srxraylib.plot.gol import plot, set_qt
        set_qt()

        plot(a1[:, 0], a1[:, -1],
             a2[:, 0], a2[:, -1],
             a3[:, 0], a3[:, -1],
             title="Crystal: " + descriptor,
             marker=[None, "o", "+"],
             legend=["current xoppy", "local xoppy_xraylib_util", "bragg_calc2 in xoppy_dabax_util"])