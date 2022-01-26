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

    if True:
        #
        # f0
        #
        Y0_xrl = xrl.f0_calc(0, "Y", 0, 3, 500)
        Y0_dbx =  dx.f0_calc(0, "Y", 0, 3, 500)
        Y3_dbx =  dx.f0_calc(0, "Y", 0, 3, 500, charge=3.0)
        B0_xrl = xrl.f0_calc(0, "B", 0, 3, 500)
        B0_dbx =  dx.f0_calc(0, "B", 0, 3, 500)
        Bf_dbx =  dx.f0_calc(0, "B", 0, 3, 500, charge=-0.045)

        if do_plot:
            plot(Y0_xrl["data"][0, :], Y0_xrl["data"][1,:],
                 Y0_dbx["data"][0, :], Y0_dbx["data"][1,:],
                 Y3_dbx["data"][0, :], Y3_dbx["data"][1, :],
                 B0_xrl["data"][0, :], B0_xrl["data"][1, :],
                 B0_dbx["data"][0, :], B0_dbx["data"][1, :],
                 Bf_dbx["data"][0, :], Bf_dbx["data"][1, :],
                 #
                 # linestyle=[None,None,None],
                 # marker=[None,'+',None,'+'],
                 # color=['r','r','b','b'],
                 legend=['Y xraylib',r'Y$^{+0}$ dabax',r'Y$^{+3}$ dabax',
                         'Y xraylib',r'B$^{+0}$ dabax',r'B$^{-0.045}$ dabax'],
                 xtitle=r'q=sin$\theta$/$\lambda$',ytitle='f$_0$ [electron units]',
                 show=0)

            plt.savefig("f0.png")
            print("File f0.png written to disk.")
            plt.show()


#
# table 3
#
    if True:


        # test Muscovite
        F200 = numpy.abs( check_structure_factor(descriptor="YB66", hh=2, kk=0, ll=0, energy=8040.0, do_assert=0, models=[0,0,1], material_constants_library=dx) )
        F400 = numpy.abs( check_structure_factor(descriptor="YB66", hh=4, kk=0, ll=0, energy=8040.0, do_assert=0, models=[0,0,1], material_constants_library=dx) )
        F600 = numpy.abs( check_structure_factor(descriptor="YB66", hh=6, kk=0, ll=0, energy=8040.0, do_assert=0, models=[0,0,1], material_constants_library=dx) )
        F800 = numpy.abs( check_structure_factor(descriptor="YB66", hh=8, kk=0, ll=0, energy=8040.0, do_assert=0, models=[0,0,1], material_constants_library=dx) )

        F200_calc_ref = 8.8 * 4
        F400_calc_ref = 137.9 * 4
        F600_calc_ref = 11.2 * 4
        F800_calc_ref = 4.8 * 4

        F200_obs_ref = 6.9 * 4
        F400_obs_ref = 150.4 * 4
        F600_obs_ref = 11.9 * 4
        F800_obs_ref = 7.2 * 4

        print("F_200 this work: F=%g, calculated(ref): Fc=%g, observed(ref): Fo=%g, |F/Fc|=%g" % (F200, F200_calc_ref, F200_obs_ref, F200/F200_calc_ref))
        print("F_400 this work: F=%g, calculated(ref): Fc=%g, observed(ref): Fo=%g, |F/Fc|=%g" % (F400, F400_calc_ref, F400_obs_ref, F400/F400_calc_ref))
        print("F_600 this work: F=%g, calculated(ref): Fc=%g, observed(ref): Fo=%g, |F/Fc|=%g" % (F600, F600_calc_ref, F600_obs_ref, F600/F600_calc_ref))
        print("F_800 this work: F=%g, calculated(ref): Fc=%g, observed(ref): Fo=%g, |F/Fc|=%g" % (F800, F800_calc_ref, F800_obs_ref, F800/F800_calc_ref))

    #
    # no prototype:
    # F_200 this work: F=39.2492, calculated(ref): Fc=35.2, observed(ref): Fo=27.6, |F/Fc|=1.11503
    # F_400 this work: F=563.74, calculated(ref): Fc=551.6, observed(ref): Fo=601.6, |F/Fc|=1.02201
    # F_600 this work: F=38.893, calculated(ref): Fc=44.8, observed(ref): Fo=47.6, |F/Fc|=0.868148
    # F_800 this work: F=24.8411, calculated(ref): Fc=19.2, observed(ref): Fo=28.8, |F/Fc|=1.29381
    #
    # with prototype ** DEFAULT **
    # F_200 this work: F=39.2639, calculated(ref): Fc=35.2, observed(ref): Fo=27.6, |F/Fc|=1.11545
    # F_400 this work: F=564.591, calculated(ref): Fc=551.6, observed(ref): Fo=601.6, |F/Fc|=1.02355
    # F_600 this work: F=39.0252, calculated(ref): Fc=44.8, observed(ref): Fo=47.6, |F/Fc|=0.871098
    # F_800 this work: F=24.9914, calculated(ref): Fc=19.2, observed(ref): Fo=28.8, |F/Fc|=1.30163
    #
    # no temperature factor (identical results useing no prototypical or with prototypical)
    # F_200 this work: F=39.3788, calculated(ref): Fc=35.2, observed(ref): Fo=27.6, |F/Fc|=1.11871
    # F_400 this work: F=571.224, calculated(ref): Fc=551.6, observed(ref): Fo=601.6, |F/Fc|=1.03558
    # F_600 this work: F=40.0644, calculated(ref): Fc=44.8, observed(ref): Fo=47.6, |F/Fc|=0.894294
    # F_800 this work: F=26.1867, calculated(ref): Fc=19.2, observed(ref): Fo=28.8, |F/Fc|=1.36389
    #
    # isotropic temperature factor
    # F_200 this work: F=39.2584, calculated(ref): Fc=35.2, observed(ref): Fo=27.6, |F/Fc|=1.1153
    # F_400 this work: F=564.274, calculated(ref): Fc=551.6, observed(ref): Fo=601.6, |F/Fc|=1.02298
    # F_600 this work: F=38.976, calculated(ref): Fc=44.8, observed(ref): Fo=47.6, |F/Fc|=0.87
    # F_800 this work: F=24.9354, calculated(ref): Fc=19.2, observed(ref): Fo=28.8, |F/Fc|=1.29872
    #
    #
    #
    #

    #
    # muscovite profile
    #
    if True:

        descriptor = 'Muscovite'
        #
        # old code
        #
        os.system("rm -f xcrystal.bra xoppy.inp")

        dic1a = bragg_calc(descriptor=descriptor,hh=1,kk=1,ll=1,temper=1.0,
                           emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra",
                           material_constants_library=xraylib)

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
             title=descriptor,
             legend=["XOPPY bragg_calc code","NEW bragg_calc2 code"], show=0,
             xtitle=r"$\theta-\theta_B$ [$\mu$rad]", ytitle="reflectivity")

        plt.savefig("muscovite.png")
        print("File muscovite.png written to disk.")
        plt.show()



    #
    # YB66
    #

    if True:
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