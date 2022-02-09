
#
# creates a list of high-dspacing crystal reflections from the list of crystals
# in DABAX files crystals.dat (original DABAX and xraylib list) and Crystals_xrayserver.dat
# (from Sergey Stepanov list).
#
# Run before table1_explore.py to create the list of selected reflections
# in file table1_explore.txt. It is edited by hand to comment some duplicated crystals
#
#

import numpy

from srxraylib.plot.gol import plot, set_qt
from xoppylib.decorators.dabax_decorated import DabaxDecorated
import scipy.constants as codata

from crystalpy.diffraction.DiffractionSetupDabax import DiffractionSetupDabax


def run(descriptor, hh, kk, ll, material_constants_library, energy, do_plot=1):
    import os
    import xraylib

    from xoppylib.crystals.tools import bragg_calc, bragg_calc2
    from xoppylib.crystals.tools import run_diff_pat



    os.system("rm -f xcrystal.bra xoppy.inp diff_pat.dat")

    dic1a = bragg_calc2(descriptor=descriptor,hh=hh,kk=kk,ll=ll,temper=1.0,
                       emin=energy-100,emax=energy+100,estep=5.0,fileout="xcrystal.bra",
                       material_constants_library=material_constants_library)

    print("KEYS: ",dic1a.keys())


    run_diff_pat(
        dic1a,
        preprocessor_file="xcrystal.bra",
        MOSAIC = 0,
        GEOMETRY = 0,
        SCAN = 2,
        UNIT = 1,
        SCANFROM = -1000.0,
        SCANTO = 2000.0,
        SCANPOINTS = 200,
        ENERGY = energy,
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

    if do_plot:
        plot(aa[:,0],aa[:,-1], title="%s (%d%d%d) @ %g eV" % (descriptor,h,k,l,energy))
    return aa[:,-1].max()



if __name__ == "__main__":

    set_qt()

    dx0 = DabaxDecorated(file_Crystals="Crystals.dat")
    dx1 = DabaxDecorated(file_Crystals="Crystals_xrayserver.dat")

    do_plot = 1



    with open('table1_explore.txt') as f:
        lines = f.readlines()
    f = open("table1.txt", 'w')

    for line in lines:
        if line[0] != "#":
            var = line.split("  ")
            var = ' '.join(var).split()
            mylib_flag = int(var[0])
            if  mylib_flag == 0:
                mylib = dx0
            else:
                mylib = dx1

            descriptor = var[1]
            h = int(var[2])
            k = int(var[3])
            l = int(var[4])



            c0 = mylib.Crystal_GetCrystal(descriptor)
            ds = dx0.Crystal_dSpacing(c0,h,k,l)

            lambda_emin = 2 * ds
            emin = codata.h * codata.c / codata.e / (lambda_emin * 1e-10)

            lambda_emax = 2 * ds * numpy.sin(5 * numpy.pi/180)
            emax = codata.h * codata.c / codata.e / (lambda_emax * 1e-10)


            cpy = DiffractionSetupDabax(geometry_type=0, crystal_name=descriptor, thickness=1e-5,
                 miller_h=h, miller_k=k, miller_l=l,
                 asymmetry_angle=0.0,
                 azimuthal_angle=0.0,
                                        dabax=mylib)

            darwin_width = 2 * cpy.darwinHalfwidthS(emin*1.1)
            bragg_angle = numpy.arcsin(lambda_emin / 1.1 / 2 / ds)
            DE_E = darwin_width  / numpy.tan(bragg_angle)

            peak_intensity = run(descriptor,h,k,l,mylib,emin*1.1,do_plot=0)


            if (round(100 * peak_intensity) > 1):
                print("%d  %20s (%d%d%d)  %5.2f  %3.1f-%3.1f (%d) theta=%g DE/E=%4.2f peak=%d" %
                      (mylib_flag, descriptor, h, k, l,
                       (2 * ds), emin * 1e-3, emax * 1e-3, darwin_width * 1e6,
                       bragg_angle * 180 / numpy.pi, DE_E * 1e3,
                       100 * peak_intensity
                       ))

                f.write("%d  %20s (%d%d%d)  %5.2f  %3.1f-%3.1f (%5d) %4.2f %d\n" %
                      (mylib_flag, descriptor, h, k, l,
                                                     (2*ds), emin*1e-3, emax*1e-3, darwin_width*1e6,
                                                     DE_E*1e3,
                                                     round(100 * peak_intensity)
                                                                      ))


    f.close()




