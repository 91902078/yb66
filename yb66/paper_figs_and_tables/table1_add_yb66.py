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

    do_plot = 1


#
# table 3
#

    list0 = dx0.Crystal_GetCrystalsList()


    n0 = len(list0)


    mylist = []
    for i in range(n0): mylist.append(list0[i])

    print(mylist)


    TXT = ""


    for i in range(n0):
        mlib = dx0
        # print(">>>>>", mlib)
        c0 = mlib.Crystal_GetCrystal(mylist[i])
        if mylist[i] == "YB66" or mylist[i] == "Si":
            for h in [6,5,4,3,2,1,0]:
                for k in [1,0]:
                    for l in [1,0]:
                        ds = mlib.Crystal_dSpacing(c0, h, k, l)
                        ee = codata.h * codata.c / codata.e / (2 * ds * 1e-10)

                        cpy = DiffractionSetupDabax(geometry_type=0, crystal_name=mylist[i], thickness=1e-5,
                             miller_h=h,miller_k=k, miller_l=l,
                             asymmetry_angle=0.0,
                             azimuthal_angle=0.0,
                             dabax=mlib)

                        FH = cpy.FH(ee * 1.1)

                        # print(mylist[i], h, k, l, ds, numpy.abs(FH) )

                        if ( (h+k+l) > 0):
                            if numpy.abs(FH) > 1:
                                print(">>>>>> FOUND!", mylist[i], h,k,l)
                                TXT += "%s %d %d %d " % (mylist[i], h,k,l)

                                lambda_emin = 2 * ds
                                emin0 = codata.h * codata.c / codata.e / (lambda_emin * 1e-10)

                                lambda_emax = 2 * ds * numpy.sin(5 * numpy.pi / 180)
                                emax0 = codata.h * codata.c / codata.e / (lambda_emax * 1e-10)

                                cpy = DiffractionSetupDabax(geometry_type=0, crystal_name=mylist[i], thickness=1e-5,
                                                            miller_h=h, miller_k=k, miller_l=l,
                                                            asymmetry_angle=0.0,
                                                            azimuthal_angle=0.0,
                                                            dabax=dx0)

                                darwin_width0 = 2 * cpy.darwinHalfwidthS(emin0 * 1.1)
                                bragg_angle0 = numpy.arcsin(lambda_emin / 1.1 / 2 / ds)
                                DE_E = darwin_width0 / numpy.tan(bragg_angle0)

                                peak_intensity = run(mylist[i], h, k, l, mlib, emin0 * 1.1, do_plot=0)


                                TXT += "2d: %g, emin: %g, emax: %g, DarwinWidth@emin*1.1: %d, DE/E: %g, E/DE: %d, R(p.c.): %d\n" % \
                                    (
                                        numpy.round(2 * ds, 2),
                                        numpy.round(emin0 * 1e-3, 1),
                                        numpy.round(emax0 * 1e-3, 1),
                                        int(darwin_width0 * 1e6),
                                        numpy.round(DE_E * 1e3,3),
                                        1 / (DE_E),
                                        int(100 * peak_intensity))



    print(TXT)





