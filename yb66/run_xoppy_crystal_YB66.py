import numpy
from crystal_util import bragg_calc2, crystal_fh2
from run_diff_pat import run_diff_pat
from srxraylib.plot.gol import plot

if __name__ == "__main__":

    descriptor = 'YB66'
    SCANFROM = 0 # in microradiants
    SCANTO = 100 # in microradiants
    MILLER_INDEX_H = 4
    MILLER_INDEX_K = 0
    MILLER_INDEX_L = 0
    TEMPER = 1.0
    ENERGY = 8040.0
    SCANPOINTS = 200


    print("Using crystal descriptor: ",descriptor)

    bragg_dictionary = bragg_calc2(descriptor=descriptor,
                                            hh=MILLER_INDEX_H,kk=MILLER_INDEX_K,ll=MILLER_INDEX_L,
                                            temper=TEMPER,
                                            emin=ENERGY-100.0,emax=ENERGY+100.0,
                                            estep=(SCANTO-SCANFROM)/SCANPOINTS,fileout="xcrystal.bra")


    run_diff_pat(
        MOSAIC = 0,
        GEOMETRY = 0,
        SCAN = 2,
        UNIT = 1,
        SCANFROM = SCANFROM,
        SCANTO = SCANTO,
        SCANPOINTS = SCANPOINTS,
        ENERGY = ENERGY,
        ASYMMETRY_ANGLE = 0.0,
        THICKNESS = 0.7,
        MOSAIC_FWHM = 0.1,
        RSAG = 125.0,
        RMER = 1290.0,
        ANISOTROPY = 0,
        POISSON = 0.22,
        CUT = "2 -1 -1 ; 1 1 1 ; 0 0 0",
        FILECOMPLIANCE = "mycompliance.dat")

    a = numpy.loadtxt("diff_pat.dat",skiprows=5)

    #
    # plot
    #
    plot(a[:, 0], a[:, -1])




