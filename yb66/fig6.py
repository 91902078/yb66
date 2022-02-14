def run004():
    #
    # script to calculate crystal diffraction profiles (created by XOPPY:crystal)
    #

    import numpy
    from xoppylib.crystals.tools import bragg_calc2, run_diff_pat
    import xraylib
    from dabax.dabax_xraylib import DabaxXraylib

    #
    # run bragg_calc (preprocessor) and create file xcrystal.bra
    #
    bragg_dictionary = bragg_calc2(
        descriptor="YB66",
        hh=0,
        kk=0,
        ll=4,
        temper=1.0,
        emin=1285.6,
        emax=1485.6,
        estep=1.0,
        ANISO_SEL=2,
        fileout="xcrystal.bra",
        do_not_prototype=0,  # 0=use site groups (recommended), 1=use all individual sites
        verbose=False,
        material_constants_library=DabaxXraylib(),
    )

    #
    # run external (fortran) diff_pat (note that some parameters may not be used)
    #
    run_diff_pat(
        bragg_dictionary,
        preprocessor_file="xcrystal.bra",
        MOSAIC=0,
        GEOMETRY=0,
        SCAN=2,
        UNIT=1,
        SCANFROM=200.0,
        SCANTO=800.0,
        SCANPOINTS=500,
        ENERGY=1385.6,
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

    #
    # example plot
    #
    # from srxraylib.plot.gol import plot
    data = numpy.loadtxt("diff_pat.dat", skiprows=5)
    # plot(data[:, 0], data[:, -1])
    return data

    #
    # end script
    #

def run006():
    #
    # script to calculate crystal diffraction profiles (created by XOPPY:crystal)
    #

    import numpy
    from xoppylib.crystals.tools import bragg_calc2, run_diff_pat
    import xraylib
    from dabax.dabax_xraylib import DabaxXraylib

    #
    # run bragg_calc (preprocessor) and create file xcrystal.bra
    #
    bragg_dictionary = bragg_calc2(
        descriptor="YB66",
        hh=0,
        kk=0,
        ll=6,
        temper=1.0,
        emin=1980.0,
        emax=2180.0,
        estep=1.0,
        ANISO_SEL=2,
        fileout="xcrystal.bra",
        do_not_prototype=0,  # 0=use site groups (recommended), 1=use all individual sites
        verbose=False,
        material_constants_library=DabaxXraylib(),
    )

    #
    # run external (fortran) diff_pat (note that some parameters may not be used)
    #
    run_diff_pat(
        bragg_dictionary,
        preprocessor_file="xcrystal.bra",
        MOSAIC=0,
        GEOMETRY=0,
        SCAN=2,
        UNIT=1,
        SCANFROM=0.0,
        SCANTO=400.0,
        SCANPOINTS=500,
        ENERGY=2080.0,
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

    #
    # example plot
    #
    # from srxraylib.plot.gol import plot
    data = numpy.loadtxt("diff_pat.dat", skiprows=5)
    # plot(data[:, 0], data[:, -1])
    return data
    #
    # end script
    #

if __name__ == "__main__":
    import numpy
    from srxraylib.plot.gol import set_qt
    set_qt()

    r1 = run004()
    r2 = run006()

    from oasys.util.oasys_util import get_fwhm

    fwhm1, tmp, tmp = get_fwhm(r1[:, -1], r1[:, 0])
    fwhm2, tmp, tmp = get_fwhm(r2[:, -1], r2[:, 0])
    energy1 = 1385.6
    energy2 = 2080.0
    peak1 = r1[:, -1].max()
    peak2 = r2[:, -1].max()

    from srxraylib.plot.gol import plot

    p = plot(r1[:, 0], r1[:, -1], r2[:, 0], r2[:, -1], xtitle=r"$\theta-\theta_B [\mu rad]$", ytitle="Reflectivity",
             legend=["E=%g eV FWHM=%3.1f Peak=%3.2f" % (energy1, fwhm1, peak1),
                     "E=%g eV FWHM=%3.1f Peak=%3.2f" % (energy2, fwhm2, peak2)], show=0, yrange=[0, 0.3])
    import matplotlib.pylab as plt

    plt.savefig("fig6.png")
    plt.show()

