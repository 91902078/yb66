#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy

def run_004(ISTAR1=0):
    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.BENER = 3.0
    oe0.EPSI_X = 1.25e-05
    oe0.EPSI_Z = 3.927e-07
    oe0.FDISTR = 4
    oe0.FSOURCE_DEPTH = 4
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.002
    oe0.HDIV2 = 0.002
    oe0.ISTAR1 = ISTAR1
    oe0.NCOL = 0
    oe0.NPOINT = 50000
    oe0.N_COLOR = 0
    oe0.PH1 = 1385.0
    oe0.PH2 = 1386.35
    oe0.POL_DEG = 0.0
    oe0.R_ALADDIN = 1191.3003399934003
    oe0.R_MAGNET = 11.913003399934002
    oe0.SIGDIX = 0.0
    oe0.SIGDIZ = 0.0
    oe0.SIGMAX = 0.0807
    oe0.SIGMAY = 0.0
    oe0.SIGMAZ = 0.0231
    oe0.VDIV1 = 0.0003
    oe0.VDIV2 = 0.0003
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0

    oe1.DUMMY = 1.0
    oe1.FMIRR = 3
    oe1.FWRITE = 1
    oe1.F_DEFAULT = 0
    oe1.SIMAG = 1400.0
    oe1.SSOUR = 1400.0
    oe1.THETA = 89.0
    oe1.T_IMAGE = 1100.0
    oe1.T_INCIDENCE = 89.0
    oe1.T_REFLECTION = 89.0
    oe1.T_SOURCE = 1400.0

    oe2.ALPHA = 180.0
    oe2.DUMMY = 1.0
    oe2.FILE_REFL = b'/nobackup/gurb1/srio/Oasys/YB66_004_1300_1500_test.dat'
    oe2.FWRITE = 1
    oe2.F_CENTRAL = 1
    oe2.F_CRYSTAL = 1
    oe2.PHOT_CENT = 1385.6
    oe2.R_LAMBDA = 5000.0
    oe2.T_INCIDENCE = 40.1985084472
    oe2.T_REFLECTION = 40.1985084472
    oe2.T_SOURCE = 0.0

    oe3.ALPHA = 180.0
    oe3.DUMMY = 1.0
    oe3.FILE_REFL = b'/nobackup/gurb1/srio/Oasys/YB66_004_1300_1500_test.dat'
    oe3.FWRITE = 1
    oe3.F_CENTRAL = 1
    oe3.F_CRYSTAL = 1
    oe3.PHOT_CENT = 1385.6
    oe3.R_LAMBDA = 5000.0
    oe3.T_IMAGE = 280.0
    oe3.T_INCIDENCE = 40.1985084472
    oe3.T_REFLECTION = 40.1985084472
    oe3.T_SOURCE = 0.0

    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    #
    # run optical element 1
    #
    print("    Running optical element: %d" % (1))
    if iwrite:
        oe1.write("start.01")

    beam.traceOE(oe1, 1)

    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")

    #
    # run optical element 2
    #
    print("    Running optical element: %d" % (2))
    if iwrite:
        oe2.write("start.02")

    beam.traceOE(oe2, 2)

    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")

    #
    # run optical element 3
    #
    print("    Running optical element: %d" % (3))
    if iwrite:
        oe3.write("start.03")

    beam.traceOE(oe3, 3)

    if iwrite:
        oe3.write("end.03")
        beam.write("star.03")

    # Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    return beam

def run_006(ISTAR1=0):
    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.BENER = 3.0
    oe0.EPSI_X = 1.25e-05
    oe0.EPSI_Z = 3.927e-07
    oe0.FDISTR = 4
    oe0.FSOURCE_DEPTH = 4
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.002
    oe0.HDIV2 = 0.002
    oe0.ISTAR1 = ISTAR1
    oe0.NCOL = 0
    oe0.NPOINT = 50000
    oe0.N_COLOR = 0
    oe0.PH1 = 2079.2
    oe0.PH2 = 2081.0
    oe0.POL_DEG = 0.0
    oe0.R_ALADDIN = 1191.3003399934003
    oe0.R_MAGNET = 11.913003399934002
    oe0.SIGDIX = 0.0
    oe0.SIGDIZ = 0.0
    oe0.SIGMAX = 0.0807
    oe0.SIGMAY = 0.0
    oe0.SIGMAZ = 0.0231
    oe0.VDIV1 = 0.0003
    oe0.VDIV2 = 0.0003
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0

    oe1.DUMMY = 1.0
    oe1.FMIRR = 3
    oe1.FWRITE = 1
    oe1.F_DEFAULT = 0
    oe1.SIMAG = 1400.0
    oe1.SSOUR = 1400.0
    oe1.THETA = 89.0
    oe1.T_IMAGE = 1100.0
    oe1.T_INCIDENCE = 89.0
    oe1.T_REFLECTION = 89.0
    oe1.T_SOURCE = 1400.0

    oe2.ALPHA = 180.0
    oe2.DUMMY = 1.0
    oe2.FILE_REFL = b'/nobackup/gurb1/srio/Oasys/YB66_006_2.dat'
    oe2.FWRITE = 1
    oe2.F_CENTRAL = 1
    oe2.F_CRYSTAL = 1
    oe2.PHOT_CENT = 2080.0
    oe2.R_LAMBDA = 5000.0
    oe2.T_INCIDENCE = 40.2678363404
    oe2.T_REFLECTION = 40.2678363404
    oe2.T_SOURCE = 0.0

    oe3.ALPHA = 180.0
    oe3.DUMMY = 1.0
    oe3.FILE_REFL = b'/nobackup/gurb1/srio/Oasys/YB66_006_2.dat'
    oe3.FWRITE = 1
    oe3.F_CENTRAL = 1
    oe3.F_CRYSTAL = 1
    oe3.PHOT_CENT = 2080.0
    oe3.R_LAMBDA = 5000.0
    oe3.T_IMAGE = 280.0
    oe3.T_INCIDENCE = 40.2678363404
    oe3.T_REFLECTION = 40.2678363404
    oe3.T_SOURCE = 0.0

    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    #
    # run optical element 1
    #
    print("    Running optical element: %d" % (1))
    if iwrite:
        oe1.write("start.01")

    beam.traceOE(oe1, 1)

    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")

    #
    # run optical element 2
    #
    print("    Running optical element: %d" % (2))
    if iwrite:
        oe2.write("start.02")

    beam.traceOE(oe2, 2)

    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")

    #
    # run optical element 3
    #
    print("    Running optical element: %d" % (3))
    if iwrite:
        oe3.write("start.03")

    beam.traceOE(oe3, 3)

    if iwrite:
        oe3.write("end.03")
        beam.write("star.03")

    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    return beam


if __name__ == "__main__":
    from srxraylib.plot.gol import plot

    npoints = 101
    nbins = 101


    ###########################################################################################################
    FWHM1= numpy.zeros(npoints)
    INTENSITY1 = numpy.zeros(npoints)

    for i in range(npoints):
        # beam = run_004(ISTAR1=7123 + i * 5)
        # tkt = beam.histo1(11, nbins=101, ref=23, nolost=1, xrange=[1385, 1386.35], calculate_widths=1)

        beam = run_006(ISTAR1=7123 + i * 5)
        tkt = beam.histo1(11, nbins=nbins, ref=23, nolost=1, xrange=[2079.2, 2081.0], calculate_widths=1)
        # Shadow.ShadowTools.histo1(beam, 11, nbins=101, ref=23, nolost=1)
        FWHM1[i] = tkt["fwhm"]
        INTENSITY1[i] = tkt["intensity"]
        x = tkt["bin_path"]
        y = tkt["histogram_path"]
        if i == 0:
            yacc = y
        else:
            yacc += y

    int1 = npoints * INTENSITY1.mean()
    error1 = npoints * INTENSITY1.std()

    ###########################################################################################################
    FWHM2= numpy.zeros(npoints)
    INTENSITY2 = numpy.zeros(npoints)

    for i in range(npoints):
        beam = run_004(ISTAR1=7123 + i * 5)
        tkt = beam.histo1(11, nbins=nbins, ref=23, nolost=1, xrange=[1385, 1386.35], calculate_widths=1)
        # Shadow.ShadowTools.histo1(beam, 11, nbins=101, ref=23, nolost=1)
        FWHM2[i] = tkt["fwhm"]
        INTENSITY2[i] = tkt["intensity"]
        x = tkt["bin_path"]
        y = tkt["histogram_path"]
        if i == 0:
            yacc = y
        else:
            yacc += y

    int2 = npoints * INTENSITY2.mean()
    error2 = npoints * INTENSITY2.std()




    print("INTENSITY1: ", int1, error1 )
    print("FWHM1: ", FWHM1.mean(), FWHM1.std())
    print("INTENSITY2: ", int2, error2 )
    print("FWHM2: ", FWHM2.mean(), FWHM2.std())


    # plot(x,yacc)
    # print()
    #
    # int1 = 20437.81
    # int2 = 9214.87
    # error1 = 414.1
    # error2 = 165.30
    mean = (1.01 * 0.086 * 0.59 * 0.027 * int1) / (1.0 * 0.097 * 0.77 * 0.045 * int2)

    print( mean, "+/-", mean * (error1/int1 + error2/int2))