import os
import numpy

from orangecontrib.xoppy.util.xoppy_util import locations, XoppyPhysics
from xoppy_dabax_util import bragg_calc2

from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_calc, crystal_fh
import platform


# copied and adapted from xcrystal.py
def run_diff_pat(
    # bragg_dictionary,
    # CRYSTAL_MATERIAL = self.CRYSTAL_MATERIAL
    # MILLER_INDEX_H = self.MILLER_INDEX_H
    # MILLER_INDEX_K = self.MILLER_INDEX_K
    # MILLER_INDEX_L = self.MILLER_INDEX_L
    # TEMPER = self.TEMPER
    MOSAIC = 0,
    GEOMETRY = 0,
    SCAN = 2,
    UNIT = 1,
    SCANFROM = -100.0,
    SCANTO = 100.0,
    SCANPOINTS = 200,
    ENERGY = 8000.0,
    ASYMMETRY_ANGLE = 0.0,
    THICKNESS = 0.7,
    MOSAIC_FWHM = 0.1,
    RSAG = 125.0,
    RMER = 1290.0,
    ANISOTROPY = 0,
    POISSON = 0.22,
    CUT = "2 -1 -1 ; 1 1 1 ; 0 0 0",
    FILECOMPLIANCE = "mycompliance.dat",
    ):

    for file in ["xoppy.inp", "diff_pat.dat", "diff_pat.gle", "diff_pat.par", "diff_pat.xop"]: #, "xcrystal.bra"]:
        try:
            os.remove(os.path.join(locations.home_bin_run(), file))
        except:
            pass

    if (GEOMETRY == 1) or (GEOMETRY == 3):
        if ASYMMETRY_ANGLE == 0.0:
            print(
                "xoppy_calc_xcrystal: WARNING: In xcrystal the asymmetry angle is the angle between Bragg planes and crystal surface," +
                "in BOTH Bragg and Laue geometries.")

    with open("./xoppy.inp", "wt") as f:
        f.write("xcrystal.bra\n")
        f.write("%d\n" % MOSAIC)
        f.write("%d\n" % GEOMETRY)

        if MOSAIC == 1:
            f.write("%g\n" % MOSAIC_FWHM)
            f.write("%g\n" % THICKNESS)
        else:
            f.write("%g\n" % THICKNESS)
            f.write("%g\n" % ASYMMETRY_ANGLE)

        scan_flag = 1 + SCAN

        f.write("%d\n" % scan_flag)

        f.write("%19.9f\n" % ENERGY)

        if scan_flag <= 3:
            f.write("%d\n" % UNIT)

        f.write("%g\n" % SCANFROM)
        f.write("%g\n" % SCANTO)
        f.write("%d\n" % SCANPOINTS)

        if MOSAIC > 1:  # bent
            f.write("%g\n" % RSAG)
            f.write("%g\n" % RMER)
            f.write("0\n")

            f.write("%d\n" % ANISOTROPY)

            if ANISOTROPY == 0:
                f.write("%g\n" % POISSON)


    if platform.system() == "Windows":
        command = "diff_pat.exe < xoppy.inp"
    elif platform.system() == "Darwin":
        command = "./diff_pat.darwin < xoppy.inp"
    elif platform.system() == "Linux":
        command = "./diff_pat.linux  < xoppy.inp"

    print("Running command '%s' in directory: %s " % (command, locations.home_bin_run()))
    print("\n--------------------------------------------------------\n")
    os.system(command)
    print("\n--------------------------------------------------------\n")

    # show calculated parameters in standard output
    txt_info = open("diff_pat.par").read()
    for line in txt_info:
        print(line, end="")


if __name__ == "__main__":

    #
    # old code
    #
    descriptor = "Si"

    dic1a = bragg_calc(descriptor=descriptor,hh=1,kk=1,ll=1,temper=1.0,emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra")
    print("KEYS: ",dic1a.keys())
    print(dic1a)

    dic1b = crystal_fh(dic1a,8000.0)
    print(dic1b["info"])
    print("KEYS: ",dic1b.keys())

    run_diff_pat()

    a1 = numpy.loadtxt("diff_pat.dat",skiprows=5)

    #
    # New code
    #

    dic2a = bragg_calc2(descriptor=descriptor,hh=1,kk=1,ll=1,temper=1.0,emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra")
    print("KEYS: ",dic2a.keys())
    print(dic2a)

    dic2b = crystal_fh(dic1a,8000.0)
    print(dic2b["info"])
    print("KEYS: ",dic2b.keys())


    for key in dic1a.keys():
        print(key,dic1a[key],dic2a[key])


    run_diff_pat()

    a2 = numpy.loadtxt("diff_pat.dat",skiprows=5)

    #
    # comparison
    #

    from srxraylib.plot.gol import plot
    plot(a1[:, 0], a1[:, -1],
         a2[:, 0], a2[:, -1],
         legend=["OLD code","NEW code"],
         title=descriptor)
