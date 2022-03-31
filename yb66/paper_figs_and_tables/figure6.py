
#
# test numeric f1f2 values with default DABAX f1f2_Windt.dat file (that uses Henke/CXRO data for E<30 keV)
#
import time
from dabax.dabax_xraylib import DabaxXraylib

f_W = DabaxXraylib()
print(f_W.info())
energy_in_eV_W, f1_W, f2_W = f_W.f1f2_extract("Y")

for i,energy in enumerate(energy_in_eV_W):
    if energy > 1300 and energy < 1400:
        print(">>>", energy, f1_W[i], f2_W[i])
    elif energy > 2000 and energy < 2100:
        print(">>>", energy, f1_W[i], f2_W[i])



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
import os
os.system("rm xcrystal.bra")

t0 = time.time()
bragg_dictionary = bragg_calc2(
    descriptor = "YB66",
    hh         = 0,
    kk         = 0,
    ll         = 4,
    temper     = 1.0,
    emin       = 1370.0,
    emax       = 1410.0,
    estep      = 0.4, #1.0,
    ANISO_SEL  = 2, #0,
    fileout    = "xcrystal.bra",
    do_not_prototype = 0,  # 0=use site groups (recommended), 1=use all individual sites
    verbose = False,
    material_constants_library = DabaxXraylib(),
    )

#
# run external (fortran) diff_pat (note that some parameters may not be used)
#
run_diff_pat(
    bragg_dictionary,
    preprocessor_file = "xcrystal.bra",
    MOSAIC             = 0,
    GEOMETRY           = 0,
    SCAN               = 2,
    UNIT               = 1,
    SCANFROM           = 200.0,
    SCANTO             = 800.0,
    SCANPOINTS         = 2000,
    ENERGY             = 1385.6,
    ASYMMETRY_ANGLE    = 0.0,
    THICKNESS          = 0.7,
    MOSAIC_FWHM        = 0.1,
    RSAG               = 125.0,
    RMER               = 1290.0,
    ANISOTROPY         = 0, #0,
    POISSON            = 0.22,
    CUT                = "2 -1 -1 ; 1 1 1 ; 0 0 0",
    FILECOMPLIANCE     = "mycompliance.dat",
    )

# #
# # example plot
# #
# from srxraylib.plot.gol import plot
data1 = numpy.loadtxt("diff_pat.dat", skiprows=5)
# plot(data[:,0], data[:,-1], title="FWHM: %g" % fwhm(data[:,0], data[:,-1]))

#
# end script
#


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
import os
os.system("rm xcrystal.bra")
bragg_dictionary = bragg_calc2(
    descriptor = "YB66",
    hh         = 0,
    kk         = 0,
    ll         = 6,
    temper     = 1.0,
    emin       = 2060.0,
    emax       = 2100.0,
    estep      = 0.25, #1.0,
    ANISO_SEL  = 2, #0,
    fileout    = "xcrystal.bra",
    do_not_prototype = 0,  # 0=use site groups (recommended), 1=use all individual sites
    verbose = False,
    material_constants_library = DabaxXraylib(),
    )

#
# run external (fortran) diff_pat (note that some parameters may not be used)
#
run_diff_pat(
    bragg_dictionary,
    preprocessor_file = "xcrystal.bra",
    MOSAIC             = 0,
    GEOMETRY           = 0,
    SCAN               = 2,
    UNIT               = 1,
    SCANFROM           = 0.0,
    SCANTO             = 400.0,
    SCANPOINTS         = 2000,
    ENERGY             = 2080.0,
    ASYMMETRY_ANGLE    = 0.0,
    THICKNESS          = 0.7,
    MOSAIC_FWHM        = 0.1,
    RSAG               = 125.0,
    RMER               = 1290.0,
    ANISOTROPY         = 0,
    POISSON            = 0.22,
    CUT                = "2 -1 -1 ; 1 1 1 ; 0 0 0",
    FILECOMPLIANCE     = "mycompliance.dat",
    )


# #
# # example plot
# #
# from srxraylib.plot.gol import plot
data2 = numpy.loadtxt("diff_pat.dat", skiprows=5)
# plot(data[:,0], data[:,-1])
#
# #
# # end script
# #


print("Calculation time: ", time.time() - t0)
import numpy





r1 = data1
r2 = data2

from oasys.util.oasys_util import get_fwhm

fwhm1, tmp, tmp = get_fwhm(r1[:,-1],r1[:,0])
fwhm2, tmp, tmp = get_fwhm(r2[:,-1],r2[:,0])
energy1=1385.6
energy2=2080.0
peak1=r1[:,-1].max()
peak2=r2[:,-1].max()
integrated_intensity1 = numpy.trapz(r1[:,-1],r1[:,0])
integrated_intensity2 = numpy.trapz(r2[:,-1],r2[:,0])


from srxraylib.plot.gol import plot
import matplotlib.pylab as plt
import matplotlib
font = {'family' : 'normal',
        'size'   : 22}

matplotlib.rc('font', **font)

p = plot(r1[:,0],r1[:,-1],r2[:,0],r2[:,-1],
         xtitle=r"$\theta-\theta_B [\mu rad]$",
         ytitle="Reflectivity",
         legend=["YB66 (004) E=%4.1f eV " %(energy1),
                 "YB66 (006) E=%4.1f eV " %(energy2) ],
         show=0,yrange=[0,0.5],figsize=(16,8))

plt.savefig("/tmp_14_days/srio/fig6.png")
plt.show()

