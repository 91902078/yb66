import os
import numpy

from dabax_util import get_dabax_file, get_f0_from_f0coeff, get_f0_coeffs_from_dabax_file

if __name__ == "__main__":
    from srxraylib.plot.gol import plot, set_qt

    set_qt()

    filename = "f0_InterTables.dat"

    coeffs_Yplus3 = get_f0_coeffs_from_dabax_file(entry_name="Y3+", filename=filename)
    coeffs_Y = get_f0_coeffs_from_dabax_file(entry_name="Y", filename=filename)
    coeffs_Kr = get_f0_coeffs_from_dabax_file(entry_name="Kr", filename=filename)

    print(coeffs_Yplus3, coeffs_Kr)


    ratio = numpy.linspace(0,3,1000)

    f0_Yplus3 = get_f0_from_f0coeff(coeffs_Yplus3, ratio)
    f0_Y = get_f0_from_f0coeff(coeffs_Y, ratio)
    f0_Kr = get_f0_from_f0coeff(coeffs_Kr, ratio)


    from srxraylib.plot.gol import plot
    plot(ratio, f0_Yplus3,
         ratio, f0_Y,
         ratio, f0_Kr,
         legend=["Y3+", "Y", "Z=39-3"],
         xtitle=r"q (sin $\theta$ / $\lambda$)", ytitle="f0 [electron units]",
         title=filename)

