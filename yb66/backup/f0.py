from orangecontrib.xoppy.util.xoppy_xraylib_util import f0_xop
import numpy

def get_f0_from_f0coeff(f0coeff, ratio):

    icentral = len(f0coeff) // 2
    F0 = f0coeff[icentral]
    for i in range(icentral):
        F0 += f0coeff[i] * numpy.exp(-1.0 * f0coeff[i + icentral + 1] * ratio ** 2)
    return F0


if __name__ == "__main__":
    # Z_B = 5
    # Z_Y = 39


    f0coeff_B = f0_xop(39)
    f0coeff_Y = f0_xop(39)
    f0coeff_Yion = f0_xop(36)

    #'B-.'
    f0coeff_BB = [1.493, 1.0472, 0.7776, 0.64929, 1.0233, 0.050981, 21.37, 65.436, 0.36215, 21.354, 1.1387]
    # 'Y3+'
    f0coeff_YY = [6.3697, 10.29, 4.3719, 5.9527, 4.3852, 4.6028, 1.28, 13.169, 0.41449, 1.967, 1.2664]


    ratio = numpy.linspace(0,8,1000)

    print(">>>>", len(f0coeff_B), len(f0coeff_BB))
    f0_B = get_f0_from_f0coeff(f0coeff_B, ratio)
    f0_BB = get_f0_from_f0coeff(f0coeff_BB, ratio)
    f0_Y = get_f0_from_f0coeff(f0coeff_Y, ratio)
    f0_YY = get_f0_from_f0coeff(f0coeff_YY, ratio)
    f0_Yion = get_f0_from_f0coeff(f0coeff_Yion, ratio)

    from srxraylib.plot.gol import plot
    plot(ratio, f0_B,
         ratio, f0_BB,
         ratio, f0_Y,
         ratio, f0_YY,
         ratio, f0_Yion,
         legend=["B", "B-.", "Y", "Y3+", "Z=39-3"],
         xtitle=r"q (sin $\theta$ / $\lambda$", ytitle="f0 [electron units]")





