import numpy
from dabax_access_f0 import get_f0_coeffs_from_dabax_file, get_f0_from_f0coeff


if __name__ == "__main__":
    from srxraylib.plot.gol import plot
    from scipy.optimize import curve_fit

    filename = "f0_InterTables.dat"
    coeffs_B = get_f0_coeffs_from_dabax_file(entry_name="B", filename=filename)


    #
    # f0 data
    #
    q = numpy.array([0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9])

    f0_B = get_f0_from_f0coeff(coeffs_B, q)
    f0_B3plus = numpy.array([2,1.995,1.979,1.954,1.919,1.875,1.824,1.766,1.703,1.566,1.42,1.274,1.132,0.999,0.877,0.767,0.669,0.582,0.507,0.441,0.384,0.335,0.293,0.256])
    f0_Bdot = f0_B + (-0.0455) / 3 * (f0_B3plus - f0_B)

    #
    # fit
    #
    def func(q, a1, a2, a3, a4, a5, a6, a7, a8, a9):
        return get_f0_from_f0coeff([a1, a2, a3, a4, a5, a6, a7, a8, a9], q)


    popt_B3plus, pcov_B3plus = curve_fit(func, q, f0_B3plus, p0=coeffs_B)
    print("fitted B3+: ", popt_B3plus)

    popt_Bdot, pcov_Bdot = curve_fit(func, q, f0_Bdot, p0=coeffs_B)
    print("fitted Bdot: ", popt_Bdot)



    #
    # plot
    #
    from srxraylib.plot.gol import plot

    coeff_Bdot = numpy.array([])
    plot(q, f0_B3plus,
         q, get_f0_from_f0coeff(popt_B3plus, q),
         xtitle=r"q (sin $\theta$ / $\lambda$)", ytitle="f0 [electron units]",
         legend=["B3plus original", "B3plus from srio fit"],
         title=filename)


    coeff_Bdot = numpy.array([0.858,0.89669,1.0756,2.118,0.095903,0.46461,1.2126,61.273,23.55])
    plot(q, f0_Bdot,
         q, get_f0_from_f0coeff(coeff_Bdot, q),
         q, get_f0_from_f0coeff(popt_Bdot, q),
         xtitle=r"q (sin $\theta$ / $\lambda$)", ytitle="f0 [electron units]",
         legend=["Bdot original", "Bdot from Xiaojiang fit","Bdot from srio fit",],
         title=filename)

    print("fitted Bdot Xiaojiang: ", coeff_Bdot)
    print("fitted srio Xiaojiang: ", popt_Bdot)