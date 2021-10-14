from orangecontrib.xoppy.util.xoppy_xraylib_util import f0_calc
from xoppy_dabax_util import f0_calc_dabax

import matplotlib.pylab as plt
from srxraylib.plot.gol import plot

from dabax_util import crystal_parser


import numpy

if __name__ == "__main__":
    dabax_repository = "/scisoft/DABAX/data"  # "http://ftp.esrf.fr/pub/scisoft/DabaxFiles/"
    do_plot = 1

    if True:
        #
        # f0
        #
        Y0_xrl = f0_calc      (0, "Y", 0, 3, 500)
        Y0_dbx = f0_calc_dabax(0, "Y", 0, 3, 500)
        Y3_dbx = f0_calc_dabax(0, "Y", 0, 3, 500, charge=3.0)

        B0_xrl = f0_calc      (0, "B", 0, 3, 500)
        B0_dbx = f0_calc_dabax(0, "B", 0, 3, 500)
        Bf_dbx = f0_calc_dabax(0, "B", 0, 3, 500, charge=-0.045)

        if do_plot:
            plot(Y0_xrl["data"][0, :], Y0_xrl["data"][1,:],
                 Y0_dbx["data"][0, :], Y0_dbx["data"][1,:],
                 Y3_dbx["data"][0, :], Y3_dbx["data"][1, :],
                 B0_xrl["data"][0, :], B0_xrl["data"][1, :],
                 B0_dbx["data"][0, :], B0_dbx["data"][1, :],
                 Bf_dbx["data"][0, :], Bf_dbx["data"][1, :],
                 #
                 # linestyle=[None,None,None],
                 # marker=[None,'+',None,'+'],
                 # color=['r','r','b','b'],
                 legend=['Y xraylib',r'Y$^{+0}$ dabax',r'Y$^{+3}$ dabax',
                         'Y xraylib',r'B$^{+0}$ dabax',r'B$^{-0.045}$ dabax'],
                 xtitle=r'q=sin$\theta$/$\lambda$',ytitle='f$_0$ [electron units]',
                 show=0)

            plt.savefig("f0.png")
            print("File f0.png written to disk.")
            plt.show()

    #
    # f0 another test
    #
    if False:

        #
        # test f0 data for B3+
        #
        q = numpy.array(
            [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
             1.7, 1.8, 1.9])
        f0_B3plus = numpy.array(
            [2, 1.995, 1.979, 1.954, 1.919, 1.875, 1.824, 1.766, 1.703, 1.566, 1.42, 1.274, 1.132, 0.999, 0.877, 0.767,
             0.669, 0.582, 0.507, 0.441, 0.384, 0.335, 0.293, 0.256])

        #
        # plot
        #
        from srxraylib.plot.gol import plot

        coeff_Bdot = numpy.array([])
        plot(q, f0_B3plus,
             q, calculate_f0_from_f0coeff(f0_with_fractional_charge(5, 3.0, dabax_repository=dabax_repository), q),
             q, calculate_f0_from_f0coeff(f0_with_fractional_charge(5, 2.8, dabax_repository=dabax_repository), q),
             xtitle=r"q (sin $\theta$ / $\lambda$)", ytitle="f0 [electron units]",
             legend=["B3plus original",
                     "B3plus from f0_with_fractional_charge(5,+3)",
                     "B3plus from f0_with_fractional_charge(5,+2.8)"],
             title="", show=1)







