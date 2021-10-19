from orangecontrib.xoppy.util.xoppy_xraylib_util import f0_calc
from xoppy_dabax_util import f0_calc_dabax

import matplotlib.pylab as plt
from srxraylib.plot.gol import plot

from dabax_util import crystal_parser


import numpy

if __name__ == "__main__":
    dabax_repository = "/scisoft/DABAX/data"  # "http://ftp.esrf.fr/pub/scisoft/DabaxFiles/"
    do_plot = 1

    if False:
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
# table 3
#
    if True:
        from xoppy_dabax_util import check_structure_factor

        # test Muscovite
        F200 = check_structure_factor(descriptor="YB66", hh=2, kk=0, ll=0, energy=8040.0, do_assert=0, models=[0,0,1])
        F400 = check_structure_factor(descriptor="YB66", hh=4, kk=0, ll=0, energy=8040.0, do_assert=0, models=[0,0,1])
        F600 = check_structure_factor(descriptor="YB66", hh=6, kk=0, ll=0, energy=8040.0, do_assert=0, models=[0,0,1])
        F800 = check_structure_factor(descriptor="YB66", hh=8, kk=0, ll=0, energy=8040.0, do_assert=0, models=[0,0,1])


        print("F_200: ",numpy.abs(F200), numpy.abs(F200) / 4)
        print("F_400: ",numpy.abs(F400), numpy.abs(F400) / 4)
        print("F_600: ",numpy.abs(F600), numpy.abs(F600) / 4)
        print("F_800: ",numpy.abs(F800), numpy.abs(F800) / 4)