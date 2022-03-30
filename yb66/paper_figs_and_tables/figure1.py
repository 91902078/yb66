import matplotlib.pylab as plt
from srxraylib.plot.gol import plot, set_qt
from xoppylib.decorators.xraylib_decorated import XraylibDecorated
from xoppylib.decorators.dabax_decorated import DabaxDecorated


if __name__ == "__main__":

    set_qt()

    xrl = XraylibDecorated()
    dx = DabaxDecorated()

    print(dx.info())

    do_plot = 1

    #
    # f0
    #
    Y0_xrl = xrl.f0_calc(0, "Y", 0, 3, 50)
    Y0_dbx =  dx.f0_calc(0, "Y", 0, 3, 50)
    Y3_dbx =  dx.f0_calc(0, "Y", 0, 3, 50, charge=3.0)
    B0_xrl = xrl.f0_calc(0, "B", 0, 3, 50)
    B0_dbx =  dx.f0_calc(0, "B", 0, 3, 50)
    Bf_dbx =  dx.f0_calc(0, "B", 0, 3, 50, charge=-0.045)

    if do_plot:
        plot(Y0_xrl["data"][0, :], Y0_xrl["data"][1,:],
             Y0_dbx["data"][0, :], Y0_dbx["data"][1,:],
             Y3_dbx["data"][0, :], Y3_dbx["data"][1, :],
             B0_xrl["data"][0, :], B0_xrl["data"][1, :],
             B0_dbx["data"][0, :], B0_dbx["data"][1, :],
             Bf_dbx["data"][0, :], Bf_dbx["data"][1, :],
             legend=['Y xraylib',r'Y$^{+0}$ dabax',r'Y$^{+3}$ dabax',
                     'Y xraylib',r'B$^{+0}$ dabax',r'B$^{-0.045}$ dabax'],
             xtitle=r'q=sin$\theta$/$\lambda$',ytitle='f$_0$ [electron units]',
             linestyle=["-.",None,"dashed","-.",None,""],
             marker=[".","None","None",'+',"None","o"],
             show=0)

        filename = "/tmp_14_days/srio/fig1.png"
        plt.savefig(filename)
        print("File %s written to disk." % filename)
        plt.show()

