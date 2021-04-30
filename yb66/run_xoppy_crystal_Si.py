import numpy
from crystal_util import bragg_calc2, crystal_fh2
from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_calc, crystal_fh
from run_diff_pat import run_diff_pat
from srxraylib.plot.gol import plot

if __name__ == "__main__":

    #
    # old code
    #
    dic1a = bragg_calc(descriptor="Si",hh=1,kk=1,ll=1,temper=1.0,emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra")
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

    dic2a = bragg_calc2(descriptor="Si",hh=1,kk=1,ll=1,temper=1.0,emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra")
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

    plot(a1[:, 0], a1[:, -1],
         a2[:, 0], a2[:, -1],
         legend=["OLD code","NEW code"])
