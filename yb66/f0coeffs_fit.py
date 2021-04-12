import numpy
from dabax_access_f0 import get_f0_coeffs_from_dabax_file, get_f0_from_f0coeff
from symbol_to_from_atomic_number import symbol_to_from_atomic_number
from scipy.optimize import curve_fit

"""
M. Sanchez del Rio srio@esrf.eu, X.J. YU, xiaojiang@nus.edu.sg
Interpolation of f0 coefficients for a fractional charged atom
"""

def func(q, a1, a2, a3, a4, a5, a6, a7, a8, a9):
    return get_f0_from_f0coeff([a1, a2, a3, a4, a5, a6, a7, a8, a9], q)

def get_f0_coeffs(atoms, list_Zatom):
    """
    Return a Dict {"B-0.0455": [f0 coefficients], ..., "Y+3":[f0 coefficients],...}
    """
    AtomicChargeList = {}
    #first row is atomic number, it is integer
    UniqueAtomicNumber = list(sorted(set(list_Zatom)))
    charge = [ atoms[i]['charge'] for i in range(len(atoms))]         
    for x in  UniqueAtomicNumber:
        AtomicChargeList[str(x)]= [] 
    for i,x in enumerate(list_Zatom):
        if charge[i] not in AtomicChargeList[str(int(x))]:
            AtomicChargeList[str(int(x))].append(charge[i])      #Charge value
    return crystal_get_f0_coeffs(AtomicChargeList.items())

def crystal_get_f0_coeffs(AtomicList):
    """
    Input: AtomicList, a list of tuple {(5,[-0.0455,]), (39,[3,])}, same atom allows with different charge
    Out:   A Dict {"B-0.0455": [f0 coefficients], ..., "Y+3":[f0 coefficients],...}
    """
    f0coeffs = {}
    searchChargeNameNeg = ['1-','2-','3-','4-','5-','6-','7-']
    searchChargeNamePos = ['1+','2+','3+','4+','5+','6+','7+']
    qq = numpy.linspace(0,2,1000)       #q = 0 to 2
    for x in AtomicList:
        n = int(x[0])       #atomic number
        sym = symbol_to_from_atomic_number(n)
        f0 = get_f0_coeffs_from_dabax_file(entry_name=sym)
        if len(f0) == 0:
                raise("cannot find f0 coefficients for '" + sym + "'")
        for charge in x[1]: #may have multiple valences for same atom, B1+, B2+, etc
            k = int(charge)
            f01 = []
            if charge < 0:
                if k == charge:  #integer charge
                    f01 = get_f0_coeffs_from_dabax_file(entry_name=sym + searchChargeNameNeg[abs(k)-1])
                if len(f01) == 0:
                    ff = []
                    for i,s in enumerate(searchChargeNameNeg):
                        f01 = get_f0_coeffs_from_dabax_file(entry_name=sym + s)
                        if len(f01) > 0:
                            ff.append((-i-1,f01))
                            if (i+1) > abs(k):      #already find one with valence higher than atom carried charge
                                break
                    if len(ff) > 0:
                        f01 = ff[-1]

            if len(f01) == 0 and 0 != charge:  #not get a f0 in negative charge direction
                ff = []
                for i,s in enumerate(searchChargeNamePos):  #try to find one with positive charge
                    f01 = get_f0_coeffs_from_dabax_file(entry_name=sym + s)
                    if len(f01) > 0:
                        ff.append((i+1,f01))
                        if (i+1) > abs(k) or charge < 0:
                            break
                if len(ff) > 0:
                    f01 = ff[-1]

            if charge == 0: #always no fit for neutral atom
                f0coeffs[sym] = f0
                continue
            #following for charged atom
            if len(f01) == 0:
                raise("No 2nd atom found for linear fit f0 coefficients")
            if charge == f01[0]: #if charged atom already listed, just get it, no fit
                f0coeffs[sym+f'%+g'%charge] = f01[1]
                continue

            #do fitting here
            f0_1 = get_f0_from_f0coeff(f0, qq)
            f0_2 = get_f0_from_f0coeff(f01[1], qq)
            f00  = f0_1 + charge / f01[0] * (f0_2 - f0_1)
            p0 = f0         #neutral f0 for p0
            #if 2nd atom with valence closer to charge, use it instead of neutral atom
            if abs(charge-f01[0]) < abs(charge):    
                p0 = f01[1]
            f00_fit, pcov_fit = curve_fit(func, qq, f00, p0=p0)
            f0coeffs[sym+f'%+g'%charge] = f00_fit
    return  f0coeffs

if __name__ == "__main__":
    from srxraylib.plot.gol import plot, set_qt
    #from scipy.optimize import curve_fit

    set_qt()

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
    #def func(q, a1, a2, a3, a4, a5, a6, a7, a8, a9):
    #    return get_f0_from_f0coeff([a1, a2, a3, a4, a5, a6, a7, a8, a9], q)


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
    print("fitted Bdot srio: ", popt_Bdot)


    #
    # add this block to f0_InterTables.dat
    #
    print("\n#S  5  B3+\n#N 9\n#L a1  a2  a3  a4  c  b1  b2  b3  b4\n"+"%g "*9 % (tuple(popt_B3plus)))

    #
    # test remote B3+
    #
    try:
        import os
        os.remove("f0_InterTables.dat")
    except:
        pass

    filename = "f0_InterTables.dat"
    coeffs_B3plus_remote = get_f0_coeffs_from_dabax_file(entry_name="B3+", filename=filename)
    coeff_Bdot = numpy.array([])
    plot(q, f0_B3plus,
         q, get_f0_from_f0coeff(popt_B3plus, q),
         xtitle=r"q (sin $\theta$ / $\lambda$)", ytitle="f0 [electron units]",
         legend=["B3plus original", "B3plus from remote f0_InterTables.dat"],
         title=filename)