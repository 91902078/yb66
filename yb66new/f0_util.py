import numpy
from dabax_util import get_f0_coeffs_from_dabax_file, get_dabax_file
from scipy.optimize import curve_fit
from silx.io.specfile import SpecFile


def f0_with_fractional_charge(Z, charge=0.0, filename="f0_InterTables.dat"):
    symbol = __symbol_to_from_atomic_number(Z)

    if charge == 0.0:
        return get_f0_coeffs_from_dabax_file(entry_name=symbol, filename=filename)
    else:
        # retrieve all entries

        error_flag = get_dabax_file(filename)
        if error_flag == False:
            raise (FileNotFoundError)
        sf = SpecFile(filename)

        # retrieve all entries
        entries = []
        for index in range(len(sf)):
            s1 = sf[index]
            name = s1.scan_header_dict["S"]
            entries.append(name.split('  ')[1])
        print(entries)

        # identify the entries that match the symbol
        interesting_entries = []
        charge_list = []
        index_list = []
        for i,entry in enumerate(entries):
            if entry.find(symbol) == 0:
                if entry == symbol:
                    interesting_entries.append(entry)
                    charge_list.append(0.0)
                    index_list.append(i)
                else:
                    entry2 = entry.replace(symbol,'')
                    try:
                        charge_item = int(entry2[::-1]) # convert to integer the reversed string
                        charge_list.append(charge_item)
                        interesting_entries.append(entry)
                        index_list.append(i)
                    except:
                        pass

        # retrieve coefficients from these interesting entries
        coefficient_list = []
        for i in index_list:
            coefficient_list.append((sf[i].data)[:, 0])

        return __f0_interpolate_coefficients(charge, interesting_entries, charge_list, coefficient_list)


def calculate_f0_from_f0coeff(f0coeff, ratio):

    icentral = len(f0coeff) // 2
    F0 = f0coeff[icentral]
    for i in range(icentral):
        F0 += f0coeff[i] * numpy.exp(-1.0 * f0coeff[i + icentral + 1] * ratio ** 2)
    return F0

#
# __* are auxiliar routines not to be exported outside
#

def __f0_interpolate_coefficients(charge, interesting_entries, charge_list, coefficient_list):
    #
    # f0 data
    #

    nitems = len(interesting_entries)

    if nitems == 1:
        print("Warning: no interpolating of charge: only one value available for ", interesting_entries[0])
        return coefficient_list[0]

    charge_list_difference = []
    for i in range(nitems):
        charge_list_difference.append(charge_list[i] - charge)

    charge_list_difference = numpy.array(charge_list_difference)  # convert to numpy array

    if numpy.abs(charge_list_difference).min() == 0:
        idx = numpy.abs(charge_list_difference).argmin()
        print("No interpolating needed: returning value for ", interesting_entries[idx])
        return coefficient_list[idx]

    # get the closer charge values, no matter of the sign

    sorted_indices = numpy.argsort(numpy.abs(charge_list_difference))
    sorted_index_0 = sorted_indices[0]
    sorted_index_1 = sorted_indices[1]
    delta_data = charge_list[sorted_index_1] - charge_list[sorted_index_0]
    delta_charge = charge - charge_list[sorted_index_0]
    delta = delta_charge / delta_data
    print("Interpolating charge %g = %s + %g (%s - %s)" % (charge,
                                                           interesting_entries[sorted_index_0],
                                                           delta,
                                                           interesting_entries[sorted_index_1],
                                                           interesting_entries[sorted_index_0]))

    # interpolate to get the f0 for the wanted charge

    q = numpy.linspace(0.0, 2.0, 100)
    f0_0 = calculate_f0_from_f0coeff(coefficient_list[sorted_index_0], q)
    f0_1 = calculate_f0_from_f0coeff(coefficient_list[sorted_index_1], q)
    f0 = f0_0 + delta * (f0_1 - f0_0)

    #
    # fit
    #
    try:
        popt, pcov = curve_fit(__func, q, f0, p0=coefficient_list[sorted_index_0], maxfev=20000)
        print("fitted: ", popt)

        return popt
    except:
        print("Error: failed to fit coefficients for fractional charge. Returning the ones of ",
              interesting_entries[sorted_index_0])
        return coefficient_list[sorted_index_0]



#
# tools
#

def __symbol_to_from_atomic_number(ATOM):
    atoms = [ [1 ,"H"] ,[2 ,"He"] ,[3 ,"Li"] ,[4 ,"Be"] ,[5 ,"B"] ,[6 ,"C"] ,[7 ,"N"] ,[8 ,"O"] ,[9 ,"F"] ,[10 ,"Ne"], \
                 [11 ,"Na"] ,[12 ,"Mg"] ,[13 ,"Al"] ,[14 ,"Si"] ,[15 ,"P"] ,[16 ,"S"] ,[17 ,"Cl"] ,[18 ,"Ar"] ,[19 ,"K"]
             ,[20 ,"Ca"], \
                 [21 ,"Sc"] ,[22 ,"Ti"] ,[23 ,"V"] ,[24 ,"Cr"] ,[25 ,"Mn"] ,[26 ,"Fe"] ,[27 ,"Co"] ,[28 ,"Ni"] ,[29 ,"Cu"]
             ,[30 ,"Zn"], \
                 [31 ,"Ga"] ,[32 ,"Ge"] ,[33 ,"As"] ,[34 ,"Se"] ,[35 ,"Br"] ,[36 ,"Kr"] ,[37 ,"Rb"] ,[38 ,"Sr"] ,[39 ,"Y"]
             ,[40 ,"Zr"], \
                 [41 ,"Nb"] ,[42 ,"Mo"] ,[43 ,"Tc"] ,[44 ,"Ru"] ,[45 ,"Rh"] ,[46 ,"Pd"] ,[47 ,"Ag"] ,[48 ,"Cd"] ,[49 ,"In"]
             ,[50 ,"Sn"], \
                 [51 ,"Sb"] ,[52 ,"Te"] ,[53 ,"I"] ,[54 ,"Xe"] ,[55 ,"Cs"] ,[56 ,"Ba"] ,[57 ,"La"] ,[58 ,"Ce"] ,[59 ,"Pr"]
             ,[60 ,"Nd"], \
                 [61 ,"Pm"] ,[62 ,"Sm"] ,[63 ,"Eu"] ,[64 ,"Gd"] ,[65 ,"Tb"] ,[66 ,"Dy"] ,[67 ,"Ho"] ,[68 ,"Er"] ,[69 ,"Tm"]
             ,[70 ,"Yb"], \
                 [71 ,"Lu"] ,[72 ,"Hf"] ,[73 ,"Ta"] ,[74 ,"W"] ,[75 ,"Re"] ,[76 ,"Os"] ,[77 ,"Ir"] ,[78 ,"Pt"] ,[79 ,"Au"]
             ,[80 ,"Hg"], \
                 [81 ,"Tl"] ,[82 ,"Pb"] ,[83 ,"Bi"] ,[84 ,"Po"] ,[85 ,"At"] ,[86 ,"Rn"] ,[87 ,"Fr"] ,[88 ,"Ra"] ,[89 ,"Ac"]
             ,[90 ,"Th"], \
                 [91 ,"Pa"] ,[92 ,"U"] ,[93 ,"Np"] ,[94 ,"Pu"] ,[95 ,"Am"] ,[96 ,"Cm"] ,[97 ,"Bk"] ,[98 ,"Cf"] ,[99 ,"Es"]
             ,[100 ,"Fm"], \
                 [101 ,"Md"] ,[102 ,"No"] ,[103 ,"Lr"] ,[104 ,"Rf"] ,[105 ,"Db"] ,[106 ,"Sg"] ,[107 ,"Bh"] 	]

    if isinstance(ATOM ,int):
        for a in atoms:
            if a[0] == ATOM:
                return a[1]
    for a in atoms:
        if a[1] == ATOM:
            return int(a[0])

    raise Exception("Why are you here?")

def __func(q, a1, a2, a3, a4, a5, a6, a7, a8, a9):
    return calculate_f0_from_f0coeff([a1, a2, a3, a4, a5, a6, a7, a8, a9], q)

#
# TODO: delete after switching to f0_with_fractional_charge()
#
# def __get_f0_coeffs(atoms, list_Zatom):
#     """
#     Return a Dict {"B-0.0455": [f0 coefficients], ..., "Y+3":[f0 coefficients],...}
#     """
#     AtomicChargeList = {}
#     #first row is atomic number, it is integer
#     UniqueAtomicNumber = list(sorted(set(list_Zatom)))
#     charge = [ atoms[i]['charge'] for i in range(len(atoms))]
#     for x in  UniqueAtomicNumber:
#         AtomicChargeList[str(x)]= []
#     for i,x in enumerate(list_Zatom):
#         if charge[i] not in AtomicChargeList[str(int(x))]:
#             AtomicChargeList[str(int(x))].append(charge[i])      #Charge value
#     return __crystal_get_f0_coeffs(AtomicChargeList.items())
#
# def __crystal_get_f0_coeffs(AtomicList):
#     """
#     Input: AtomicList, a list of tuple {(5,[-0.0455,]), (39,[3,])}, same atom allows with different charge
#     Out:   A Dict {"B-0.0455": [f0 coefficients], ..., "Y+3":[f0 coefficients],...}
#     """
#     f0coeffs = {}
#     searchChargeNameNeg = ['1-','2-','3-','4-','5-','6-','7-']
#     searchChargeNamePos = ['1+','2+','3+','4+','5+','6+','7+']
#     qq = numpy.linspace(0,2,1000)       #q = 0 to 2
#     for x in AtomicList:
#         n = int(x[0])       #atomic number
#         sym = __symbol_to_from_atomic_number(n)
#         f0 = get_f0_coeffs_from_dabax_file(entry_name=sym)
#         if len(f0) == 0:
#                 raise("cannot find f0 coefficients for '" + sym + "'")
#         for charge in x[1]: #may have multiple valences for same atom, B1+, B2+, etc
#             k = int(charge)
#             f01 = []
#             if charge < 0:
#                 if k == charge:  #integer charge
#                     f01 = get_f0_coeffs_from_dabax_file(entry_name=sym + searchChargeNameNeg[abs(k)-1])
#                 if len(f01) == 0:
#                     ff = []
#                     for i,s in enumerate(searchChargeNameNeg):
#                         f01 = get_f0_coeffs_from_dabax_file(entry_name=sym + s)
#                         if len(f01) > 0:
#                             ff.append((-i-1,f01))
#                             if (i+1) > abs(k):      #already find one with valence higher than atom carried charge
#                                 break
#                     if len(ff) > 0:
#                         f01 = ff[-1]
#
#             if len(f01) == 0 and 0 != charge:  #not get a f0 in negative charge direction
#                 ff = []
#                 for i,s in enumerate(searchChargeNamePos):  #try to find one with positive charge
#                     f01 = get_f0_coeffs_from_dabax_file(entry_name=sym + s)
#                     if len(f01) > 0:
#                         ff.append((i+1,f01))
#                         if (i+1) > abs(k) or charge < 0:
#                             break
#                 if len(ff) > 0:
#                     f01 = ff[-1]
#
#             if charge == 0: #always no fit for neutral atom
#                 f0coeffs[sym] = f0
#                 continue
#             #following for charged atom
#             if len(f01) == 0:
#                 raise("No 2nd atom found for linear fit f0 coefficients")
#             if charge == f01[0]: #if charged atom already listed, just get it, no fit
#                 f0coeffs[sym+f'%+g'%charge] = f01[1]
#                 continue
#
#             #do fitting here
#             f0_1 = calculate_f0_from_f0coeff(f0, qq)
#             f0_2 = calculate_f0_from_f0coeff(f01[1], qq)
#             f00  = f0_1 + charge / f01[0] * (f0_2 - f0_1)
#             p0 = f0         #neutral f0 for p0
#             #if 2nd atom with valence closer to charge, use it instead of neutral atom
#             if abs(charge-f01[0]) < abs(charge):
#                 p0 = f01[1]
#             f00_fit, pcov_fit = curve_fit(__func, qq, f00, p0=p0)
#             f0coeffs[sym+f'%+g'%charge] = f00_fit
#     return  f0coeffs


if __name__ == "__main__":
    # print(f0_with_fractional_charge(39,0.0) )
    # print(f0_with_fractional_charge(39, 3.5) )
    # print(f0_with_fractional_charge(39, 2.5))
    # print(f0_with_fractional_charge(39, 0.1))
    # print(f0_with_fractional_charge(39, -0.1))


    #
    # test f0 data for B3+
    #
    q = numpy.array([0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9])
    f0_B3plus = numpy.array([2,1.995,1.979,1.954,1.919,1.875,1.824,1.766,1.703,1.566,1.42,1.274,1.132,0.999,0.877,0.767,0.669,0.582,0.507,0.441,0.384,0.335,0.293,0.256])

    #
    # plot
    #
    from srxraylib.plot.gol import plot

    coeff_Bdot = numpy.array([])
    plot(q, f0_B3plus,
         q, calculate_f0_from_f0coeff(f0_with_fractional_charge(5, 3.0), q),
         q, calculate_f0_from_f0coeff(f0_with_fractional_charge(5, 2.8), q),
         xtitle=r"q (sin $\theta$ / $\lambda$)", ytitle="f0 [electron units]",
         legend=["B3plus original",
                 "B3plus from f0_with_fractional_charge(5,+3)",
                 "B3plus from f0_with_fractional_charge(5,+2.8)"],
         title="")




