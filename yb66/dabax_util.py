import os
import numpy
from urllib.request import urlretrieve
from silx.io.specfile import SpecFile

from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_metrictensor
import scipy.constants as codata

from scipy.optimize import curve_fit

"""
X.J. YU, xiaojiang@nus.edu.sg, M. Sanchez del Rio srio@esrf.eu
"""

# global dabax_repository
dabax_repository="http://ftp.esrf.eu/pub/scisoft/DabaxFiles/"


#
# common access tools
#
def get_dabax_file(filename, dabax_repository=dabax_repository):

    #
    # file exists in current directory
    #
    if os.path.exists(filename):
        print("Dabax file exists in local directory: %s " % filename)
        return filename

    #
    # download remote file
    #
    if dabax_repository[0:3] == "htt" or dabax_repository[0:3] == "ftp":
        try:
            filepath, http_msg = urlretrieve(dabax_repository + filename,
                        filename=filename,
                        reporthook=None,
                        data=None)
            print("Dabax file %s downloaded from %s" % (filepath, dabax_repository + filename))
            return filename
        except:
            raise Exception("Failed to download file %s from %s" % (filename, dabax_repository))

    #
    # file exists in local repository
    #
    f1 = os.path.join(dabax_repository,filename)
    if os.path.exists(f1):
        print("Dabax file exists in local directory: %s " % f1)
        return f1

    raise Exception(FileNotFoundError)


#
# crystal
#
def crystal_parser(filename='Crystals.dat', entry_name='YB66',
                   dabax_repository=dabax_repository):
    """
    parse a complex crystal structure file into a dictionary (like xraylib.Crystal_GetCrystal('Si'))

    it has an additional fiels for each atom: the charge

    return a dictionary containing crystal infomation
    """

    file1 = get_dabax_file(filename, dabax_repository=dabax_repository)

    sf = SpecFile(file1)

    flag_found = False

    for index in range(len(sf)):
        s1 = sf[index]
        name = s1.scan_header_dict["S"]
        if name.split(' ')[1] == entry_name:
            flag_found = True
            index_found = index

    if not flag_found:
        raise (Exception("Entry name not found: %s" % entry_name))



    cryst = {'name':entry_name}     #returned dictionary like that one created by xraylib.Crystal_GetCrystal(descriptor)

    cell_parameters = sf[index_found].scan_header_dict["UCELL"]
    cell_parameters = ' '.join(cell_parameters.split()) # remove multiple blanks


    a = cell_parameters.split(' ')
    cryst['a'] =     float(a[0])
    cryst['b'] =     float(a[1])
    cryst['c'] =     float(a[2])
    cryst['alpha'] = float(a[3])
    cryst['beta'] =  float(a[4])
    cryst['gamma'] = float(a[5])

    volume = bragg_metrictensor(float(a[0]), float(a[1]), float(a[2]),
                                float(a[3]), float(a[4]), float(a[5]), RETURN_VOLUME=1)

    cryst['volume'] = volume

    cell_data = numpy.array(sf[index_found].data)

    cryst['n_atom'] = cell_data.shape[1]
    atom = []

    for i in range(cell_data.shape[1]):
        if cell_data.shape[0] == 5: # standard 5 columns
            # not here, this info is not in the dabax file
            # s = symbol_to_from_atomic_number(int(cell_data[0,i]))
            atom.append({
                         # 'AtomicName': s,
                         'Zatom':int(cell_data[0,i]),
                         'fraction':cell_data[1,i],
                         'x': cell_data[2,i],
                         'y': cell_data[3, i],
                         'z': cell_data[4, i],
                         'charge': 0.0,})
        else: # 6 columns (charge)
            #'AtomicName' required to compatible my current code
            # s = symbol_to_from_atomic_number(int(cell_data[0,i]))
            # if cell_data[5, i] != 0:  #charged
            #     s = s + f'%+.6g'%cell_data[5, i]
            atom.append({
                         # 'AtomicName': s,
                         'Zatom':int(cell_data[0,i]),
                         'fraction':cell_data[1,i],
                         'x': cell_data[2,i],
                         'y': cell_data[3, i],
                         'z': cell_data[4, i],
                         'charge': cell_data[5, i],})

    cryst['atom'] = atom
    cryst['cpointer'] = None

    ANISO_KEY = "UANISO_COFF"   #prefix for a line with anisotropic coefficients
    d = sf[index_found].scan_header_dict
    AnisoItem = {'Name': '       ',
                 'start': 0,
                 'end': 0,
                 'beta11': 0.0,
                 'beta22': 0.0,
                 'beta33': 0.0,
                 'beta12': 0.0,
                 'beta13': 0.0,
                 'beta23': 0.0}

    a=[ (x, d[x].split()) for x in d if x[:len(ANISO_KEY)] == ANISO_KEY]
    if len(a) >0:       #found Anisotropic coefficients in the header, process it
        a=sorted(a,key=lambda x:int(x[1][0]),reverse=False)     #sort 'Start' ascendant, avoid order changed by the SpecFile
        n = 0
        Aniso = []
        for x in a: #tuple('UANISO_COFF_B1',[1 96 0.00038 0.00044 0.00039 0 0 0])
            AnisoItem['Name']=   x[0][len(ANISO_KEY)+1:]      #get site atom name starting from 13th character 'B1', etc
            AnisoItem['start']=  int(x[1][0])
            AnisoItem['end']=    int(x[1][1])
            AnisoItem['beta11']= float(x[1][2])
            AnisoItem['beta22']= float(x[1][3])
            AnisoItem['beta33']= float(x[1][4])
            AnisoItem['beta12']= float(x[1][5])
            AnisoItem['beta13']= float(x[1][6])
            AnisoItem['beta23']= float(x[1][7])
            Aniso.append(AnisoItem.copy())
            n = n + 1
        cryst['Aniso'] = Aniso      #if having key 'Ansio' when there is anisotropic data,otherwise no
        cryst['n_aniso']= n
    else:       #create a dummy Aniso to compatible with xraylib
        cryst['Aniso'] = [AnisoItem.copy()]
        cryst['n_aniso']= 1
 
    return cryst

def Crystal_GetCrystalsList(dabax_repository=dabax_repository):
    """
    get crystal names from crystals.dat
    """
    file1 = get_dabax_file('Crystals.dat', dabax_repository=dabax_repository)
    sf = SpecFile(file1)
    crystals = []
    for index in range(len(sf)):
        s1 = sf[index]
        name = s1.scan_header_dict["S"]
        crystals.append(name.split(' ')[1])
    
    return crystals

#
# dabax crystal functions with the same interface as xraylib
#

#
#
#   TODO:
#          F_0 = xraylib.Crystal_F_H_StructureFactor(_crystal, E_keV, h, k, l, _debyeWaller, 1.0)
#
#          F_H = xraylib.Crystal_F_H_StructureFactor(_crystal, E_keV, h, k, l, _debyeWaller, 1.0)
#

def Crystal_GetCrystal(descriptor, dabax_repository=dabax_repository):
    return crystal_parser(filename='Crystals.dat', entry_name=descriptor, dabax_repository=dabax_repository)

def Crystal_dSpacing(cryst, h, k, l):
    return bragg_metrictensor(cryst['a'], cryst['b'], cryst['c'],
                              cryst['alpha'], cryst['beta'], cryst['gamma'],
                              HKL=[h, k, l])

def Bragg_angle(cryst, E_keV, h, k, l):
    dspacing = Crystal_dSpacing(cryst, h, k, l)  # in A
    wavelength = codata.h * codata.c / codata.e / (E_keV * 1e3) * 1e10 # in A
    return numpy.arcsin(wavelength / 2 / dspacing)



#
# f0
#

def get_f0_coeffs_from_dabax_file(entry_name="Y3+", filename="f0_InterTables.dat",
                                  dabax_repository=dabax_repository):

    # if getattr(get_f0_coeffs_from_dabax_file,'sf') is not None:
    #     sf = getattr(get_f0_coeffs_from_dabax_file,'sf')
    # else:
    file1 = get_dabax_file(filename, dabax_repository=dabax_repository)
    sf = SpecFile(file1)
        # get_f0_coeffs_from_dabax_file.sf = sf

    flag_found = False

    for index in range(len(sf)):
        s1 = sf[index]
        name = s1.scan_header_dict["S"]

        if name.split('  ')[1] == entry_name:
            flag_found = True
            index_found = index

    if flag_found:
        return (sf[index_found].data)[:,0]
    else:
        return []
#setattr(get_f0_coeffs_from_dabax_file,'sf',None)


def f0_with_fractional_charge(Z, charge=0.0, filename="f0_InterTables.dat", dabax_repository=dabax_repository):
    symbol = __symbol_to_from_atomic_number(Z)

    if charge == 0.0:
        return get_f0_coeffs_from_dabax_file(entry_name=symbol, filename=filename, dabax_repository=dabax_repository)
    else:
        # retrieve all entries

        filename = get_dabax_file(filename, dabax_repository=dabax_repository)
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
        popt, pcov = curve_fit(__f0func, q, f0, p0=coefficient_list[sorted_index_0], maxfev=20000)
        print("fitted: ", popt)

        return popt
    except:
        print("Error: failed to fit coefficients for fractional charge. Returning the ones of ",
              interesting_entries[sorted_index_0])
        return coefficient_list[sorted_index_0]


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

def __f0func(q, a1, a2, a3, a4, a5, a6, a7, a8, a9):
    return calculate_f0_from_f0coeff([a1, a2, a3, a4, a5, a6, a7, a8, a9], q)



if __name__ == "__main__":
    # dabax_repository = "http://ftp.esrf.eu/pub/scisoft/DabaxFiles/"
    dabax_repository = "/scisoft/DABAX/data" # "http://ftp.esrf.fr/pub/scisoft/DabaxFiles/"

    #
    # crystal tests
    #
    if True:
        print(get_dabax_file("Crystals.dat", dabax_repository=dabax_repository))

        print(get_f0_coeffs_from_dabax_file(entry_name="Y3+",
                                            filename="f0_InterTables.dat",
                                            dabax_repository=dabax_repository))

        print(Crystal_GetCrystalsList(dabax_repository=dabax_repository))

        yb = crystal_parser(filename='Crystals.dat', entry_name='YB66', dabax_repository=dabax_repository)

        si = Crystal_GetCrystal("Si", dabax_repository=dabax_repository)
        print("Si 111 d-spacing: ", Crystal_dSpacing(si,1,1,1))
        print("Si 111 bragg angle at 10 keV [deg]: ", 180 / numpy.pi * Bragg_angle(si,10, 1,1,1))

    #
    # crystal vs xraylib tests
    #

    if True:
        print(crystal_parser(filename='Crystals.dat', entry_name='YB66', dabax_repository=dabax_repository))

        # compare with xraylib
        xdabax = crystal_parser(filename='Crystals.dat', entry_name='Si', dabax_repository=dabax_repository)

        import xraylib
        xxraylib = xraylib.Crystal_GetCrystal('Si')

        for key in xxraylib.keys():
            tmp = xxraylib[key]

            if isinstance(tmp, list):
                for i, element in enumerate(tmp):
                    print(key, i, xdabax[key][i], xxraylib[key][i])
            else:
                print(key, xdabax[key], xxraylib[key])

    #
    # f0
    #
    if True:
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
             q, calculate_f0_from_f0coeff(f0_with_fractional_charge(5, 3.0, dabax_repository=dabax_repository), q),
             q, calculate_f0_from_f0coeff(f0_with_fractional_charge(5, 2.8, dabax_repository=dabax_repository), q),
             xtitle=r"q (sin $\theta$ / $\lambda$)", ytitle="f0 [electron units]",
             legend=["B3plus original",
                     "B3plus from f0_with_fractional_charge(5,+3)",
                     "B3plus from f0_with_fractional_charge(5,+2.8)"],
             title="")

