import os
import numpy
from urllib.request import urlretrieve
from silx.io.specfile import SpecFile

from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_metrictensor, interface_reflectivity
from xoppy_xraylib_util import f1f2_calc
import scipy.constants as codata

from scipy.optimize import curve_fit

import scipy.constants as codata

toangstroms = codata.h * codata.c / codata.e * 1e10

import xraylib # TODO: remove this dependency, used for density etc.

"""
X.J. YU, xiaojiang@nus.edu.sg, M. Sanchez del Rio srio@esrf.eu
"""

# global dabax_repository
dabax_repository="http://ftp.esrf.eu/pub/scisoft/DabaxFiles/"


#
# common access tools
#
def get_dabax_file(filename, dabax_repository=dabax_repository, verbose=True):

    #
    # file exists in current directory
    #
    if os.path.exists(filename):
        if verbose: print("Dabax file exists in local directory: %s " % filename)
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
            if verbose: print("Dabax file %s downloaded from %s" % (filepath, dabax_repository + filename))
            return filename
        except:
            raise Exception("Failed to download file %s from %s" % (filename, dabax_repository))

    #
    # file exists in local repository
    #
    f1 = os.path.join(dabax_repository,filename)
    if os.path.exists(f1):
        if verbose: print("Dabax file exists in local directory: %s " % f1)
        return f1

    raise Exception(FileNotFoundError)


#
# crystal
#
def crystal_parser(filename='Crystals.dat', entry_name='YB66',
                   dabax_repository=dabax_repository, verbose=True):
    """
    parse a complex crystal structure file into a dictionary (like xraylib.Crystal_GetCrystal('Si'))

    it has an additional fiels for each atom: the charge

    return a dictionary containing crystal infomation
    """

    file1 = get_dabax_file(filename, dabax_repository=dabax_repository, verbose=verbose)

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


def f0_with_fractional_charge(Z, charge=0.0, filename="f0_InterTables.dat", dabax_repository=dabax_repository,
                              verbose=True):
    symbol = __symbol_to_from_atomic_number(Z)

    if charge == 0.0:
        return get_f0_coeffs_from_dabax_file(entry_name=symbol, filename=filename, dabax_repository=dabax_repository)
    else:
        # retrieve all entries

        filename = get_dabax_file(filename, dabax_repository=dabax_repository, verbose=verbose)
        sf = SpecFile(filename)

        # retrieve all entries
        entries = []
        for index in range(len(sf)):
            s1 = sf[index]
            name = s1.scan_header_dict["S"]
            entries.append(name.split('  ')[1])
        # print(entries)

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

        return __f0_interpolate_coefficients(charge, interesting_entries, charge_list, coefficient_list, verbose=verbose)


def calculate_f0_from_f0coeff(f0coeff, ratio):

    icentral = len(f0coeff) // 2
    F0 = f0coeff[icentral]
    for i in range(icentral):
        F0 += f0coeff[i] * numpy.exp(-1.0 * f0coeff[i + icentral + 1] * ratio ** 2)
    return F0

#
# __* are auxiliar routines not to be exported outside
#

def __f0_interpolate_coefficients(charge, interesting_entries, charge_list, coefficient_list, verbose=True):
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
        if verbose: print("No interpolating needed: returning value for ", interesting_entries[idx])
        return coefficient_list[idx]

    # get the closer charge values, no matter of the sign

    sorted_indices = numpy.argsort(numpy.abs(charge_list_difference))
    sorted_index_0 = sorted_indices[0]
    sorted_index_1 = sorted_indices[1]
    delta_data = charge_list[sorted_index_1] - charge_list[sorted_index_0]
    delta_charge = charge - charge_list[sorted_index_0]
    delta = delta_charge / delta_data
    if verbose: print("Interpolating charge %g = %s + %g (%s - %s)" % (charge,
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
        if verbose: print("fitted: ", popt)

        return popt
    except:
        if verbose: print("Error: failed to fit coefficients for fractional charge. Returning the ones of ",
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

    raise Exception("Why are you here?", ATOM)

def __f0func(q, a1, a2, a3, a4, a5, a6, a7, a8, a9):
    return calculate_f0_from_f0coeff([a1, a2, a3, a4, a5, a6, a7, a8, a9], q)

#
# f1f2
#

def f1f2_calc_dabax(descriptor,
                    energy,
                    theta=3.0e-3,
                    F=0,
                    density=None,
                    rough=0.0,
                    verbose=True,
                    filename="f1f2_Windt.dat",
                    dabax_repository=dabax_repository,
                    interpolation_log=False):
    """
    calculate the elastic Photon-Atom anonalous f1 and f2  coefficients as a function of energy.
    It also gives the refractive index components delta and beta (n=1-delta - i beta),
    the absorption photoelectric coefficient and the reflectivities (s,p and unpolarized).
    :param descriptor: string with the element symbol or integer with Z
    :param energy: array with energies (eV)
    :param theta: array with grazing angles (rad)
    :param F: calculation flag:

           F=0 (default) returns a 2-col array with f1 and f2
           F=1  returns f1
           F=2  returns f2
           F=3  returns delta  [n = 1 -delta -i beta]
           F=4  returns betaf  [n = 1 -delta -i beta]
           F=5  returns Photoelectric linear absorption coefficient
           F=6  returns Photoelectric mass absorption coefficient
           F=7  returns Photoelectric Cross Section
           F=8  returns s-polarized reflectivity
           F=9  returns p-polarized reflectivity
           F=10  returns unpolarized reflectivity
           F=11  returns delta/betaf
    :param density: the density to be used for some calculations. If None, get it from xraylib
    :param rough: the roughness RMS in Angstroms for reflectivity calculations
    :return: a numpy array with results
    """

    energy = numpy.array(energy,dtype=float).reshape(-1)
    theta = numpy.array(theta,dtype=float).reshape(-1)

    if isinstance(descriptor,str):
        Z = xraylib.SymbolToAtomicNumber(descriptor)
        symbol = descriptor
    else:
        Z = descriptor
        symbol = xraylib.AtomicNumberToSymbol(descriptor)

    if density is None:
        density = xraylib.ElementDensity(Z)

    if verbose:
        print("f1f2_calc: using density: %f g/cm3" % density)

    # access spec file

    file1 = get_dabax_file(filename, dabax_repository=dabax_repository, verbose=verbose)

    sf = SpecFile(file1)

    flag_found = False
    for index in range(len(sf)):
        s1 = sf[index]
        name = s1.scan_header_dict["S"]
        line = " ".join(name.split())
        # print(">>>>>>%s  **%s**" % ( line, line.split(' ')[1]))
        if (line.split(' ')[1]) == descriptor:
            flag_found = True
            index_found = index

    if not flag_found:
        raise (Exception("Entry name not found: %s" % descriptor))

    zadd = sf[index_found].scan_header_dict["UF1ADD"]


    data = sf[index_found].data
    L = sf.labels(index_found)

    photon_energy = data[0, :].copy()

    if filename in ['f1f2_asf_Kissel.dat','f1f2_Chantler.dat']:
        photon_energy *= 1e3

    if filename == 'f1f2_asf_Kissel.dat':
        f1 = data[4,:].copy()
        f2 =  numpy.abs(data[1,:].copy())
    else:
        f1 = data[1, :].copy()
        f2 = data[2, :].copy()

    if interpolation_log:
        f1_interpolated = 10 ** numpy.interp(numpy.log10(energy), numpy.log10(photon_energy), numpy.log10(numpy.abs(f1)))
        f2_interpolated = 10 ** numpy.interp(numpy.log10(energy), numpy.log10(photon_energy), numpy.log10(numpy.abs(f2)))
    else:
        f1_interpolated = numpy.interp(energy, photon_energy, f1)
        f2_interpolated = numpy.interp(energy, photon_energy, f2)

    if zadd != 0: # adds Z if not included in the data
        f1_interpolated += float(zadd)

    if F == 0:   # F=0 (default) returns a 2-col array with f1 and f2
        out = numpy.zeros((2,energy.size))
        out[0,:] = f1_interpolated
        out[1,:] = f2_interpolated
        return out
    elif F == 1: # F=1  returns f1
        return f1_interpolated
    elif F == 2: # F=2  returns f2
        return f2_interpolated

    atwt = xraylib.AtomicWeight(Z)
    avogadro = codata.Avogadro
    toangstroms = codata.h * codata.c / codata.e * 1e10
    re = codata.e ** 2 / codata.m_e / codata.c ** 2 / (4 * numpy.pi * codata.epsilon_0) * 1e2  # in cm

    molecules_per_cc = density * avogadro / atwt
    wavelength = toangstroms / energy * 1e-8  # in cm
    k = molecules_per_cc * re * wavelength * wavelength / 2.0 / numpy.pi

    delta = k * f1_interpolated
    betaf = k * f2_interpolated
    mu = 4.0 * numpy.pi * betaf / wavelength

    if F == 3:   # F=3  returns delta  [n = 1 -delta -i beta]
        return delta
    elif F == 4: # F=4  returns betaf  [n = 1 -delta -i beta]
        return betaf

    elif F == 5: # F=5  returns Photoelectric linear absorption coefficient
        return mu
    elif F == 6: # F=6  returns Photoelectric mass absorption coefficient
        return mu / density
    elif F == 7: # F=7  returns Photoelectric Cross Section
        return mu / molecules_per_cc * 1e24
    elif F == 11: # F=11  returns delta/betaf
        return delta / betaf
    #
    # mirror reflectivity
    #
    alpha = 2.0 * k * f1_interpolated
    gamma = 2.0 * k * f2_interpolated
    #
    rs,rp,runp = interface_reflectivity(alpha,gamma,theta)


    if rough != 0:
        rough *= 1e-8 # to cm
        debyewaller = numpy.exp( -( 4.0 * numpy.pi * numpy.sin(theta) * rough / wavelength)**2)
    else:
        debyewaller = 1.0

    if F == 8:   # returns s-polarized reflectivity
        return rs * debyewaller
    elif F == 9: # returns p-polarized reflectivity
        return rp * debyewaller
    elif F == 10: # returns unpolarized reflectivity
        return runp * debyewaller

    raise Exception("Invalid F=%g" % F)

def atomic_weights_dabax(descriptor,
                    filename="AtomicWeights.dat",
                    return_mode=0,
                    dabax_repository=dabax_repository,
                    verbose=True,
                    ):
    """
    ; ATOMIC_WEIGHTS
    ;
    ; PURPOSE:
    ;       Returns atomic weights from DABAX.
    ;
    ; CATEGORY:
    ;       X-Ray optics. DABAX data base.
    ;
    ; CALLING SEQUENCE:
    ;       out = atomic_constants(id,file,return=return)
    ; INPUTS:
    ;       id: an identifier string (i.e. 'Si', '70Ge)
    ;
    ;       If descriptor is the symbol (e.g., Ge),
    ;         the averaged atomic mass is returned.
    ;       If descriptor contains the isotope (number of nucleons) (e.g., 70Ge),
    ;         the atomic mass for the isotope is returned.
    ;
    ;       filename = the DABAX  inout file (default AtomicWeights.dat)


    """

    if isinstance(descriptor,str):
        descriptor = [descriptor]
        descriptor_is_string = 1
    else: # is list
        descriptor_is_string = 0

    # Z = []
    # for idescriptor in descriptor:
    #     Z.append( xraylib.SymbolToAtomicNumber(idescriptor) )

    # access spec file

    file1 = get_dabax_file(filename, dabax_repository=dabax_repository, verbose=verbose)

    sf = SpecFile(file1)

    out = []

    for idescriptor in descriptor:
        flag_found = False
        index_found = []
        for index in range(len(sf)):
            s1 = sf[index]
            name = s1.scan_header_dict["S"]
            line = " ".join(name.split())
            scan_name = line.split(' ')[1]
            # print(">>>>>>%s  **%s**" % ( line, scan_name))
            if scan_name[-len(idescriptor):] == idescriptor:
                flag_found = True
                index_found.append(index)

        if not flag_found:
            raise (Exception("Entry name not found: %s" % idescriptor))

        data = sf[index_found[0]].data

        if idescriptor[0].isdigit():
            out.append(data[0,0])
        else:
            out.append(data[2,0])

    if descriptor_is_string:
        return out[0]
    else:
        return out


def atomic_symbols_dabax():
    return [
    'Vacuum','H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
    'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
    'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
    'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
    'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
    'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
    'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
    'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
    'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
    'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
    'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds',
    'Rg','Uub','Uut','Uuq','Uup','Uuh','Uus','Uuo']

def atomic_names_dabax():
    return [
            'Vacuum',
            'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron',
            'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine', 'Neon',
            'Sodium', 'Magnesium', 'Aluminum', 'Silicon', 'Phosphorus',
            'Sulfur', 'Chlorine', 'Argon', 'Potassium', 'Calcium',
            'Scandium', 'Titanium', 'Vanadium', 'Chromium', 'Manganese',
            'Iron', 'Cobalt', 'Nickel', 'Copper', 'Zinc',
            'Gallium', 'Germanium', 'Arsenic', 'Selenium', 'Bromine',
            'Krypton', 'Rubidium', 'Strontium', 'Yttrium', 'Zirconium',
            'Niobium', 'Molybdenum', 'Technetium', 'Ruthenium', 'Rhodium',
            'Palladium', 'Silver', 'Cadmium', 'Indium', 'Tin',
            'Antimony', 'Tellurium', 'Iodine', 'Xenon', 'Cesium',
            'Barium', 'Lanthanum', 'Cerium', 'Praseodymium', 'Neodymium',
            'Promethium', 'Samarium', 'Europium', 'Gadolinium', 'Terbium',
            'Dysprosium', 'Holmium', 'Erbium', 'Thulium', 'Ytterbium',
            'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium',
            'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury',
            'Thallium', 'Lead', 'Bismuth', 'Polonium', 'Astatine',
            'Radon', 'Francium', 'Radium', 'Actinium', 'Thorium',
            'Protactinium', 'Uranium', 'Neptunium', 'Plutonium', 'Americium',
            'Curium', 'Berkelium', 'Californium', 'Einsteinium', 'Fermium',
            'Mendelevium', 'Nobelium', 'Lawrencium', 'Rutherfordium', 'Dubnium',
            'Seaborgium', 'Bohrium', 'Hassium', 'Meitnerium', 'Darmstadtium',
            'Roentgenium', 'Ununbium', 'Ununtrium', 'Ununquadium', 'Ununpentium',
            'Ununhexium', 'Ununseptium', 'Ununoctium']




if __name__ == "__main__":

    #
    # at ESRF use one of these. Otherwise comment (use then the default at the top of this file)
    #
    # dabax_repository = "http://ftp.esrf.fr/pub/scisoft/DabaxFiles/"
    dabax_repository = "/scisoft/DABAX/data"

    #
    # crystal tests
    #
    if False:
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

    if False:
        print(crystal_parser(filename='Crystals.dat', entry_name='YB66', dabax_repository=dabax_repository))

        # compare with xraylib
        xdabax = crystal_parser(filename='Crystals.dat', entry_name='Si', dabax_repository=dabax_repository)

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
    if False:
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

    #
    # f0 another test
    #
    if False:
        from dabax_util import calculate_f0_from_f0coeff, f0_with_fractional_charge
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

    #
    # f1f2 tests
    #

    if False:
        files_f1f2 = [
            "f1f2_asf_Kissel.dat",
            # "f1f2_BrennanCowan.dat",
            # "f1f2_Chantler.dat",
            # "f1f2_CromerLiberman.dat",
            # "f1f2_EPDL97.dat",
            # "f1f2_Henke.dat",
            # "f1f2_Sasaki.dat",
            # "f1f2_Windt.dat",
        ]

        for file_f1f2 in files_f1f2:
            for F in range(12):
                print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>F=", F)
                a_dabax = f1f2_calc_dabax("Si", 10000.0, F=F, theta=2e-3, verbose=0,
                                      filename=file_f1f2, dabax_repository=dabax_repository )
                a_xraylib = f1f2_calc      ("Si", 10000.0, F=F, theta=2e-3, verbose=0)
                diff = (numpy.array(a_dabax) - numpy.array(a_xraylib)) / numpy.array(a_xraylib)
                print("dabax: ", file_f1f2, a_dabax)
                print("xraylib: ", a_xraylib)
                print("diff: ", numpy.abs( diff.sum()))
                assert (numpy.abs( diff.sum()) < 0.11 )

    if True:

        print("Ge, Si: ", atomic_weights_dabax(["Ge","Si"],dabax_repository=dabax_repository))
        print("70Ge: ", atomic_weights_dabax("70Ge",dabax_repository=dabax_repository))
        # print("40Ge: ", atomic_weights_dabax("40Ge",dabax_repository=dabax_repository))

        print(atomic_symbols_dabax()[14], atomic_names_dabax()[14])