import os
import numpy
from urllib.request import urlretrieve
from silx.io.specfile import SpecFile
from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_metrictensor
from scipy.optimize import curve_fit

"""
X.J. YU, xiaojiang@nus.edu.sg, M. Sanchez del Rio srio@esrf.eu
"""

def get_dabax_file(filename, url="http://ftp.esrf.eu/pub/scisoft/DabaxFiles/"):

    try:
        if os.path.exists(filename):
            print("File exists: %s " % filename)
        else:
            filepath, http_msg = urlretrieve(url + filename,
                        filename=filename,
                        reporthook=None,
                        data=None)

            print("File %s downloaded from %s" % (filepath, url + filename))
        return True
    except:
        return False

#
# f0
#

def func(q, a1, a2, a3, a4, a5, a6, a7, a8, a9):
    return get_f0_from_f0coeff([a1, a2, a3, a4, a5, a6, a7, a8, a9], q)

def get_f0_coeffs_from_dabax_file(entry_name="Y3+", filename="f0_InterTables.dat"):

    if getattr(get_f0_coeffs_from_dabax_file,'sf') is not None:
        sf = getattr(get_f0_coeffs_from_dabax_file,'sf')
    else:    
        error_flag = get_dabax_file(filename)
        if error_flag == False:
            raise(FileNotFoundError)
        sf = SpecFile(filename)
        get_f0_coeffs_from_dabax_file.sf = sf

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
setattr(get_f0_coeffs_from_dabax_file,'sf',None)

def get_f0_from_f0coeff(f0coeff, ratio):

    icentral = len(f0coeff) // 2
    F0 = f0coeff[icentral]
    for i in range(icentral):
        F0 += f0coeff[i] * numpy.exp(-1.0 * f0coeff[i + icentral + 1] * ratio ** 2)
    return F0

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

#
# crystal
#
def crystal_parser(filename='Crystals.dat', entry_name='YB66'):
    """
    parse a complex crystal structure file into a dictionary (like xraylib)
    return a dictionary containing crystal infomation
    """

    error_flag = get_dabax_file(filename)


    sf = SpecFile(filename)

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
    alpha = float(a[3]) * numpy.pi / 180
    beta =  float(a[4]) * numpy.pi / 180
    gamma = float(a[5]) * numpy.pi / 180


    # I do not know is this is valid for all crystal systems...
    # volume = float(a[1]) * float(a[2]) * float(a[3]) * \
    #          numpy.sqrt((1 - numpy.cos(alpha) ** 2 - numpy.cos(beta) ** 2 - numpy.cos(gamma) ** 2) + \
    #                     2 * numpy.cos(alpha) * numpy.cos(beta) * numpy.cos(gamma))  # for cubic=a*b*c

    volume = bragg_metrictensor(float(a[0]), float(a[1]), float(a[2]), float(a[3]), float(a[4]), float(a[5]), RETURN_VOLUME=1)

    cryst['volume'] = volume

    cell_data = numpy.array(sf[index_found].data)

    cryst['n_atom'] = cell_data.shape[1]
    atom = []

    for i in range(cell_data.shape[1]):
        if cell_data.shape[0] == 5: # standard 5 columns
            atom.append({'Zatom':int(cell_data[0,i]),
                         'fraction':cell_data[1,i],
                         'x': cell_data[2,i],
                         'y': cell_data[3, i],
                         'z': cell_data[4, i],})
        else: # 6 columns (charge)
            #'AtomicName' required to compatible my current code
            s = symbol_to_from_atomic_number(int(cell_data[0,i]))
            if cell_data[5, i] != 0:  #charged
                s = s + f'%+g'%cell_data[5, i]
            atom.append({'AtomicName': s,  
                         'Zatom':int(cell_data[0,i]),
                         'fraction':cell_data[1,i],
                         'x': cell_data[2,i],
                         'y': cell_data[3, i],
                         'z': cell_data[4, i],
                         'charge': cell_data[5, i],})

    cryst['atom'] = atom
    cryst['cpointer'] = None

    # TODO: Get and store anisotropic coeffs
    ANISO_KEY = "UANISO_COFF"   #prefix for a line with anisotropic coefficients
    #tmp = sf[index_found].scan_header_dict["UANISO_COFF_B1"]
    d = sf[index_found].scan_header_dict
    AnisoItem = {'Name': '       ', 'start': 0, 'end': 0, 'beta11': 0.0, 'beta22': 0.0, 'beta33': 0.0,
                 'beta12': 0.0, 'beta13': 0.0, 'beta23': 0.0}

    a=[ (x, d[x].split()) for x in d if x[:len(ANISO_KEY)] == ANISO_KEY]
    if len(a) >0:       #found Anisotropic coefficients in the header, process it
        a=sorted(a,key=lambda x:int(x[1][0]),reverse=False)     #sort 'Start' ascendant, avoid order changed by the SpecFile
        n = 0
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
            if n ==0:
                Aniso = numpy.array([AnisoItem.copy()])
            else:
                Aniso = numpy.append(Aniso,[AnisoItem.copy()])
            n = n + 1
        cryst['Aniso'] = Aniso      #if having key 'Ansio' when there is anisotropic data,otherwise no
        cryst['n_aniso']= n

 
    return cryst


def crystal_atnum(list_AtomicName, unique_AtomicName, unique_Zatom,list_fraction, f0coeffs):
    """
    To get the atom and fractional factor in diffierent sites
    list_AtomicName:  list of all atoms in the crystal
    unique_AtomicName:  list of unique atomicname in the list
    unique_Zatom:    list of unique atomic number
    list_fraction: list of unique fractial factor

    return: (num_e, fract, n_atom, n_name)
    (number of electrons, fraction, atomic sites, Unique name)
    """
    import re
    from orangecontrib.xoppy.util.xoppy_xraylib_util import f0_xop

    num_e = []
    fract = []
    n_atom = []
    n_ATUM = []
    n_name = []
    for k,x in enumerate(unique_AtomicName):
        tmp1 = re.search('(^[a-zA-Z]*)',x)
        if tmp1.group(0) == x:   #AtomicName only, without valence info (i.e., B, Y, O)
            f0 = f0_xop(unique_Zatom[k])
        else:   
            #f0 = f0_xop(0,AtomicName=x)
            f0 = f0coeffs[x]

        icentral = int(len(f0)/2)
        F000 = f0[icentral]
        for i in range(icentral):
            F000 += f0[i]  
        a=[list_fraction[i] for i,v in enumerate(list_AtomicName) if v==x]
        fac = list(set(a))
        for y in fac:
            n = a.count(y)
            num_e.append(F000)
            fract.append(y)
            n_atom.append(n)
            n_ATUM.append(unique_Zatom[k])
            n_name.append(x)

    return num_e.copy(), fract.copy(), n_atom.copy(),n_ATUM.copy(), n_name.copy()

def Crystal_GetCrystalsList():
    """
    get crystal names from crystals.dat
    """
    sf = SpecFile('Crystals.dat')
    crystals = []
    for index in range(len(sf)):
        s1 = sf[index]
        name = s1.scan_header_dict["S"]
        crystals.append(name.split(' ')[1])
    
    return crystals


#
# tools
#

def symbol_to_from_atomic_number(ATOM):
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

