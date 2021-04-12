import os
import numpy
from urllib.request import urlretrieve
from silx.io.specfile import SpecFile
from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_metrictensor
from f0coeffs_fit import crystal_get_f0_coeffs
from symbol_to_from_atomic_number import symbol_to_from_atomic_number

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

def get_f0_coeffs_from_dabax_file(entry_name="Y3+", filename="f0_InterTables.dat"):
    error_flag = get_dabax_file(filename)
    if error_flag == False:
        raise(FileNotFoundError)

    sf = SpecFile(filename)

    flag_found = False

    for index in range(len(sf)):
        s1 = sf[index]
        name = s1.scan_header_dict["S"]

        if name.split('  ')[1] == entry_name:
            flag_found = True
            index_found = index

    if flag_found:
        return numpy.array(sf[index_found].data)[:,0]
    else:
        raise(Exception("Entry name not found: %s" % entry_name))


def get_f0_from_f0coeff(f0coeff, ratio):

    icentral = len(f0coeff) // 2
    F0 = f0coeff[icentral]
    for i in range(icentral):
        F0 += f0coeff[i] * numpy.exp(-1.0 * f0coeff[i + icentral + 1] * ratio ** 2)
    return F0

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

    return: num_e, fract, n_atom, list of number of electrons for atom with same fractional factor, and  corresponding fractional factor, atomic number
    """
    import re
    from orangecontrib.xoppy.util.xoppy_xraylib_util import f0_xop

    num_e = []
    fract = []
    n_atom = []
    n_ATUM = []
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

    return num_e.copy(), fract.copy(), n_atom.copy(),n_ATUM.copy()

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
    