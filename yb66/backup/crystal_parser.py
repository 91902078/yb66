import numpy
from dabax_util import get_dabax_file
from silx.io.specfile import SpecFile
from f0coeffs_fit import Crystal_get_f0coeffs
from SymbolToFromAtomicNumber import SymbolToFromAtomicNumber

from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_metrictensor


def crystal_parser(filename='Crystals.dat', entry_name='YB66'):
    """
    X.J. YU, xiaojiang@nus.edu.sg, M. Sanchez del Rio srio@esrf.eu

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
            atom.append({'AtomicName': SymbolToFromAtomicNumber(int(cell_data[0,i]))+f'%+g'%cell_data[5, i],  
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

    #process charged f0 coefficients, with 6 column, f0 coefficients is not tabulated
    #stroe a f0 coefficients in a key f0coeffs, no effect for 5 column
    if cell_data.shape[0] == 6: # standard 6 columns
        AtomicChargeList = {}
        #first row is atomic number, it is integer
        UniqueAtomicNumber = [int(x) for x in list(sorted(set(cell_data[0,:])))]         
        for x in  UniqueAtomicNumber:
            AtomicChargeList[str(x)]= [] 
        for i,x in enumerate(cell_data[0,:]):
            if cell_data[5,i] not in AtomicChargeList[str(int(x))]:
                AtomicChargeList[str(int(x))].append(cell_data[5,i])      #Charge value
        cryst['f0coeffs'] =Crystal_get_f0coeffs(AtomicChargeList.items())
    

    return cryst


if __name__ == "__main__":


    import xraylib

    if True:
        #
        # testing new parser against xraylib...
        #
        descriptor = "Ge"
        dict_xraylib = xraylib.Crystal_GetCrystal(descriptor)
        dict_parser = crystal_parser(filename='Crystals.dat', entry_name=descriptor)


        for key in dict_xraylib.keys():
            if isinstance(dict_xraylib[key],list):
                for i in range(len(dict_xraylib[key])):
                    print("    ", i, dict_xraylib[key][i], dict_parser[key][i])
            else:
                print(key, dict_xraylib[key], dict_parser[key])


    # #
    # # new entry
    # #
    #
    # # local
    # dict = crystal_parser(filename='YB66_2.dat', entry_name='YB66')
    # print(dict)
    #
    #
    # # remote
    #
    # # delete old file if exists
    # try:
    #     import os
    #     os.remove("Crystals.dat")
    # except:
    #     pass
    #
    # dict = crystal_parser(filename='Crystals.dat', entry_name='YB66')
    # print(dict)
