import numpy
from dabax_util import get_dabax_file
from silx.io.specfile import SpecFile

from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_metrictensor


def crystal_parser(filename='YB66_1.dat', entry_name='YB66'):
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
            atom.append({'Zatom':cell_data[0,i],
                         'fraction':cell_data[1,i],
                         'x': cell_data[2,i],
                         'y': cell_data[3, i],
                         'z': cell_data[4, i],})
        else: # 6 columns (charge)
            atom.append({'Zatom':cell_data[0,i],
                         'fraction':cell_data[1,i],
                         'x': cell_data[2,i],
                         'y': cell_data[3, i],
                         'z': cell_data[4, i],
                         'charge': cell_data[5, i],})

    cryst['atom'] = atom
    cryst['cpointer'] = None

    # TODO: Get and store anisotropic coeffs
    try:
        tmp = sf[index_found].scan_header_dict["UANISO_COFF_B1"]

        AnisoItem = {'Name': '       ', 'start': 0, 'end': 0, 'beta11': 0.0, 'beta22': 0.0, 'beta33': 0.0,
                     'beta12': 0.0, 'beta13': 0.0, 'beta23': 0.0}


        # AnisoItem['Name'] = str(a[0][13:])  # get site atom name starting from 13th character 'B1', etc
        # AnisoItem['start'] = int(a[1])
        # AnisoItem['end'] = int(a[2])
        # AnisoItem['beta11'] = float(a[3])
        # AnisoItem['beta22'] = float(a[4])
        # AnisoItem['beta33'] = float(a[5])
        # AnisoItem['beta12'] = float(a[6])
        # AnisoItem['beta13'] = float(a[7])
        # AnisoItem['beta23'] = float(a[8])
    except:
        pass


    return cryst


if __name__ == "__main__":


    import xraylib

    if False:
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


    #
    # new entry
    #

    # local
    dict = crystal_parser(filename='YB66_2.dat', entry_name='YB66')
    print(dict)


    # remote

    # delete old file if exists
    try:
        import os
        os.remove("Crystals.dat")
    except:
        pass

    dict = crystal_parser(filename='Crystals.dat', entry_name='YB66')
    print(dict)
