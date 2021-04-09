import os
import numpy
from urllib.request import urlretrieve
from silx.io.specfile import SpecFile


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

