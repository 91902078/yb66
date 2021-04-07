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


if __name__ == "__main__":
    from srxraylib.plot.gol import plot

    filename = "f0_InterTables.dat"

    coeffs_Yplus3 = get_f0_coeffs_from_dabax_file(entry_name="Y3+", filename=filename)
    coeffs_Y = get_f0_coeffs_from_dabax_file(entry_name="Y", filename=filename)
    coeffs_Kr = get_f0_coeffs_from_dabax_file(entry_name="Kr", filename=filename)

    print(coeffs_Yplus3, coeffs_Kr)


    ratio = numpy.linspace(0,3,1000)

    f0_Yplus3 = get_f0_from_f0coeff(coeffs_Yplus3, ratio)
    f0_Y = get_f0_from_f0coeff(coeffs_Y, ratio)
    f0_Kr = get_f0_from_f0coeff(coeffs_Kr, ratio)


    from srxraylib.plot.gol import plot
    plot(ratio, f0_Yplus3,
         ratio, f0_Y,
         ratio, f0_Kr,
         legend=["Y3+", "Y", "Z=39-3"],
         xtitle=r"q (sin $\theta$ / $\lambda$)", ytitle="f0 [electron units]",
         title=filename)

