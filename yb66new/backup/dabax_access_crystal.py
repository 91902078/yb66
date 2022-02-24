import numpy
from silx.io.specfile import SpecFile

from dabax_util import get_dabax_file


def get_crystal_data_from_dabax_file(entry_name="Si", filename="Crystal.dat"):
    error_flag = get_dabax_file(filename)
    if error_flag == False:
        raise(FileNotFoundError)

    sf = SpecFile(filename)

    flag_found = False

    for index in range(len(sf)):
        s1 = sf[index]
        name = s1.scan_header_dict["S"]

        if name.split(' ')[1] == entry_name:
            flag_found = True
            index_found = index

    if flag_found:
        cell_parameters = sf[index_found].scan_header_dict["UCELL"]
        cell_data = numpy.array(sf[index_found].data)
        return cell_parameters, cell_data
    else:
        raise(Exception("Entry name not found: %s" % entry_name))




if __name__ == "__main__":

    filename = "Crystals.dat"

    cell_parameters, cell_data = get_crystal_data_from_dabax_file(entry_name="Si", filename=filename)

    print("Cell parameters: ", cell_parameters)
    print("Cell data: ", cell_data)
