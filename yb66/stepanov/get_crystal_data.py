from dabax.common_tools import atomic_number

def __parse_line(line, remove=[]):
    if len(remove) > 0:
        for str1 in remove:
            line = line.replace(str1, " ")
    line = " ".join(line.split())
    variables = line.split(" ")
    return variables

def get_cell(Syngony_number, line):
    var = __parse_line(line, remove=[","])

    if Syngony_number == 0:
        a = float(var[0])
        b = float(var[0])
        c = float(var[0])
        alpha = 90.0
        beta  = 90.0
        gamma = 90.0
    elif Syngony_number == 1:
        a = float(var[0])
        b = float(var[0])
        c = float(var[0])
        alpha = 90.0
        beta  = 90.0
        gamma = 90.0
    elif Syngony_number == 2:
        a = float(var[0])
        b = float(var[0])
        c = float(var[1])
        alpha = 90.0
        beta  = 90.0
        gamma = 90.0
    elif Syngony_number == 3:
        a = float(var[0])
        b = float(var[0])
        c = float(var[1])
        alpha = 90.0
        beta  = 90.0
        gamma = 90.0
    elif Syngony_number == 4:
        a = float(var[0])
        b = float(var[0])
        c = float(var[0])
        alpha = float(var[1])
        beta = 90.0
        gamma = 90.0
    elif Syngony_number == 5:
        a = float(var[0])
        b = float(var[1])
        c = float(var[2])
        alpha = 90.0
        beta  = 90.0
        gamma = 90.0
    elif Syngony_number == 6:
        a = float(var[0])
        b = float(var[1])
        c = float(var[2])
        alpha = 90.0
        beta = float(var[3])
        gamma = 90.0
    elif Syngony_number == 7:
        a = float(var[0])
        b = float(var[1])
        c = float(var[2])
        alpha = float(var[3])
        beta  = float(var[4])
        gamma = float(var[5])
    elif Syngony_number == 8:
        a     = None
        b     = None
        c     = None
        alpha = None
        beta  = None
        gamma = None

    return a,b,c, alpha,beta, gamma

def get_crystal_data(filename, verbose=0):

    # see https://x-server.gmca.aps.anl.gov/structure_submission_example.html

    if verbose: print(">>>>> file: ", filename)

    f = open(filename, 'r')
    lines = f.read().splitlines()
    f.close()

    nlines = len(lines)

    istart = -1
    header = []
    comments = []
    symbols = []
    natoms = []
    fractions = []
    COORDINATES = []
    name = ""
    a,b,c,alpha,beta,gamma = 0,0,0,0,0,0

    is_header = 1

    for i in range(nlines):
        line = lines[i]
        if len(line) > 0:
            if line[0] == ";":
                if is_header:
                    header.append(line)
                else:
                    comments.append(line)
            if line[0] == "#":
                istart = i
                name = line[1:]
                is_header = 0

    #
    # get name and dimensions
    #
    if lines[istart+1][0] == ";":
        itmp = istart + 2
    else:
        itmp = istart + 1

    var = __parse_line(lines[itmp], remove=[","])

    Nr_of_components  = int(var[0])
    Syngony_number = int(var[1])
    if len(var) == 3:
        rho = float(var[2])
    else:
        rho = None

    if verbose: print(">>>Nr_of_components, Syngony_number, rho: ", Nr_of_components, Syngony_number, rho)

    if (Nr_of_components > 0):
        # get a,b,c,alpha,beta,gamma

        if Syngony_number < 8:
            istart = itmp + 1
            while lines[istart][0] == ";":
                istart += 1
            a,b,c,alpha,beta,gamma = get_cell(Syngony_number, lines[istart])



        lines_without_comments = []
        for i,line in enumerate(lines):
            if i > istart and line[0] != ";":
                lines_without_comments.append(line)


        istart = 0
        while istart >= 0:
            try:
                if verbose: print(">>>> line SYMBOL: ", lines_without_comments[istart])
                symbol = (lines_without_comments[istart][0:2]).strip()
                if len(symbol) == 2:
                    if symbol[1] == "-":
                        symbol = symbol[0]
                # symbol = lines_without_comments[istart]
                # if symbol.find(";") > 0:
                #     symbol = symbol[:symbol.find(";")].strip()
                if verbose: print(">>>> SYMBOL: ", symbol)
                symbols.append(symbol)
                if verbose: print(">>>> line NATOM: ", lines_without_comments[istart+1])
                var = __parse_line(lines_without_comments[istart+1], remove=[","])
                natoms.append(int(var[0]))
                if len(var) == 2:
                    fractions.append(float(var[1]))
                else:
                    fractions.append(1.0)
                istart = istart+2
                coordinates = []
                for i in range(natoms[-1]):
                    coordinates.append(lines_without_comments[istart+i])
                COORDINATES.append(coordinates)

                istart += natoms[-1]
            except:
                istart = -1
            # variables = __parse_line(lines[line_index])

    return {"name":name,
            "header":header,
            "symbols":symbols,
            "natoms":natoms,
            "fractions":fractions,
            "Nr_of_components":Nr_of_components,
            "Syngony_number":Syngony_number,
            "rho":rho,
            "COORDINATES":COORDINATES,
            "a":a, "b":b, "c":c, "alpha":alpha, "beta":beta, "gamma":gamma}

def add_crystal(out, filename):
    """
#S 19 Muscovite
#UCELL 5.18900      8.99500      20.0970      90.0000      95.1800 90.0000
#USYSTEM Monoclinic
#UREF R.W.G. Wyckoff, Crystal Structures, Interscience Publ. 1965, 4, p.  346
#UCOMMENT KAl2(AlSi3)O10(OH)2
#UCOMMENT The muscovite or mica is composed by K Al Si O and H atoms in 10
#UCOMMENT different sites  making 13 independent element-sites. It belongs
#UCOMMENT to the monoclinic system with abc=5.189 8.995 20.097 beta=95.18
#N  5
#L  AtomicNumber  Fraction  X  Y  Z
19 1.00    0.0000  0.1016  0.2500
19 1.00    0.5000  0.6016  0.2500
19 1.00   -0.0000 -0.1016 -0.2500
19 1.00   -0.5000 -0.6016 -0.2500
13 1.00    0.2484  0.0871  0.0016
    """
    f = open(filename,'a')
    f.write("\n#S %d %s\n" % (atomic_number(out["symbols"][0]), out["name"]))
    f.write("#UCELL %g  %g  %g  %g  %g  %g\n" % (
        out["a"],
        out["b"],
        out["c"],
        out["alpha"],
        out["beta"],
        out["gamma"],
    ))
    f.write("#UCOMMENT: data URL: https://x-server.gmca.aps.anl.gov/cgi/x0h_dump.pl?code=%s\n" % out["name"])
    for line in out["header"]:
        f.write("#UCOMMENT %s\n" % line)

    f.write("#N  5\n#L  AtomicNumber  Fraction  X  Y  Z\n")
    for i in range(len(out["symbols"])):
        for j in range(len(out["COORDINATES"][i])):
            f.write("%d  %g  %s\n" % (atomic_number(out["symbols"][i]), out["fractions"][i], out["COORDINATES"][i][j] ))
    f.close()
    print(">>>> added to file: ", filename)

def add_header(filename):
    txt = """#F Crystals_xrayserver.dat
#UT Crystals_xrayserver: Crystal structures with atomic coordinates in unit cell and cell dimensions
#UIDL xfh
#UD
#UD  This file contains the list of crystals available Sergey Stepanov's X-Ray Server
#UD  https://x-server.gmca.aps.anl.gov/
#UD S. Stepanov, "X-ray server: an online resource for simulations of X-ray diffraction and scattering".
#UD In: "Advances in Computational Methods for X-ray and Neutron Optics", Ed. M.Sanches del Rio; Proceedings SPIE, v.5536, p.16-26, (2004).
#UD
#UD The crystal structures are arranged using as a scan
#UD identifier (#S) the atomic number of the heavier atom
#UD present in the crystal. Several crystal are possible with
#UD the same scan id, but they do not conflict.
#UD  
#UD The following keyword contain other information:
#UD #UCOMMENT  comment 
#UD #UCELL a b c alpha beta gamma  
#UD        The unit cell dimensions (A and deg) (*MANDATORY, IT MUST EXIST*)
#UD
#UD  
#UD  Data columns:  
#UD  AtomicNumber  Fraction  X  Y  Z Biso
#UD  
#UD
#UD  This file has been created with the python code that can be found at: 
#UD     https://github.com/91902078/yb66/tree/main/yb66/stepanov
#UD
    """

    f = open(filename,'a')
    f.write(txt)
    f.close()

if __name__ == "__main__":

    import os

    f = open("crystal_list.dat", 'r')
    crystals = f.read().splitlines()
    f.close()

    # crystals = ["Mica"]

    crystal_file = "Crystals_xrayserver.dat"
    os.system("rm %s" % crystal_file)

    add_header(crystal_file)
    for crystal in crystals:
        crystal = crystal.replace("(","").replace(")","")
        print("------------------------------------ crystal: ", crystal)
        out = get_crystal_data("downloads/%s" % crystal, verbose=0)
        if (out["Nr_of_components"] > 0) and (out["Syngony_number"] < 8):
            # print (out)
            add_crystal(out, crystal_file)
        else:
            print("           Not used!")

    print("File written to disk: %s" % crystal_file)