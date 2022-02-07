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

    is_header = 1
    for i in range(nlines):
        line = lines[i]
        if line[0] == ";":
            if is_header:
                header.append(line)
            else:
                comments.append(line)
        if line[0] == "#":
            istart = i
            name = line[1:]
            is_header = 0

    iend = i

    #
    # get name and dimensions
    #
    if line[istart+1][0] == ";":
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

    # get a,b,c,alpha,beta,gamma
    istart = itmp + 1
    a,b,c,alpha,beta,gamma = get_cell(Syngony_number, lines[istart])



    lines_without_comments = []
    for i,line in enumerate(lines):
        if i > istart and line[0] != ";":
            lines_without_comments.append(line)


    istart = 0
    while istart >= 0:
        try:
            if verbose: print(">>>> line SYMBOL: ", lines_without_comments[istart])
            symbols.append(lines_without_comments[istart])
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
    f.write("#S %d %s\n" % (atomic_number(out["symbols"][0]), out["name"]))
    f.write("#UCELL %g  %g  %g  %g  %g  %g\n" % (
        out["a"],
        out["b"],
        out["c"],
        out["alpha"],
        out["beta"],
        out["gamma"],
    ))
    for line in out["header"]:
        f.write("#UCOMMENT %s\n" % line)

    f.write("#N  5\n#L  AtomicNumber  Fraction  X  Y  Z\n")
    for i in range(len(out["symbols"])):
        for j in range(len(out["COORDINATES"][i])):
            f.write("%d  %g  %s\n" % (atomic_number(out["symbols"][i]), out["fractions"][i], out["COORDINATES"][i][j] ))
    f.close()
    print(">>>> added to file: ", filename)


if __name__ == "__main__":

    import os

    f = open("crystal_list.dat", 'r')
    crystals = f.read().splitlines()
    f.close()

    crystals = ["Mica"]

    os.system("rm tmp.dat")
    for crystal in crystals:
        out = get_crystal_data("downloads/%s" % crystal, verbose=1)
        print (out)
        add_crystal(out, "tmp.dat")