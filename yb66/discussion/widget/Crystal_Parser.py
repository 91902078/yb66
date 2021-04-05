import re
import numpy
#import xraylib
import sys
import os

def Crystal_Parser(filename='YB66.dat'):
    """
    X.J. YU, xiaojiang@nus.edu.sg
    parse a complex crystal structure file, into a dictionary
    return a dictionary containing crystal infomation
    """

    header_end = False          #True, end of header parsing
    cryst = {'name':'YB66'}     #returned dictionary like that one created by xraylib.Crystal_GetCrystal(descriptor)
    AnisoItem = {'Name':'       ','start':0,'end':0,'beta11':0.0,'beta22':0.0,'beta33':0.0,'beta12':0.0,'beta13':0.0,'beta23':0.0}
    Aniso = numpy.array([AnisoItem.copy() for i in range(100)])
    AtomItem = {'AtomicName': '       ','Zatom':0,'fraction':0.0,'x':0.0,'y':0.0,'z':0.0,'SeqNo':0,'charge':0.0}
    n = 0           #index anisotropic lines
    n_atom = 0      #index atom list lines
    # open file
    if not os.path.isfile(filename):
        filename = os.path.dirname(os.path.abspath(__file__)) + '\\YB66.dat'
    print(filename)
    fo = open(filename, "r")
    while True:
        line = fo.readline().strip()
        if line[-1:] == '\n':   
            line = line[0:-1]   #discard last \n character
        if len(line) == 0:      #end of file, exit
            break                          
        a = line.split()
        if header_end == False:
            #process header info
            if a[0][:2]== '#L': #L  AtomicNumber  Fraction  X  Y  Z
                header_end = True
                tot_entry = 0
                curr_pos = fo.tell()
                for x in fo.readlines():
                    if x[0][:2] == '#S': #next crystal begin
                        break
                    if len(x) == 1 and x[-1:] == '\n':  #not allow empty line in the atom list
                        break
                    tot_entry = tot_entry +1
                fo.seek(curr_pos,0)     #rewind to begin of first line of atom list
                Atom = numpy.array([AtomItem.copy() for i in range(tot_entry)])     # create atom dictionary

            elif a[0][:2] == '#S':  #S 5 YB66
                cryst['name'] = str(a[2])
                
            elif a[0][:6] == '#UCELL':  #UCELL 23.44 23.44 23.44 90.0000 90.0000 90.0000
                cryst['a'] = float(a[1])
                cryst['b'] = float(a[2])
                cryst['c'] = float(a[3])
                cryst['alpha'] = float(a[4])
                cryst['beta'] =  float(a[5])
                cryst['gamma'] = float(a[6])
                alpha = float(a[4]) * numpy.pi/180
                beta =  float(a[5]) * numpy.pi/180
                gamma = float(a[6]) * numpy.pi/180
                cryst['volume'] = float(a[1]) * float(a[2]) * float(a[3]) * \
                  numpy.sqrt( (1 - numpy.cos(alpha)**2 - numpy.cos(beta)**2 - numpy.cos(gamma)**2) +  \
                         2 * numpy.cos(alpha) * numpy.cos(beta) * numpy.cos(gamma))      #for cubic=a*b*c
            elif  a[0][:12] == '#UANISO_COFF':  #UANISO_COFF_B1 1 96 0.00038 0.00044 0.00039 0 0 0 
                Aniso[n]['Name']=   str(a[0][13:])
                Aniso[n]['start']=  int(a[1])
                Aniso[n]['end']=    int(a[2])
                Aniso[n]['beta11']= float(a[3])
                Aniso[n]['beta22']= float(a[4])
                Aniso[n]['beta33']= float(a[5])
                Aniso[n]['beta12']= float(a[6])
                Aniso[n]['beta13']= float(a[7])
                Aniso[n]['beta23']= float(a[8])
                n = n + 1                       #increase anisotropic index

        else:
            #B-. 1 0.5 0.94077 0.96284
            if len(a) < 5:          #min 5 column required, end of atom list or new header section start
                break
            tmp1 = re.search('(^[0-9]*)',a[0])
            if len(tmp1.group(0)) > 0:   #numerical atomic number
                Atom[n_atom]['Zatom'] = int(a[0])
                Atom[n_atom]['AtomicName'] = ''
            else:   #atomic name string
                tmp1 = re.search('(^[a-zA-Z]*)',a[0])   #get atomic name
                Atom[n_atom]['AtomicName'] = str(a[0])
#                Atom[n_atom]['Zatom'] = int(xraylib.SymbolToAtomicNumber(tmp1.group(0)))    #atomic name to atomic number
                Atom[n_atom]['Zatom'] = SymbolToAtomicNumber(tmp1.group(0))    #atomic name to atomic number
            Atom[n_atom]['fraction'] =  float(a[1])
            Atom[n_atom]['x'] =         float(a[2])
            Atom[n_atom]['y'] =         float(a[3])
            Atom[n_atom]['z'] =         float(a[4])
            Atom[n_atom]['SeqNo'] = int(n_atom)
            if len(a) == 6:     #6 colume list, last one is carried charge
                Atom[n_atom]['charge'] = float(a[5])
            n_atom = n_atom + 1
 
    # close file
    fo.close()
    cryst['atom'] = Atom
    cryst['n_atom']= n_atom
    if n >0:
        cryst['Aniso'] = Aniso[:n]
        cryst['n_aniso']= n
    else:
        AnisoItem['name']=''
        cryst['Aniso'] =  numpy.array([AnisoItem])
        cryst['n_aniso']= 1
    return cryst

def SymbolToAtomicNumber(ATOM):
    atoms = [ 	[1,"H"],[2,"He"],[3,"Li"],[4,"Be"],[5,"B"],[6,"C"],[7,"N"],[8,"O"],[9,"F"],[10,"Ne"], 	\
        [11,"Na"],[12,"Mg"],[13,"Al"],[14,"Si"],[15,"P"],[16,"S"],[17,"Cl"],[18,"Ar"],[19,"K"],[20,"Ca"], 	\
        [21,"Sc"],[22,"Ti"],[23,"V"],[24,"Cr"],[25,"Mn"],[26,"Fe"],[27,"Co"],[28,"Ni"],[29,"Cu"],[30,"Zn"], \
        [31,"Ga"],[32,"Ge"],[33,"As"],[34,"Se"],[35,"Br"],[36,"Kr"],[37,"Rb"],[38,"Sr"],[39,"Y"],[40,"Zr"], \
        [41,"Nb"],[42,"Mo"],[43,"Tc"],[44,"Ru"],[45,"Rh"],[46,"Pd"],[47,"Ag"],[48,"Cd"],[49,"In"],[50,"Sn"], 	\
        [51,"Sb"],[52,"Te"],[53,"I"],[54,"Xe"],[55,"Cs"],[56,"Ba"],[57,"La"],[58,"Ce"],[59,"Pr"],[60,"Nd"], 	\
        [61,"Pm"],[62,"Sm"],[63,"Eu"],[64,"Gd"],[65,"Tb"],[66,"Dy"],[67,"Ho"],[68,"Er"],[69,"Tm"],[70,"Yb"], 	\
        [71,"Lu"],[72,"Hf"],[73,"Ta"],[74,"W"],[75,"Re"],[76,"Os"],[77,"Ir"],[78,"Pt"],[79,"Au"],[80,"Hg"], 	\
        [81,"Tl"],[82,"Pb"],[83,"Bi"],[84,"Po"],[85,"At"],[86,"Rn"],[87,"Fr"],[88,"Ra"],[89,"Ac"],[90,"Th"], 	\
        [91,"Pa"],[92,"U"],[93,"Np"],[94,"Pu"],[95,"Am"],[96,"Cm"],[97,"Bk"],[98,"Cf"],[99,"Es"],[100,"Fm"], 	\
        [101,"Md"],[102,"No"],[103,"Lr"],[104,"Rf"],[105,"Db"],[106,"Sg"],[107,"Bh"] 	]
    for a in atoms:
        if a[1] == ATOM:
            return int(a[0])

    raise Exception("Why are you here?")    