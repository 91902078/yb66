import re
import numpy
import xraylib

def Crystal_Parser(filename='YB66.dat'):
    """
    X.J. YU, xiaojiang@nus.edu.sg
    parse a complex crystal structure file, into a dictionary
    return a dictionary containing crystal infomation
    """

    header_end = False          #True, end of header parsing
    cryst = {}                  #returned dictionary like that one created by xraylib.Crystal_GetCrystal(descriptor)
    AnisoItem = {'name':'       ','start':0,'end':0,'beta11':0.0,'beta22':0.0,'beta33':0.0,'beta12':0.0,'beta13':0.0,'beta23':0.0}
    Aniso = numpy.array([AnisoItem.copy() for i in range(100)])
    AtomItem = {'AtomicName': '       ','Zatom':0,'fraction':0.0,'x':0.0,'y':0.0,'z':0.0,'SeqNo':0,'charge':0.0}
    n = 0           #index anisotropic lines
    n_atom = 0      #index atom list lines
    # open file
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
                pass
            elif a[0][:6] == '#UCELL':  #UCELL 23.44 23.44 23.44 90.0000 90.0000 90.0000
                cryst['a'] = float(a[1])
                cryst['b'] = float(a[2])
                cryst['c'] = float(a[3])
                cryst['alpha'] = float(a[4])
                cryst['beta'] =  float(a[5])
                cryst['gamma'] = float(a[6])
                cryst['volume'] = float(a[1])*float(a[2])*float(a[3])       #for cubic=a*b*c
            elif  a[0][:12] == '#UANISO_COFF':  #UANISO_COFF_B1 1 96 0.00038 0.00044 0.00039 0 0 0 
                Aniso[n]['name']=a[0][14:]
                Aniso[n]['start']=    int(a[1])
                Aniso[n]['end']=      int(a[2])
                Aniso[n]['beta11']= float(a[3])
                Aniso[n]['beta22']= float(a[4])
                Aniso[n]['beta33']= float(a[5])
                Aniso[n]['beta12']= float(a[6])
                Aniso[n]['beta13']= float(a[7])
                Aniso[n]['beta23']= float(a[8])
                n = n + 1                       #increase anisotropic index

        else:
            #B-. 1 0.5 0.94077 0.96284
            if len(a) < 5:          #min 5 column required
                break
            tmp1 = re.search('(^[0-9]*)',a[0])
            if len(tmp1.group(0)) > 0:   #numerical atomic number
                Atom[n_atom]['Zatom'] = int(a[0])
                Atom[n_atom]['AtomicName'] = ''
            else:   #atomic name string
                tmp1 = re.search('(^[a-zA-Z]*)',a[0])   #get atomic name
                Atom[n_atom]['AtomicName'] = a[0]
                Atom[n_atom]['Zatom'] = int(xraylib.SymbolToAtomicNumber(tmp1.group(0)))    #atomic name to atomic number
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
    return cryst
