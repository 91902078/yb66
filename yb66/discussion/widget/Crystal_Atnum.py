from orangecontrib.xoppy.util.xoppy_xraylib_util import f0_xop
import re

def Crystal_Atnum(list_AtomicName, unique_AtomicName, unique_Zatom,list_fraction):
    """
    XJ. Yu, xiaojing@nus.edu.sg

    To get the atom and fractional factor in diffierent sites
    list_AtomicName:  list of all atoms in the crystal
    unique_AtomicName:  list of unique atomicname in the list
    unique_Zatom:    list of unique atomic number
    list_fraction: list of unique fractial factor

    return: num_e, fract, n_atom, list of number of electrons for atom with same fractional factor, and  corresponding fractional factor, atom number
    """

    num_e = []
    fract = []
    n_atom = []
    n_ATUM = []
    for k,x in enumerate(unique_AtomicName):
        tmp1 = re.search('(^[a-zA-Z]*)',x)
        if tmp1.group(0) == x:   #AtomicName only, without valence info (i.e., B, Y, O)
            f0 = f0_xop(unique_Zatom[k])
        else:   
            f0 = f0_xop(0,AtomicName=x)

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

