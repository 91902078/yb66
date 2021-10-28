import numpy
import os
import scipy.constants as codata

from orangecontrib.xoppy.util.xoppy_xraylib_util import f0_xop
from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_metrictensor

from dabax_util import calculate_f0_from_f0coeff, f0_with_fractional_charge
from dabax_util import Crystal_GetCrystal
from dabax_util import atomic_symbols_dabax # __symbol_to_from_atomic_number
from dabax_util import f0_with_fractional_charge
from dabax_util import CompoundParser

# to be removed...  TODO: move the f1 f2 routines from xraylib to dabax.
import xraylib


#################################################################################
# crystal tools
#################################################################################


def bragg_calc2(descriptor="YB66", hh=1, kk=1, ll=1, temper=1.0, emin=5000.0, emax=15000.0, estep=100.0, ANISO_SEL=0,
                fileout=None,
                sourceCryst=2, # 0=xraylib, 1=dabax, 2=auto
                sourceF0=2,    # 0=xraylib, 1=dabax, 2=auto
                do_not_prototype=0, # 0=use site groups (recommended), 1=use all individual sites
                verbose=True,
                dabax_repository="",
                ):
    """
    Preprocessor for Structure Factor (FH) calculations. It calculates the basic ingredients of FH.

    :param descriptor: crystal name (as in xraylib)
    :param hh: miller index H
    :param kk: miller index K
    :param ll: miller index L
    :param temper: temperature factor (scalar <=1.0 )
    :param emin:  photon energy minimum
    :param emax: photon energy maximum
    :param estep: photon energy step
    :param ANISO_SEL: source of temperature factor:
                0: use scalar value defined in temper
                1: use isotropic value calculated from keyword UNIANISO_COFF in dabax Crystal.dat file
                2: use anisotropic value calculated from keyword UNIANISO_COFF in dabax Crystal.dat file
    :param fileout: name for the output file (default=None, no output file)
    :return: a dictionary with all ingredients of the structure factor.
    """

    output_dictionary = {}

    codata_e2_mc2 = codata.e ** 2 / codata.m_e / codata.c ** 2 / (4 * numpy.pi * codata.epsilon_0)  # in m

    # f = open(fileout,'w')

    txt = ""
    txt += "# Bragg version, Data file type\n"
    txt += "2.6 1\n"

    if sourceCryst == 0:
        cryst = xraylib.Crystal_GetCrystal(descriptor)
        if cryst is None:
            raise Exception("Crystal descriptor %s not found in xraylib" % descriptor)
    elif sourceCryst == 1:
        cryst = Crystal_GetCrystal(entry_name=descriptor, dabax_repository=dabax_repository, verbose=verbose)
    elif sourceCryst == 2:
        try:
            cryst = xraylib.Crystal_GetCrystal(descriptor)
            if cryst is None:
                raise Exception("Crystal descriptor %s not found in xraylib" % descriptor)
            sourceCryst = 0
        except:
            try:
                cryst = Crystal_GetCrystal(entry_name=descriptor, dabax_repository=dabax_repository, verbose=verbose)
                sourceCryst = 1
            except:
                raise Exception("Crystal descriptor %s not found in xraylib nor in dabax" % descriptor)

    volume = cryst['volume']

    # test crystal data - not needed
    itest = 0
    if itest:

        print("  Unit cell dimensions are %f %f %f" % (cryst['a'], cryst['b'], cryst['c']))
        print("  Unit cell angles are %f %f %f" % (cryst['alpha'], cryst['beta'], cryst['gamma']))
        print("  Unit cell volume is %f A^3" % volume)
        print("  Atoms at:")
        print("     Z  fraction    X        Y        Z")
        for i in range(cryst['n_atom']):
            atom = cryst['atom'][i]
            print("    %3i %f %f %f %f" % (atom['Zatom'], atom['fraction'], atom['x'], atom['y'], atom['z']))
        print("  ")

    volume = volume * 1e-8 * 1e-8 * 1e-8  # in cm^3
    rn = (1e0 / volume) * (codata_e2_mc2 * 1e2)

    dspacing = bragg_metrictensor(cryst['a'], cryst['b'], cryst['c'], cryst['alpha'], cryst['beta'], cryst['gamma'], HKL=[hh, kk, ll])
    dspacing *= 1e-8  # in cm

    txt += "# RN = (e^2/(m c^2))/V) [cm^-2], d spacing [cm]\n"
    txt += "%e %e \n" % (rn, dspacing)

    output_dictionary["rn"] = rn
    output_dictionary["dspacing"] = dspacing

    atom = cryst['atom']
    list_Zatom = [atom[i]['Zatom'] for i in range(len(atom))]
    number_of_atoms = len(list_Zatom)
    list_fraction = [atom[i]['fraction'] for i in range(len(atom))]
    try:
        list_charge = [atom[i]['charge'] for i in range(len(atom))]
    except:
        list_charge = [0.0] * number_of_atoms
    list_x = [atom[i]['x'] for i in range(len(atom))]
    list_y = [atom[i]['y'] for i in range(len(atom))]
    list_z = [atom[i]['z'] for i in range(len(atom))]

    #
    # calculate indices of prototypical atoms (the different atomic sites)
    #
    # AtomicNames contains a string with the atomic symbol with charge appended (if not zero)
    list_AtomicName = []
    for i in range(len(atom)):
        s = atomic_symbols_dabax()[atom[i]['Zatom']] #__symbol_to_from_atomic_number(atom[i]['Zatom'])
        if sourceCryst == 1: # charge is not available in xraylib
            if atom[i]['charge'] != 0.0:  # if charge is 0, s is symbol only, not B0, etc
                s = s + f'%+.6g' % atom[i]['charge']
        list_AtomicName.append(s)

    # calculate indices of unique_AtomicName's (including charge) sorted by Z
    unique_indexes1 = numpy.unique(list_AtomicName, return_index=True)[1]
    # sort by Z
    unique_Zatom1 = [list_Zatom[i] for i in unique_indexes1]

    ii = numpy.argsort(unique_Zatom1)
    unique_indexes = unique_indexes1[ii]

    # calculate indices of sites (using unique_AtomicName's (including charge) and fraction)
    unique_indexes_with_fraction = []
    for k in unique_indexes:
        unique_indexes_with_fraction.append(k)

        itmp = [k]
        ifraction = [list_fraction[k]]

        for i in range(len(list_AtomicName)):
            if list_AtomicName[i] == list_AtomicName[k]:
                if list_fraction[i] in ifraction:
                    pass
                else:
                    itmp.append([i])
                    ifraction.append(list_fraction[i])
                    unique_indexes_with_fraction.append(i)

    # do not prototype....
    if do_not_prototype:
        unique_indexes_with_fraction = numpy.arange(number_of_atoms)
    output_dictionary["unique_indexes_with_fraction"] = unique_indexes_with_fraction  # different with diff_pat for complex crystal

    #
    # get f0 coefficients
    #
    if True:
        f0coeffs = []
        if sourceF0 == 0:
            for i in unique_indexes_with_fraction:
                f0coeffs.append(f0_xop(atom[i]['Zatom']))
        elif sourceF0 == 1:
            for i in unique_indexes_with_fraction:
                    f0coeffs.append(f0_with_fractional_charge(atom[i]['Zatom'], atom[i]['charge'],
                                                            dabax_repository=dabax_repository, verbose=verbose) )
        elif sourceF0 == 2:
            total_charge_flag = numpy.abs(numpy.array(list_charge)).sum() # note the abs(): to be used as flag...

            if total_charge_flag != 0: # Use dabax
                for i in unique_indexes_with_fraction:
                    f0coeffs.append(f0_with_fractional_charge(atom[i]['Zatom'], atom[i]['charge'],
                                                              dabax_repository=dabax_repository, verbose=verbose))
            else: # use xraylib
                if 'AtomicName' not in atom[0].keys():
                    for i in unique_indexes_with_fraction:  #normal case come in here
                        f0coeffs.append(f0_xop(atom[i]['Zatom']))
                else:   #for case with like 'Y3+' entries in f0_xop
                    import re
                    for i in unique_indexes_with_fraction:
                        x = atom[i]['AtomicName']
                        tmp_x = re.search('(^[a-zA-Z]*)',x)
                        if tmp_x.group(0) == x:
                            f0coeffs.append(f0_xop(atom[i]['Zatom']))  #neutral atom
                        else:    
                            f0coeffs.append(f0_xop(0,AtomicName=x))    #charged atom

    else:
        libmethod=1 # 0=dabax, 1=combined(Dabax for charge!=0, xraylib for charge=0)
        if libmethod == 0:
            f0coeffs = []
            for i in unique_indexes_with_fraction:
                    f0coeffs.append(f0_with_fractional_charge(atom[i]['Zatom'], atom[i]['charge'],
                                                            dabax_repository=dabax_repository, verbose=verbose) )
        else:
            f0coeffs = []
            total_charge_flag = numpy.abs(numpy.array(list_charge)).sum() # note the abs(): to be used as flag...

            if total_charge_flag != 0: # Use dabax
                for i in unique_indexes_with_fraction:
                    f0coeffs.append(f0_with_fractional_charge(atom[i]['Zatom'], atom[i]['charge'],
                                                              dabax_repository=dabax_repository, verbose=verbose))
            else: # use xraylib
                for i in unique_indexes_with_fraction:
                    f0coeffs.append(f0_xop(atom[i]['Zatom']))


    nbatom_new = (len(unique_indexes_with_fraction))
    txt += "# Number of different element-sites in unit cell NBATOM:\n%d \n" % nbatom_new
    output_dictionary["nbatom"] = nbatom_new  # different with diff_pat for complex crystal

    txt += "# for each element-site, the number of scattering electrons (Z_i + charge_i)\n"
    for i in unique_indexes_with_fraction:
        txt += "%f " % (list_Zatom[i] + list_charge[i])
    txt += "\n"

    atnum_list = []
    for i in unique_indexes_with_fraction:
        atnum_list.append (list_Zatom[i] + list_charge[i])
    output_dictionary["atnum"] = atnum_list


    if do_not_prototype:
        txt += "# for each element-site, the occupation factor\n"
        unique_fraction = [list_fraction[i] for i in range(number_of_atoms)]
        for z in unique_fraction:
            txt += "%g " % (1)
        txt += "\n"
        output_dictionary["fraction"] = unique_fraction
    else:
        txt += "# for each element-site, the occupation factor\n"
        unique_fraction = [list_fraction[i] for i in unique_indexes_with_fraction]
        for z in unique_fraction:
            txt += "%g " % (z)
        txt += "\n"
        output_dictionary["fraction"] = unique_fraction

    #
    # Consider anisotropic temperature factor
    # X.J. Yu, slsyxj@nus.edu.sg
    # A dummy dictionary Aniso with start =0 if no aniso temperature factor input
    # start
    if 'Aniso' in cryst.keys() and cryst['Aniso'][0]['start'] > 0:  # most crystals have no Anisotropic input
        TFac = TemperFactor(1.0 / (2.0 * dspacing * 1e8), cryst['Aniso'], Miller={'h': hh, 'k': kk, 'l': ll}, \
                            cell={'a': cryst['a'], 'b': cryst['b'], 'c': cryst['c']}, n=len(atom))
        B_TFac = 1
    else:
        B_TFac = 0

    txt += "# for each element-site, the temperature factor\n"  # temperature parameter

    list_temper = []
    if ANISO_SEL == 0:
        for i in unique_indexes_with_fraction:
            txt += "%5.3f " % temper
            list_temper.append(temper)
    elif ANISO_SEL == 1:
        if B_TFac:
            for i in unique_indexes_with_fraction:
                txt += "%5.3f " %  TFac[0, i]
                list_temper.append(TFac[0, i])
        else:
            raise Exception("No crystal data to calculate isotropic temperature factor for crystal %s" % descriptor)
    elif ANISO_SEL == 2:
        if B_TFac:
            for i in unique_indexes_with_fraction:
                txt += "%5.3f " %  TFac[1, i]
                list_temper.append(TFac[1, i])
        else:
            raise Exception("No crystal data to calculate anisotropic temperature factor for crystal %s" % descriptor)

    txt += "\n"
    output_dictionary["temper"] = list_temper  # not necessary same with diff_pat

    #
    # Geometrical part of structure factor:  G and G_BAR
    #
    txt += "# for each type of element-site, COOR_NR=G_0\n"
    list_multiplicity = []
    for i in unique_indexes_with_fraction:
        #    txt += "%d "%list_AtomicName.count(z)
        zz = list_AtomicName[i]
        fraction = list_fraction[i]
        count = 0
        for j in range(len(list_Zatom)):
            if (list_AtomicName[j] == zz) and (list_fraction[j] == fraction): count += 1
        if do_not_prototype:
            txt += "%d " % 1
            list_multiplicity.append(1)
        else:
            txt += "%d " % count
            list_multiplicity.append(count)
    txt += "\n"
    output_dictionary["G_0"] = list_multiplicity


    txt += "# for each type of element-site, G and G_BAR (both complex)\n"
    list_g = []
    list_g_bar = []
    for i in unique_indexes_with_fraction:
        ga = 0.0 + 0j
        zz = list_AtomicName[i]
        ff = list_fraction[i]
        if do_not_prototype:
            ga_item = numpy.exp(2j * numpy.pi * (hh * list_x[i] + kk * list_y[i] + ll * list_z[i]))
            ga += ga_item
        else:
            for j in range(number_of_atoms):
                if list_AtomicName[j] == zz and list_fraction[j] == ff:
                    ga_item = numpy.exp(2j * numpy.pi * (hh * list_x[j] + kk * list_y[j] + ll * list_z[j]))
                    ga += ga_item

        txt += "(%g,%g) \n" % (ga.real, ga.imag)
        txt += "(%g,%g) \n" % (ga.real, -ga.imag)
        list_g.append(ga)
        list_g_bar.append(ga.conjugate())
    output_dictionary["G"] = list_g
    output_dictionary["G_BAR"] = list_g_bar

    #
    # F0 part
    #
    txt += "# for each type of element-site, the number of f0 coefficients followed by them\n"
    for f0coeffs_item in f0coeffs:
        txt += "%d " % len(f0coeffs_item)
        for cc in f0coeffs_item:
            txt += "%g " % cc
        txt += "\n"
    output_dictionary["f0coeff"] = f0coeffs


    # X.J. Yu, use ceil to round up, otherwise we may get actual max energy less than emax
    npoint = int(numpy.ceil(((emax - emin) / estep + 1)))
    txt += "# The number of energy points NPOINT: \n"
    txt += ("%i \n") % npoint
    output_dictionary["npoint"] = npoint

    txt += "# for each energy point, energy, F1(1),F2(1),...,F1(nbatom),F2(nbatom)\n"
    list_energy = []
    out_f1 = numpy.zeros((len(unique_indexes_with_fraction), npoint), dtype=float)
    out_f2 = numpy.zeros((len(unique_indexes_with_fraction), npoint), dtype=float)
    out_fcompton = numpy.zeros((len(unique_indexes_with_fraction), npoint), dtype=complex)

    for i in range(npoint):
        energy = (emin + estep * i)
        txt += ("%20.11e \n") % (energy)
        list_energy.append(energy)

        for j,jj in enumerate(unique_indexes_with_fraction):
            f1a = xraylib.Fi(list_Zatom[jj], energy * 1e-3)
            f2a = -xraylib.Fii(list_Zatom[jj], energy * 1e-3)  # TODO: check the sign!!
            txt += (" %20.11e %20.11e 1.000 \n") % (f1a, f2a)
            out_f1[j, i] = f1a
            out_f2[j, i] = f2a
            out_fcompton[j, i] = 1.0

    output_dictionary["energy"] = list_energy
    output_dictionary["f1"] = out_f1
    output_dictionary["f2"] = out_f2
    output_dictionary["fcompton"] = out_fcompton

    if fileout != None:
        with open(fileout, "w") as f:
            f.write(txt)
        if verbose: print("File written to disk: %s" % fileout)

    return output_dictionary

def TemperFactor(sinTheta_lambda,anisos,Miller={'h':1,'k':1,'l':1},cell={'a':23.44,'b':23.44,'c':23.44},n=1936):
    '''
    #+
    # Singapore Synchrotron Light Source (SSLS)
    # :Author: X.J. Yu, slsyxj@nus.edu.sg
    # :Name:  TemperFactor
    # :Purpose: Calculation isotropic & anisotropic temerature factors
    # :Input:
    #     Miller: Miller indices
    #     cell:  dictionary of lattice [a,b,c] in units of Angstrom
    #     sinTheta_lambda: Sin(theta)/lambda, lambda in units of Angstrom
    #     n: number of atomic sites
    #     anisos: array of dictionary containing anisotropic coefficients
    #     Out: output results in a 2-elements list: [[sotropic],[anisotropic]]
    #-
    '''
    #0: isotropic, 1: anisotropic temerature factors
    results = numpy.zeros([2,n])
    for i,aniso in enumerate(anisos):
        s = aniso['start']-1
        e = aniso['end']
        if aniso['beta11'] >= 1:
            #if beta11>=1, then beta22 is Beq, the other fields are unused
            #if Beq specified, anisotropic temperature factor same as isotropic
            Beq = aniso['beta22']
            results[1,s:e] = numpy.exp(-sinTheta_lambda*sinTheta_lambda*Beq)
        else:
            Beq = 4.0/3.0*( aniso['beta11']*cell['a']*cell['a']+aniso['beta22']*cell['b']*cell['b']+ \
                aniso['beta33']*cell['c']*cell['c'] )
            results[1,s:e] = numpy.exp(-(aniso['beta11']*Miller['h']*Miller['h'] + \
                  aniso['beta22']*Miller['k']*Miller['k'] + aniso['beta33']*Miller['l']*Miller['l'] + \
                  2.0*Miller['h']*Miller['k']*aniso['beta12'] + 2.0*Miller['h']*Miller['l']*aniso['beta13'] + 2.0*Miller['k']*Miller['l']*aniso['beta23']))
        results[0,s:e] = numpy.exp(-sinTheta_lambda*sinTheta_lambda*Beq)

    return results



#################################################################################
# f0 tools
#################################################################################
def f0_calc_dabax(
    MAT_FLAG,
    DESCRIPTOR,
    GRIDSTART,
    GRIDEND,
    GRIDN,
    FILE_NAME="",
    charge=0.0,
    dabax_repository="",
    ):

    qscale = numpy.linspace(GRIDSTART, GRIDEND, GRIDN)

    f0 = numpy.zeros_like(qscale)

    if MAT_FLAG == 0: # element
        descriptor = DESCRIPTOR
        # for i,iqscale in enumerate(qscale):
        Z = xraylib.SymbolToAtomicNumber(descriptor)
        coeffs = f0_with_fractional_charge(Z, charge=charge,
                                          filename="f0_InterTables.dat",
                                          dabax_repository=dabax_repository)
        f0 = calculate_f0_from_f0coeff(coeffs, qscale)
    elif MAT_FLAG == 1: # formula
        tmp = CompoundParser(DESCRIPTOR, dabax_repository=dabax_repository)
        zetas = tmp["Elements"]
        multiplicity = tmp["nAtoms"]
        for j,jz in enumerate(zetas):
            coeffs = f0_with_fractional_charge(jz, charge=charge,
                                               filename="f0_InterTables.dat",
                                               dabax_repository=dabax_repository)
            f0 += multiplicity[j] * calculate_f0_from_f0coeff(coeffs, qscale)
    elif MAT_FLAG == 2: # nist

        Zarray = xraylib.GetCompoundDataNISTByName(DESCRIPTOR)
        zetas = Zarray["Elements"]
        fractions = numpy.array(Zarray["massFractions"])

        multiplicity = []
        for i in range(fractions.size):
            multiplicity.append( fractions[i] / xraylib.AtomicWeight(zetas[i]) )

        multiplicity = numpy.array(multiplicity)
        multiplicity /= multiplicity.min()

        atwt = 0.0
        for i in range(fractions.size):
            atwt += multiplicity[i] * xraylib.AtomicWeight(zetas[i])

        print("f0_calc - nist: ")
        print("    Descriptor: ", DESCRIPTOR)
        print("    Zs: ", zetas)
        print("    n: ", multiplicity)
        print("    atomic weight: ", atwt)

        for j, jz in enumerate(zetas):
            coeffs = f0_with_fractional_charge(jz, charge=charge,
                                               filename="f0_InterTables.dat",
                                               dabax_repository=dabax_repository)
            f0 += multiplicity[j] * calculate_f0_from_f0coeff(coeffs, qscale)

    else:
        raise Exception("Not implemented")

    if FILE_NAME != "":
        with open(FILE_NAME, "w") as file:
            try:
                file.write("#F %s\n"%FILE_NAME)
                file.write("\n#S 1 xoppy f0 results\n")
                file.write("#N 2\n")
                file.write("#L  q=sin(theta)/lambda [A^-1]  f0 [electron units]\n")
                for j in range(qscale.size):
                    # file.write("%19.12e  "%energy[j])
                    file.write("%19.12e  %19.12e\n"%(qscale[j],f0[j]))
                file.close()
                print("File written to disk: %s \n"%FILE_NAME)
            except:
                raise Exception("f0: The data could not be dumped onto the specified file!\n")
    #
    # return
    #
    return {"application":"xoppy","name":"f0","data":numpy.vstack((qscale,f0)),"labels":["q=sin(theta)/lambda [A^-1]","f0 [electron units]"]}


#
# test routines
#
def check_temperature_factor():
    #
    # anisotropy
    #

    cryst = Crystal_GetCrystal(filename='Crystals.dat', entry_name='YB66', dabax_repository=dabax_repository)


    if 'Aniso' in cryst.keys() and cryst['Aniso'][0]['start']>0:    #most crystals have no Anisotropic input
        print(">>> is Anisotropic", len(cryst["Aniso"]))

    for ele in cryst['Aniso']:
        print(ele)


    # Consider anisotropic temperature factor
    # X.J. Yu, slsyxj@nus.edu.sg
    # A dummy dictionary Aniso with start =0 if no aniso temperature factor input
    # start
    hh = 4
    kk = 0
    ll = 0
    dspacing = bragg_metrictensor(cryst['a'], cryst['b'], cryst['c'], cryst['alpha'], cryst['beta'], cryst['gamma'], HKL=[hh, kk, ll])
    dspacing *= 1e-8 # in cm


    if 'Aniso' in cryst.keys() and cryst['Aniso'][0]['start']>0:    #most crystals have no Anisotropic input
        print(">>> is Anisotropic")
        TFac = TemperFactor( 1.0/(2.0*dspacing*1e8),cryst['Aniso'],Miller={'h':hh,'k':kk,'l':ll}, \
            cell={'a':cryst['a'],'b':cryst['b'],'c':cryst['c']},
                             # n=len(atom),
                             n = cryst["n_atom"]
                             )
        B_TFac = 1
    else:
        B_TFac = 0

    print("TFac: ", TFac, cryst["n_atom"], len(TFac[0]), len(TFac[1]))


def check_structure_factor(descriptor="Si", hh=1, kk=1, ll=1, energy=8000,
                           do_assert=True, models=[1,1,1],
                           dabax_repository=""):

    from xoppy_xraylib_util import bragg_calc, crystal_fh


    os.system("rm -f xcrystal.bra xcrystal_0.bra xcrystal_1.bra xcrystal_2.bra")
    #
    # installed xoppy
    #
    if models[0]:
        dic0a = bragg_calc_old(descriptor=descriptor, hh=hh, kk=kk, ll=ll, temper=1.0,
                               emin=energy-100, emax=energy+100, estep=5.0,
                               fileout="xcrystal.bra")
        os.system("cp xcrystal.bra xcrystal_0.bra")
        dic0b = crystal_fh(dic0a, energy)

    #
    # new xoppy
    #
    if models[1]:
        dic1a = bragg_calc(descriptor=descriptor, hh=hh, kk=kk, ll=ll, temper=1.0,
                            emin=energy-100, emax=energy+100, estep=5.0,
                            fileout="xcrystal.bra")
        os.system("cp xcrystal.bra xcrystal_1.bra")
        dic1b = crystal_fh(dic1a, energy)

    #
    # Xiaojiang
    #
    if descriptor == "YB66":
        ANISO_SEL = 2
        do_not_prototype = 0
    else:
        ANISO_SEL = 0
        do_not_prototype = 0

    if models[2]:
        dic2a = bragg_calc2(descriptor=descriptor, hh=hh, kk=kk, ll=ll, temper=1.0,
                            emin=energy-100, emax=energy+100, estep=5.0, ANISO_SEL=ANISO_SEL,
                            do_not_prototype=do_not_prototype,sourceCryst=2,
                            fileout="xcrystal.bra", verbose=False, dabax_repository=dabax_repository)
        os.system("cp xcrystal.bra xcrystal_2.bra")

        dic2b = crystal_fh(dic2a, energy)


    if False:
        print(dic2b["info"])

        if models[0]: print("KEYS dict0a: ", dic0a.keys())
        if models[1]: print("KEYS dict1a: ", dic1a.keys())
        if models[2]: print("KEYS dict2a: ", dic2a.keys())
        if models[0]: print("KEYS dict1b: ", dic0b.keys())
        if models[1]: print("KEYS dict1b: ", dic1b.keys())
        if models[2]: print("KEYS dict2b: ", dic2b.keys())


        # if models[1] and models[2]:
        #     print(">>> COMPARING RESULT OF bragg_calc NEW  -  XIAOJIANG")
        #     for key in dic1a.keys():
        #         if key != "info":
        #             print(">>>", key, "\n   ", dic1a[key], "\n   ", dic2a[key])


        print("For Si 111 at % eV: " % energy)

        print("F0:")
        if models[0]: print(dic0b['F_0'])
        if models[1]: print(dic1b['F_0'])
        if models[2]: print(dic2b['F_0'])

        print("STRUCT:")
        if models[0]: print(dic0b['STRUCT'])
        if models[1]: print(dic1b['STRUCT'])
        if models[2]: print(dic2b['STRUCT'])


    if do_assert:
        if models[1] and models[2]:
            assert (numpy.abs(dic2b['STRUCT'] - dic1b['STRUCT']) < 1e-5)
            assert (numpy.abs(dic2b['FH'] - dic1b['FH']) < 1e-5)
            assert (numpy.abs(dic2b['FH_BAR'] - dic1b['FH_BAR']) < 1e-5)
            assert (numpy.abs(dic2b['F_0'] - dic1b['F_0']) < 1)
        if descriptor == "YB66":
            values_from = 0 # 0=Xiaojiang, 1=srio
            if values_from == 1:
                assert (numpy.abs(dic2b['STRUCT'] -  (565.7225232608029 + 35.9668881704435j))  < 1e-2)
                assert (numpy.abs(dic2b['FH'] -     (565.7225232608029 + 35.966888170443404j)) < 1e-2)
                assert (numpy.abs(dic2b['FH_BAR'] - (565.7225232608029 + 35.96688817044359j))  < 1e-2)
                assert (numpy.abs(dic2b['F_0'] -    (8846.406209552279 + 56.12593721027547j))  < 0.3)
            else:
                use_Atomic_name = True
                if not use_Atomic_name:
                    assert (numpy.abs(dic2b['STRUCT'] -  (563.4529619470779+35.8256810337139j))  < 1e-2)
                    assert (numpy.abs(dic2b['FH'] -      (563.4529619470779+35.82568103371383j)) < 1e-2)
                    assert (numpy.abs(dic2b['FH_BAR'] -  (563.4529619470779+35.82568103371397j))  < 1e-2)
                    assert (numpy.abs(dic2b['F_0'] -     (8848.638071350899+56.12049122626621j))  < 0.3)
                else:
                    assert (numpy.abs(dic2b['STRUCT'] -  (563.3428927313457+35.82567872885292j))  < 1e-2)
                    assert (numpy.abs(dic2b['FH'] -      (563.3428924224959+35.82568727860286j)) < 1e-2)
                    assert (numpy.abs(dic2b['FH_BAR'] -  (563.3428930401956+35.825670179102985j))  < 1e-2)
                    assert (numpy.abs(dic2b['F_0'] -     (8848.637396764192+56.12049119444192j))  < 0.3)

    return dic2b['STRUCT']

if __name__ == "__main__":
    import os


    # redefine the default server at ESRF because default server name is different outside and inside ESRF
    import socket
    dabax_repository = "http://ftp.esrf.fr/pub/scisoft/DabaxFiles/"
    if socket.getfqdn().find("esrf") >= 0:
        # dabax_repository = "http://ftp.esrf.fr/pub/scisoft/DabaxFiles/"
        dabax_repository = "/scisoft/DABAX/data"


    #
    # crystal
    #
    from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_calc as bragg_calc_old

    # test temperature
    # check_temperature_factor()

    # test Si
    check_structure_factor(descriptor="Si", hh=1, kk=1, ll=1, energy=8000,
                           dabax_repository=dabax_repository)

    # test Muscovite
    check_structure_factor(descriptor="Muscovite", hh=1, kk=1, ll=1, energy=8000, do_assert=1, models=[0,1,1],
                           dabax_repository=dabax_repository)

    # Test YB66
    check_structure_factor(descriptor="YB66", hh=4, kk=0, ll=0, energy=8040.0, do_assert=1, models=[0,0,1],
                           dabax_repository=dabax_repository)

    #
    # f0
    #
    if True:
        from orangecontrib.xoppy.util.xoppy_xraylib_util import f0_calc

        Si_xrl = f0_calc      (0, "Si", 0, 6, 100)
        Si_dbx = f0_calc_dabax(0, "Si", 0, 6, 100, dabax_repository=dabax_repository)

        H2O_xrl = f0_calc      (1, "H2O", 0, 6, 100)
        H2O_dbx = f0_calc_dabax(1, "H2O", 0, 6, 100, dabax_repository=dabax_repository)

        H2O_xrl = f0_calc      (2, "Water, Liquid", 0, 6, 100)
        H2O_dbx = f0_calc_dabax(2, "Water, Liquid", 0, 6, 100, dabax_repository=dabax_repository)


        from srxraylib.plot.gol import plot
        plot(Si_xrl["data"][0,:],Si_xrl["data"][1,:],
             Si_dbx["data"][0,:],Si_dbx["data"][1,:],
             H2O_xrl["data"][0, :], H2O_xrl["data"][1, :],
             H2O_dbx["data"][0, :], H2O_dbx["data"][1, :],
             linestyle=[None,'',None,''],
             marker=[None,'+',None,'+'],
             color=['r','r','b','b'],
             legend=['Si xraylib','Si dabax','H2O xraylib','H2O dabax'])

