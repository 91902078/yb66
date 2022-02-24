#
# this is the debugged code that goes in xoppy_xraylib_util!!!!!!!!!!!!!
#
# It corrects:
# differentiates element-sites with different fraction in bragg_calc
# uses the occupation fraction in crystal_fh
#

import xraylib
import numpy
import scipy.constants as codata

from orangecontrib.xoppy.util.xoppy_xraylib_util import f0_xop

#-------------------------------------------------------------------------
toangstroms = codata.h * codata.c / codata.e * 1e10

def bragg_calc(descriptor="Si",hh=1,kk=1,ll=1,temper=1.0,emin=5000.0,emax=15000.0,estep=100.0,fileout=None):
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
    :param fileout: name for the output file (default=None, no output file)
    :return: a dictionary with all ingredients of the structure factor.
    """

    output_dictionary = {}

    codata_e2_mc2 = codata.e**2 / codata.m_e / codata.c**2 / (4*numpy.pi*codata.epsilon_0)  # in m

    # f = open(fileout,'w')

    txt = ""
    txt += "# Bragg version, Data file type\n"
    txt += "2.4 1\n"

    cryst = xraylib.Crystal_GetCrystal(descriptor)
    volume = cryst['volume']

    #test crystal data - not needed
    itest = 0
    if itest:

        print ("  Unit cell dimensions are %f %f %f" % (cryst['a'],cryst['b'],cryst['c']))
        print ("  Unit cell angles are %f %f %f" % (cryst['alpha'],cryst['beta'],cryst['gamma']))
        print ("  Unit cell volume is %f A^3" % volume )
        print ("  Atoms at:")
        print ("     Z  fraction    X        Y        Z")
        for i in range(cryst['n_atom']):
            atom =  cryst['atom'][i]
            print ("    %3i %f %f %f %f" % (atom['Zatom'], atom['fraction'], atom['x'], atom['y'], atom['z']) )
        print ("  ")

    volume = volume*1e-8*1e-8*1e-8 # in cm^3
    dspacing = xraylib.Crystal_dSpacing(cryst, hh, kk, ll)
    rn = (1e0/volume)*(codata_e2_mc2*1e2)
    dspacing *= 1e-8 # in cm

    txt += "# RN = (e^2/(m c^2))/V) [cm^-2], d spacing [cm]\n"
    txt += "%e %e \n" % (rn , dspacing)

    output_dictionary["rn"] = rn
    output_dictionary["dspacing"] = dspacing

    atom = cryst['atom']
    list_Zatom = [ atom[i]['Zatom'] for i in range(len(atom))]
    number_of_atoms = len(list_Zatom)
    list_fraction = [ atom[i]['fraction'] for i in range(len(atom))]
    try:
        list_charge = [atom[i]['charge'] for i in range(len(atom))]
    except:
        list_charge = [0.0] * number_of_atoms
    list_x = [ atom[i]['x'] for i in range(len(atom))]
    list_y = [ atom[i]['y'] for i in range(len(atom))]
    list_z = [ atom[i]['z'] for i in range(len(atom))]

    # creates an is that contains Z, occupation and charge, that will
    # define the different sites.
    IDs = []
    number_of_atoms = len(list_Zatom)
    for i in range(number_of_atoms):
        IDs.append("Z:%2d-F:%g-C:%g" % (list_Zatom[i],list_fraction[i], list_charge[i]))

    unique_indexes = [IDs.index(x) for x in set(IDs)]

    unique_Zatom = [list_Zatom[i] for i in unique_indexes]
    unique_charge = [list_charge[i] for i in unique_indexes]
    unique_scattering_electrons = []
    for i, Zi in enumerate(unique_Zatom):
        unique_scattering_electrons.append(Zi - unique_charge[i])

    nbatom = (len(unique_Zatom))
    txt += "# Number of different element-sites in unit cell NBATOM:\n%d \n" % nbatom
    output_dictionary["nbatom"] = nbatom

    txt += "# for each element-site, the number of scattering electrons (Z_i - charge_i)\n"
    for i in unique_Zatom:
        txt += "%d "%i
    txt += "\n"
    output_dictionary["atnum"] = list(unique_scattering_electrons)

    txt += "# for each element-site, the occupation factor\n"
    unique_fraction = []
    for i in range(len(unique_indexes)):
        unique_fraction.append(list_fraction[unique_indexes[i]])
        txt += "%g "%(unique_fraction[i])
    txt += "\n"
    output_dictionary["fraction"] = unique_fraction


    txt += "# for each element-site, the temperature factor\n" # temperature parameter
    list_temper = []
    for i in range(len(unique_indexes)):
        txt += "%5.3f "%temper
        list_temper.append(temper)
    txt += "\n"
    output_dictionary["temper"] = list_temper

    #
    # Geometrical part of structure factor:  G and G_BAR
    #
    txt += "# for each type of element-site, COOR_NR=G_0\n"
    list_multiplicity = []
    for i in range(len(unique_indexes)):
        id = IDs[unique_indexes[i]]
        txt += "%d "%IDs.count(id)
        list_multiplicity.append(IDs.count(id))
    txt += "\n"
    output_dictionary["G_0"] = list_multiplicity

    txt += "# for each type of element-site, G and G_BAR (both complex)\n"
    list_g = []
    list_g_bar = []
    for i in range(len(unique_indexes)):
        id = IDs[unique_indexes[i]]
        ga = 0.0 + 0j
        for i,zz in enumerate(IDs):
            if zz == id:
                ga += numpy.exp(2j*numpy.pi*(hh*list_x[i]+kk*list_y[i]+ll*list_z[i]))
        txt += "(%g,%g) \n"%(ga.real,ga.imag)
        txt += "(%g,%g) \n"%(ga.real,-ga.imag)
        list_g.append(ga)
        list_g_bar.append(ga.conjugate())
    output_dictionary["G"] = list_g
    output_dictionary["G_BAR"] = list_g_bar

    #
    # F0 part
    #
    txt += "# for each type of element-site, the number of f0 coefficients followed by them\n"
    list_f0 = []
    for i in range(len(unique_indexes)):
        zeta = list_Zatom[unique_indexes[i]]
        tmp = f0_xop(zeta)
        txt += ("11 "+"%g "*11+"\n")%(tuple(tmp))
        list_f0.append(tmp.tolist())
    output_dictionary["f0coeff"] = list_f0


    npoint  = int( (emax - emin)/estep + 1 )
    txt += "# The number of energy points NPOINT: \n"
    txt +=  ("%i \n") % npoint
    output_dictionary["npoint"] = npoint
    txt += "# for each energy point, energy, F1(1),F2(1),...,F1(nbatom),F2(nbatom)\n"
    list_energy = []
    out_f1 = numpy.zeros( (len(unique_indexes),npoint), dtype=float)
    out_f2 = numpy.zeros( (len(unique_indexes),npoint), dtype=float)
    out_fcompton = numpy.zeros( (len(unique_indexes),npoint), dtype=complex)
    for i in range(npoint):
        energy = (emin+estep*i)
        txt += ("%20.11e \n") % (energy)
        list_energy.append(energy)

        for j in range(len(unique_indexes)):
            zeta = list_Zatom[unique_indexes[j]]
            f1a = xraylib.Fi(int(zeta),energy*1e-3)
            f2a = -xraylib.Fii(int(zeta),energy*1e-3) # TODO: check the sign!!
            txt +=  (" %20.11e %20.11e 1.000 \n")%(f1a, f2a)
            out_f1[j,i] = f1a
            out_f2[j,i] = f2a
            out_fcompton[j,i] = 1.0

    output_dictionary["energy"] = list_energy
    output_dictionary["f1"] = out_f1
    output_dictionary["f2"] = out_f2
    output_dictionary["fcompton"] = out_fcompton

    if fileout != None:
        with open(fileout,"w") as f:
            f.write(txt)
            print("File written to disk: %s" % fileout)

    return output_dictionary


def crystal_fh(input_dictionary,phot_in,theta=None,forceratio=0):
    """

    :param input_dictionary: as resulting from bragg_calc()
    :param phot_in: photon energy in eV
    :param theta: incident angle (half of scattering angle) in rad
    :return: a dictionary with structure factor
    """

    # outfil    = input_dictionary["outfil"]
    # fract     = input_dictionary["fract"]
    rn        = input_dictionary["rn"]
    dspacing  = numpy.array(input_dictionary["dspacing"])
    nbatom    = numpy.array(input_dictionary["nbatom"])
    atnum     = numpy.array(input_dictionary["atnum"])
    temper    = numpy.array(input_dictionary["temper"])
    G_0       = numpy.array(input_dictionary["G_0"])
    G         = numpy.array(input_dictionary["G"])
    G_BAR     = numpy.array(input_dictionary["G_BAR"])
    f0coeff   = numpy.array(input_dictionary["f0coeff"])
    npoint    = numpy.array(input_dictionary["npoint"])
    energy    = numpy.array(input_dictionary["energy"])
    fp        = numpy.array(input_dictionary["f1"])
    fpp       = numpy.array(input_dictionary["f2"])
    fraction = numpy.array(input_dictionary["fraction"])



    phot_in = numpy.array(phot_in,dtype=float).reshape(-1)

    toangstroms = codata.h * codata.c / codata.e * 1e10


    itheta = numpy.zeros_like(phot_in)
    for i,phot in enumerate(phot_in):

        if theta is None:
            itheta[i] = numpy.arcsin(toangstroms*1e-8/phot/2/dspacing)
        else:
            itheta[i] = theta

        # print("energy= %g eV, theta = %15.13g deg"%(phot,itheta[i]*180/numpy.pi))
        if phot < energy[0] or phot > energy[-1]:
            raise Exception("Photon energy %g eV outside of valid limits [%g,%g]"%(phot,energy[0],energy[-1]))

        if forceratio == 0:
            ratio = numpy.sin(itheta[i]) / (toangstroms / phot)
        else:
            ratio = 1 / (2 * dspacing * 1e8)
        # print("Ratio: ",ratio)

        F0 = numpy.zeros(nbatom)
        for j in range(nbatom):
            icentral = int(f0coeff.shape[1]/2)
            F0[j] = f0coeff[j,icentral]
            for i in range(icentral):
                F0[j] += f0coeff[j,i] * numpy.exp(-1.0*f0coeff[j,i+icentral+1]*ratio**2)


        # ;C
        # ;C Interpolate for the atomic scattering factor.
        # ;C
        for j,ienergy in enumerate(energy):
            if ienergy > phot:
                break
        nener = j - 1


        F1 = numpy.zeros(nbatom,dtype=float)
        F2 = numpy.zeros(nbatom,dtype=float)
        F = numpy.zeros(nbatom,dtype=complex)

        for j in range(nbatom):
            F1[j] = fp[j,nener] + (fp[j,nener+1] - fp[j,nener]) * \
            (phot - energy[nener]) / (energy[nener+1] - energy[nener])
            F2[j] = fpp[j,nener] + (fpp[j,nener+1] - fpp[j,nener]) * \
            (phot - energy[nener]) / (energy[nener+1] - energy[nener])

        r_lam0 = toangstroms * 1e-8 / phot
        for j in range(nbatom):
            F[j] = F0[j] + F1[j] + 1j * F2[j]
            # print("F",F)


        F_0 = 0.0 + 0.0j
        FH = 0.0 + 0.0j
        FH_BAR = 0.0 + 0.0j
        FHr = 0.0 + 0.0j
        FHi = 0.0 + 0.0j
        FH_BARr = 0.0 + 0.0j
        FH_BARi = 0.0 + 0.0j


        TEMPER_AVE = 1.0
        for j in range(nbatom):
            FH  += fraction[j] * (G[j] *   F[j] * 1.0)
            FHr += fraction[j] * (G[j] * (F0[j] + F1[j])* 1.0)
            FHi += fraction[j] * (G[j] *  F2[j] * 1.0)
            F_0 += fraction[j] * (G_0[j] * ( atnum[j] + F1[j] + 1j * F2[j] ) * 1.0)
            TEMPER_AVE *= (temper[j])**(G_0[j]/(G_0.sum()))

            FH_BAR  += fraction[j] * ((G_BAR[j] * F[j] * 1.0))
            FH_BARr += fraction[j] * ((G_BAR[j] * (F0[j]  + F1[j]) *1.0))
            FH_BARi += fraction[j] * ((G_BAR[j] *  F2[j] * 1.0))
            # print("TEMPER_AVE: ",TEMPER_AVE)


        # ;C
        # ;C multiply by the average temperature factor
        # ;C


        FH      *= TEMPER_AVE
        FHr     *= TEMPER_AVE
        FHi     *= TEMPER_AVE
        FH_BAR  *= TEMPER_AVE
        FH_BARr *= TEMPER_AVE
        FH_BARi *= TEMPER_AVE

        STRUCT = numpy.sqrt(FH * FH_BAR)

        # ;C
        # ;C   PSI_CONJ = F*( note: PSI_HBAR is PSI at -H position and is
        # ;C   proportional to fh_bar but PSI_CONJ is complex conjugate os PSI_H)
        # ;C


        psi_over_f = rn * r_lam0**2 / numpy.pi
        psi_h      = rn * r_lam0**2 / numpy.pi * FH
        psi_hr     = rn * r_lam0**2 / numpy.pi * FHr
        psi_hi     = rn * r_lam0**2 / numpy.pi * FHi
        psi_hbar   = rn * r_lam0**2 / numpy.pi * FH_BAR
        psi_hbarr  = rn * r_lam0**2 / numpy.pi * FH_BARr
        psi_hbari  = rn * r_lam0**2 / numpy.pi * FH_BARi
        psi_0      = rn * r_lam0**2 / numpy.pi * F_0
        psi_conj   = rn * r_lam0**2 / numpy.pi * FH.conjugate()

        # ;
        # ; Darwin width
        # ;
        # print(rn,r_lam0,STRUCT,itheta)
        ssvar = rn * (r_lam0**2) * STRUCT / numpy.pi / numpy.sin(2.0*itheta)
        spvar = ssvar * numpy.abs((numpy.cos(2.0*itheta)))
        ssr = ssvar.real
        spr = spvar.real

        # ;C
        # ;C computes refractive index.
        # ;C ([3.171] of Zachariasen's book)
        # ;C
        REFRAC = (1.0+0j) - r_lam0**2 * rn * F_0 / 2/ numpy.pi
        DELTA_REF = 1.0 - REFRAC.real
        ABSORP = 4.0 * numpy.pi * (-REFRAC.imag) / r_lam0


        txt = ""
        txt += '\n******************************************************'
        txt += '\n       at energy    = '+repr(phot)+' eV'
        txt += '\n                    = '+repr(r_lam0*1e8)+' Angstroms'
        txt += '\n       and at angle = '+repr(itheta*180.0/numpy.pi)+' degrees'
        txt += '\n                    = '+repr(itheta)+' rads'
        txt += '\n******************************************************'

        for j in range(nbatom):
            txt += '\n  '
            txt += '\nFor atom '+repr(j+1)+':'
            txt += '\n       fo + fp+ i fpp = '
            txt += '\n        '+repr(F0[j])+' + '+ repr(F1[j].real)+' + i'+ repr(F2[j])+" ="
            txt += '\n        '+repr(F0[j] + F1[j] + 1j * F2[j])
            txt += '\n       Z = '+repr(atnum[j])
            txt += '\n       Temperature factor = '+repr(temper[j])
        txt += '\n  '
        txt += '\n Structure factor F(0,0,0) = '+repr(F_0)
        txt += '\n Structure factor FH = '      +repr(FH)
        txt += '\n Structure factor FH_BAR = '  +repr(FH_BAR)
        txt += '\n Structure factor F(h,k,l) = '+repr(STRUCT)
        txt += '\n  '
        txt += '\n Psi_0  = '   +repr(psi_0)
        txt += '\n Psi_H  = '   +repr(psi_h)
        txt += '\n Psi_HBar  = '+repr(psi_hbar)
        txt += '\n  '
        txt += '\n Psi_H(real) Real and Imaginary parts = '   + repr(psi_hr)
        txt += '\n Psi_H(real) Modulus  = '                   + repr(numpy.abs(psi_hr))
        txt += '\n Psi_H(imag) Real and Imaginary parts = '   + repr(psi_hi)
        txt += '\n Psi_H(imag) Modulus  = '                   + repr(abs(psi_hi))
        txt += '\n Psi_HBar(real) Real and Imaginary parts = '+ repr(psi_hbarr)
        txt += '\n Psi_HBar(real) Modulus  = '                + repr(abs(psi_hbarr))
        txt += '\n Psi_HBar(imag) Real and Imaginary parts = '+ repr(psi_hbari)
        txt += '\n Psi_HBar(imag) Modulus  = '                + repr(abs(psi_hbari))
        txt += '\n  '
        txt += '\n Psi/F factor = '                           + repr(psi_over_f)
        txt += '\n  '
        txt += '\n Average Temperature factor = '             + repr(TEMPER_AVE)
        txt += '\n Refraction index = 1 - delta - i*beta'
        txt += '\n            delta = '                       + repr(DELTA_REF)
        txt += '\n             beta = '                       + repr(1.0e0*REFRAC.imag)
        txt += '\n Absorption coeff = '                       + repr(ABSORP)+' cm^-1'
        txt += '\n  '
        txt += '\n e^2/(mc^2)/V = '                           + repr(rn)+' cm^-2'
        txt += '\n d-spacing = '                              + repr(dspacing*1.0e8)+' Angstroms'
        txt += '\n SIN(theta)/Lambda = '                      + repr(ratio)
        txt += '\n  '
        txt += '\n Darwin width for symmetric s-pol [microrad] = ' + repr(2.0e6*ssr)
        txt += '\n Darwin width for symmetric p-pol [microrad] = ' + repr(2.0e6*spr)

    return {"PHOT":phot, "WAVELENGTH":r_lam0*1e-2 ,"THETA":itheta, "F_0":F_0, "FH":FH, "FH_BAR":FH_BAR,
	        "STRUCT":STRUCT, "psi_0":psi_0, "psi_h":psi_h, "psi_hbar":psi_hbar,
        	"DELTA_REF":DELTA_REF, "REFRAC":REFRAC, "ABSORP":ABSORP, "RATIO":ratio,
        	"ssr":ssr, "spr":spr, "psi_over_f":psi_over_f, "info":txt}


if __name__ == "__main__":

    from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_calc as bragg_calc_old
    from orangecontrib.xoppy.util.xoppy_xraylib_util import crystal_fh as crystal_fh_old
    from crystal_util import bragg_calc2, crystal_fh2
    import os
    from run_diff_pat import run_diff_pat


    #
    #
    #

    for descriptor in ["Muscovite", "Si"]:
        #
        # old code muscovite
        #
        dic1a = bragg_calc_old(descriptor=descriptor,hh=1,kk=1,ll=1,temper=1.0,emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra")
        os.system("cp xcrystal.bra xcrystal.bra.old")

        dic1b = crystal_fh_old(dic1a,8000.0)

        run_diff_pat(
            MOSAIC=0,
            GEOMETRY=0,
            SCAN=2,
            UNIT=1,
            SCANFROM=25.0,
            SCANTO=100.0,
            SCANPOINTS=200,
            ENERGY=8040.0,
            ASYMMETRY_ANGLE=0.0,
            THICKNESS=0.7,
            MOSAIC_FWHM=0.1,
            RSAG=125.0,
            RMER=1290.0,
            ANISOTROPY=0,
            POISSON=0.22,
            CUT="2 -1 -1 ; 1 1 1 ; 0 0 0",
            FILECOMPLIANCE="mycompliance.dat")

        a1 = numpy.loadtxt("diff_pat.dat", skiprows=5)


        #
        # New code muscovite
        #

        dic2a = bragg_calc(descriptor=descriptor,hh=1,kk=1,ll=1,temper=1.0,emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra")
        print("KEYS: ",dic2a.keys())
        os.system("cp xcrystal.bra xcrystal.bra.new")

        dic2b = crystal_fh(dic2a,8000.0)
        print(dic2b["info"])
        print("KEYS: ",dic2b.keys())

        run_diff_pat(
            MOSAIC=0,
            GEOMETRY=0,
            SCAN=2,
            UNIT=1,
            SCANFROM=25.0,
            SCANTO=100.0,
            SCANPOINTS=200,
            ENERGY=8040.0,
            ASYMMETRY_ANGLE=0.0,
            THICKNESS=0.7,
            MOSAIC_FWHM=0.1,
            RSAG=125.0,
            RMER=1290.0,
            ANISOTROPY=0,
            POISSON=0.22,
            CUT="2 -1 -1 ; 1 1 1 ; 0 0 0",
            FILECOMPLIANCE="mycompliance.dat")

        a2 = numpy.loadtxt("diff_pat.dat", skiprows=5)

        #
        # Xiaojiang code
        #


        dic3a = bragg_calc2(descriptor=descriptor, hh=1, kk=1, ll=1, temper=1.0, emin=7900.0, emax=8100.0, estep=5.0,
                            fileout="xcrystal.bra")
        print("KEYS: ", dic3a.keys())
        os.system("cp xcrystal.bra xcrystal.bra.xj")


        dic3b = crystal_fh2(dic3a, 8000.0)
        print(dic3b["info"])
        print("KEYS: ", dic3b.keys())

        run_diff_pat(
            MOSAIC=0,
            GEOMETRY=0,
            SCAN=2,
            UNIT=1,
            SCANFROM=25.0,
            SCANTO=100.0,
            SCANPOINTS=200,
            ENERGY=8040.0,
            ASYMMETRY_ANGLE=0.0,
            THICKNESS=0.7,
            MOSAIC_FWHM=0.1,
            RSAG=125.0,
            RMER=1290.0,
            ANISOTROPY=0,
            POISSON=0.22,
            CUT="2 -1 -1 ; 1 1 1 ; 0 0 0",
            FILECOMPLIANCE="mycompliance.dat")

        a3 = numpy.loadtxt("diff_pat.dat", skiprows=5)





        for key in dic1b.keys():
            if key != "info":
                print(">>>", key,dic1b[key],dic2b[key],dic3b[key])
                # tmp = numpy.abs(dic1b[key] - dic2b[key])
                # if tmp.size == 1:
                #     assert (tmp < 1e-6)
                # else:
                #     assert(tmp.sum() < 1e-6)

        #
        # comparison
        #

        print(a1.shape, a2.shape, a3.shape)
        from srxraylib.plot.gol import plot, set_qt
        set_qt()

        plot(a1[:, 0], a1[:, -1],
             a2[:, 0], a2[:, -1],
             a3[:, 0], a3[:, -1],
             title="Crystal: " + descriptor,
             marker=[None, "o", "+"],
             legend=["OLD code (current xoppy)", "FIXED code (this file)", "XIAOJIANG code"])