import sys
import xraylib
import numpy
import os
import scipy.constants as codata
# X.J. Yu, slsyxj@nus.edu.sg
# from temperature_anisotropy import TemperFactor
from orangecontrib.xoppy.util.xoppy_xraylib_util import f0_xop
from dabax_util import crystal_parser, crystal_atnum
from dabax_util import get_f0_coeffs
from dabax_util import symbol_to_from_atomic_number
from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_metrictensor
import re
#-------------------------------------------------------------------------
toangstroms = codata.h * codata.c / codata.e * 1e10


def bragg_calc2(descriptor="YB66",hh=1,kk=1,ll=1,temper=1.0,emin=5000.0,emax=15000.0,estep=100.0,ANISO_SEL=0,fileout=None):
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
    
    #cryst = xraylib.Crystal_GetCrystal(descriptor)
    cryst = crystal_parser(entry_name=descriptor)
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

    alpha = cryst['alpha'] * numpy.pi/180
    beta =  cryst['beta']  * numpy.pi/180
    gamma = cryst['gamma'] * numpy.pi/180
    #dspacing = xraylib.Crystal_dSpacing(cryst, hh, kk, ll)
    dspacing = bragg_metrictensor(cryst['a'], cryst['b'], cryst['c'], cryst['alpha'], cryst['beta'], cryst['gamma'], HKL=[hh, kk, ll])
    dspacing *= 1e-8 # in cm
    volume = volume*1e-8*1e-8*1e-8 # in cm^3
    rn = (1e0/volume)*(codata_e2_mc2*1e2)

    txt += "# RN = (e^2/(m c^2))/V) [cm^-2], d spacing [cm]\n"
    txt += "%e %e \n" % (rn , dspacing)

    output_dictionary["rn"] = rn
    output_dictionary["dspacing"] = dspacing

    atom = cryst['atom']
    list_Zatom = [ atom[i]['Zatom'] for i in range(len(atom))]
    list_fraction = [ atom[i]['fraction'] for i in range(len(atom))]
    list_x = [ atom[i]['x'] for i in range(len(atom))]
    list_y = [ atom[i]['y'] for i in range(len(atom))]
    list_z = [ atom[i]['z'] for i in range(len(atom))]

    unique_Zatom = set(list_Zatom)
##  ------------ XJ.YU  Singapore Synchrotorn Light Source --------------------------
##  For backward compatible
    if len(atom[0]) >= 7:  #6 column + 1 AtomicName or +1 SeqNo (xraylib)
        f0coeffs = get_f0_coeffs(atom, list_Zatom)
        list_AtomicName = [ atom[i]['AtomicName'] for i in range(len(atom))]
        unique_AtomicName = list(sorted(set(list_AtomicName)))
    else:  #usually normal 5 column  
        """         cryst['Aniso']=[{'start':0}]
                for i in range(len(atom)):
                    atom[i]['AtomicName']=''
        """    
        list_AtomicName = ['']
        unique_AtomicName = ['']
    
    #unique_AtomicName has at least one empty string
    if unique_AtomicName[0] !='':
        #now unique_Zatom is changed from set to list, allow duplicate atomic number
        #because same atom at different sites may have different valences, i.e., O2-,O1.5-
        unique_Zatom=[]
        for z in unique_AtomicName:
            tmp = re.search('(^[a-zA-Z]*)',z)
            unique_Zatom.append(symbol_to_from_atomic_number(tmp.group(0)))
##  ------------ Singapore Synchrotorn Light Source ---------------------------------
    TmpCrystal = () # for diff_pat.exe
    if unique_AtomicName[0] !='':   #Complex crystal
        TmpCrystal = crystal_atnum(list_AtomicName, unique_AtomicName, unique_Zatom,list_fraction,f0coeffs)
    nbatom = (len(unique_Zatom))
    if unique_AtomicName[0] =='':
        txt += "# Number of different element-sites in unit cell NBATOM:\n%d \n" % nbatom
    else:    
        txt += "# Number of different element-sites in unit cell NBATOM:\n%d \n" % len(TmpCrystal[0])
    output_dictionary["nbatom"] = nbatom    # different with diff_pat for complex crystal

    txt += "# for each element-site, the atomic number\n"
    if unique_AtomicName[0] !='':   #Complex crystal
        for i in TmpCrystal[0]:
            #i = int(i + 0.5)        #round to integer value, diff_pat.exe not support float ATNUM
            #txt += "%d "%i
            txt += "%f "%i
    else:    #normal crystals
        for i in unique_Zatom:
            txt += "%d "%i
    txt += "\n"
    if len(TmpCrystal) > 0:
        output_dictionary["atnum"] = list(TmpCrystal[0])
    else:    
        output_dictionary["atnum"] = list(unique_Zatom)
    #XJ.YU  Singapore Synchrotorn Light Source
    output_dictionary["zcol"] = list(list_Zatom)
    output_dictionary["unique_AtomicName"] = list(unique_AtomicName)
    output_dictionary["list_AtomicName"] = list(list_AtomicName)
    
    #TODO: manage correctly fraction, the ones in non-representative atoms are ignored.
    txt += "# for each element-site, the occupation factor\n"
    unique_fraction = []
    if len(TmpCrystal) == 0:    #normal crystal
        for i in range(len(unique_Zatom)):
    #
    #commenuts By XJ.YU, xiaojiang@nus.edu.sg
    # for Muscovite crystal (KAl2(AlSi3)O10(OH)2), Al has two occupancy (1, 0.25), 5 atomic types
    # if using the number of unique_Zatom for unique_fraction, only first five numbers used in list_fraction
            unique_fraction.append(list_fraction[i])
            txt += "%g "%(unique_fraction[i])
    else:   #complex crystal with charge
         for z in TmpCrystal[1]:  #fractional
             unique_fraction.append(z)
             txt += "%g "%(z)
    txt += "\n"
    
# coment out by XJ.YU
#    output_dictionary["fraction"] = unique_fraction
#
# because even for same kind atom in different sites could have different occupancy,Like YB66, B1,B2,etc
# so keep the original fraction list
#
    output_dictionary["fraction"] = list_fraction   #not unique_fraction, full list

    txt += "# for each element-site, the temperature factor\n" # temperature parameter
    list_temper = []
    if len(TmpCrystal) > 0:    #complex crystal
        for i in TmpCrystal[1]:
            txt += "%5.3f "%temper      #for diff_pat.exe 
    for i in range(len(unique_Zatom)):
        if len(TmpCrystal) == 0:
            txt += "%5.3f "%temper
        list_temper.append(temper)
    txt += "\n"
    output_dictionary["temper"] = list_temper   #not necessary same with diff_pat
    
    #
    # Geometrical part of structure factor:  G and G_BAR
    #
    txt += "# for each type of element-site, COOR_NR=G_0\n"
    list_multiplicity = []
    #when there are duplicate atomic numbers in unique_Zatom, it is not correct anymore
    #should use unique_AtomicName, and list_AtomicName instead
    # commnent out: X.J. Yu, slsyxj@nus.edu.sg
    if unique_AtomicName[0] =='':
        for z in unique_Zatom:
            txt += "%d "%list_Zatom.count(z)
            list_multiplicity.append(list_Zatom.count(z))
    else:
        for z in unique_AtomicName:
        #    txt += "%d "%list_AtomicName.count(z)
            list_multiplicity.append(list_AtomicName.count(z))
        for z in TmpCrystal[2]:
            txt += "%d "%z
    txt += "\n"
    output_dictionary["G_0"] = list_multiplicity
    #
    # Consider anisotropic temperature factor
    # X.J. Yu, slsyxj@nus.edu.sg
    # A dummy dictionary Aniso with start =0 if no aniso temperature factor input
    # start
    if 'Aniso' in cryst.keys() and cryst['Aniso'][0]['start']>0:    #most crystals have no Anisotropic input
        TFac = TemperFactor( 1.0/(2.0*dspacing*1e8),cryst['Aniso'],Miller={'h':hh,'k':kk,'l':ll}, \
            cell={'a':cryst['a'],'b':cryst['b'],'c':cryst['c']},n=len(atom))
        B_TFac = 1
    else:
        B_TFac = 0
    # end
    #
    txt += "# for each type of element-site, G and G_BAR (both complex)\n"
    list_g = []
    list_g_bar = []
    ss =  numpy.array([numpy.exp(2j*numpy.pi*(hh*list_x[i]+kk*list_y[i]+ll*list_z[i]))*list_fraction[i] for i in range(len(list_x))])
    if B_TFac:
        TCoff = TFac[ANISO_SEL,:]
        ss = ss*TCoff           #multiple Anisotropic factor
    if unique_AtomicName[0] =='': #normal crystal
        for z in unique_Zatom:
            ga = numpy.sum(numpy.where(numpy.array(list_Zatom)==z,ss,0))
            txt += "(%g,%g) \n"%(ga.real,ga.imag)
            txt += "(%g,%g) \n"%(ga.real,-ga.imag)
            list_g.append(ga)
            list_g_bar.append(ga.conjugate())
    else:  #complex crystal
        for z in unique_AtomicName:
            ga = numpy.sum(numpy.where(numpy.array(list_AtomicName)==z,ss,0))
            list_g.append(ga)
            list_g_bar.append(ga.conjugate())
        ss = ss/list_fraction  #fraction handled in diff_pat.exe
        for i,z in enumerate(TmpCrystal[4]): #prepare G for xcraystal.bra, z is unique name
            s2 = TmpCrystal[1][i]       #fraction
            s3 = numpy.where(numpy.array(list_AtomicName)==z,list_fraction,0)  #get fraction from unique name
            ga = numpy.sum(numpy.where(s3==s2,ss,0))
            txt += "(%g,%g) \n"%(ga.real,ga.imag)
            txt += "(%g,%g) \n"%(ga.real,-ga.imag)
    output_dictionary["G"] = list_g
    output_dictionary["G_BAR"] = list_g_bar
    #
    # F0 part
    #
    txt += "# for each type of element-site, the number of f0 coefficients followed by them\n"
    list_f0 = []
    tmp_g={}        #add for diff_pat.exe, multiple different sites with same atom    
    for i,zeta in enumerate(unique_Zatom):
    #Comment: X.J. Yu, slsyxj@nus.edu.sg 
    #For complicate compound crystal, we use unique_AtomicName instead of atomic number to get f0
    #
        if unique_AtomicName[0] !='':   #with compound crystal name input
            tmp1 = re.search('(^[a-zA-Z]*)',unique_AtomicName[i])
            if tmp1.group(0) == unique_AtomicName[i]:   
            #Atom name only, without suffix for valence info (i.e., B, Y, O)
                tmp = f0_xop(zeta)
            elif unique_AtomicName[i] in f0coeffs:   #second notation, 'B-0.045548', etc
                tmp = f0coeffs[unique_AtomicName[i]]
            else:   #first notation 'B-.', etc
                tmp = f0_xop(0,unique_AtomicName[i])
        else:   # original notation only
            tmp = f0_xop(zeta)
        # print(("%g "*11)%(tmp.tolist()))
        nn = len(tmp)
        if len(TmpCrystal) == 0: #normal crystal
            txt += (str(nn) + " "+"%g "*nn+"\n")%(tuple(tmp))
        else:  #temporaty save here
            tmp_g[str(zeta)] =  tuple(tmp)
        # By XJ.Yu, return value already changed from array to list
        #list_f0.append(tmp.tolist())
        list_f0.append(tmp)
    if len(TmpCrystal) > 0: #for diff_pat.exe
        for zeta in TmpCrystal[3]: #receive the f0 for atom at each site
            nn = len(tmp_g[str(zeta)])
            #txt += ("11 "+"%g "*11+"\n")%(tmp_g[str(zeta)])        #not work for 9 column f0 
            txt += (str(nn) + " "+"%g "*nn+"\n")%(tmp_g[str(zeta)])
    output_dictionary["f0coeff"] = list_f0

    # f.write("# -----------------------------------------------\n")


    # zetas = numpy.array([atom[0]["Zatom"],atom[7]["Zatom"]])
    # X.J. Yu, use ceil to round up, otherwise we may get actual max energy less than emax
    npoint  = int(numpy.ceil(( (emax - emin)/estep + 1 )))      
    txt += "# The number of energy points NPOINT: \n"
    txt +=  ("%i \n") % npoint
    output_dictionary["npoint"] = npoint
    txt += "# for each energy point, energy, F1(1),F2(1),...,F1(nbatom),F2(nbatom)\n"
    list_energy = []
    out_f1 =        numpy.zeros( (len(unique_Zatom),npoint), dtype=float)
    out_f2 =        numpy.zeros( (len(unique_Zatom),npoint), dtype=float)
    out_fcompton =  numpy.zeros( (len(unique_Zatom),npoint), dtype=complex)
    for i in range(npoint):
        energy = (emin+estep*i)
        txt += ("%20.11e \n") % (energy)
        list_energy.append(energy)

        if len(TmpCrystal) > 0: #for diff_pat.exe
            tmp_g = TmpCrystal[3]
        else:
            tmp_g =  list(unique_Zatom)  
        for j,zeta in enumerate(unique_Zatom):
            f1a = xraylib.Fi(int(zeta),energy*1e-3)
            f2a = -xraylib.Fii(int(zeta),energy*1e-3) # TODO: check the sign!!
            for x in range(tmp_g.count(zeta)):  #treat different occupation for same atoms with a input line
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




def crystal_fh2(input_dictionary,phot_in,theta=None,forceratio=0):
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
#X.J. Yu, slsyxj@nus.edu.sg
    ZCOL     = numpy.array(input_dictionary["zcol"])
    FCOL     = numpy.array(input_dictionary["fraction"])
    UCOL     = numpy.array(input_dictionary["unique_AtomicName"])
    LCOL     = numpy.array(input_dictionary["list_AtomicName"])
#---------------------------------------------------------
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
        #X.J. Yu, slsyxj@nus.edu.sg
        F000 = numpy.zeros(nbatom)
        for j in range(nbatom):
            #allow f0coeff contains 9 and 11 columns at the same time
            #icentral = int(f0coeff.shape[1]/2)
            #F0[j] = f0coeff[j,icentral]
            f0_item = numpy.array(f0coeff[j])
            icentral = int(len(f0_item)/2)
            F0[j] = f0_item[icentral]
            F000[j] = F0[j] #X.J. Yu, slsyxj@nus.edu.sg
            for i in range(icentral):
                #F0[j] += f0coeff[j,i] * numpy.exp(-1.0*f0coeff[j,i+icentral+1]*ratio**2)
                #F000[j] += f0coeff[j,i]  #actual number of electrons carried by each atom, X.J. Yu, slsyxj@nus.edu.sg
                F0[j] += f0_item[i] * numpy.exp(-1.0*f0_item[i+icentral+1]*ratio**2)
                F000[j] += f0_item[i]  #actual number of electrons carried by each atom, X.J. Yu, slsyxj@nus.edu.sg

            # print("F0: ",F0,xraylib.FF_Rayl(int(atnum[j]),ratio))


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
            # print("F1,F2",F1,F2)

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
        CI =  0.0 + 1.0j

        TEMPER_AVE = 1.0
        #X.J. Yu, slsyxj@nus.edu.sg
        #Occupancy for FH already included in G in Bragg_Calc function
        BOOL_UCOL = UCOL[0]==''
        for j in range(nbatom):
            FH  += G[j] *   F[j] * 1.0
            FHr += G[j] * (F0[j] + F1[j])* 1.0
            FHi += G[j] *  F2[j] * 1.0
#charged atom, the number of electrons not equal to atum anymore,while
# it is euqal to F000, and notably, fractial occupancy need consideration here
# occupancy till now, only consider in calculation of G, and G_BAR in bragg_calc
#comment out: X.J. Yu, slsyxj@nus.edu.sg
#
#            F_0 += G_0[j] * ( atnum[j] + F1[j] + 1j * F2[j] ) * 1.0
#
            FN = F000[j] + F1[j] + CI * F2[j]
            if BOOL_UCOL:   #normal crystal
                F_0 += FN*numpy.sum( numpy.where(ZCOL==atnum[j],FCOL,0.0))
            else:   #complicate compound crystals
                F_0 += FN*numpy.sum( numpy.where(LCOL==UCOL[j],FCOL,0.0))
                
            TEMPER_AVE *= (temper[j])**(G_0[j]/(G_0.sum()))

            FH_BAR  += (G_BAR[j] * F[j] * 1.0)
            FH_BARr += (G_BAR[j] * (F0[j]  + F1[j]) *1.0)
            FH_BARi += (G_BAR[j] *  F2[j] * 1.0)
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

        THETA_B =r_lam0/(1-(DELTA_REF/numpy.sin(itheta)**2))/2.0/dspacing
        THETA_B = numpy.arcsin(THETA_B)

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

    return {"PHOT":phot, "WAVELENGTH":r_lam0*1e-2 ,"THETA":itheta,"THETAcor":THETA_B, "F_0":F_0, "FH":FH, "FH_BAR":FH_BAR,
	        "STRUCT":STRUCT, "psi_0":psi_0, "psi_h":psi_h, "psi_hbar":psi_hbar,
        	"DELTA_REF":DELTA_REF, "REFRAC":REFRAC, "ABSORP":ABSORP, "RATIO":ratio,
        	"ssr":ssr, "spr":spr, "psi_over_f":psi_over_f, "info":txt}

def TemperFactor(sinTheta_lambda,anisos,Miller={'h':1,'k':1,'l':1},cell={'a':23.44,'b':23.44,'c':23.44},n=1936):
    '''
    #+
    # Singapore Synchrotron Light Source (SSLS)
    # :Author: X.J. Yu, slsyxj@nus.edu.sg
    # :Name:  TemperFactor
    # :Purpose: Calculation isotropic & anisotropic temerature factors
    # :Input:
    #     Miller: Miller indice
    #     cell:  dictionary of lattice [a,b,c] in units of Aangstrom
    #     sinTheta_lambda: Sin(theta)/lambda, lambda in units of Aangstrom
    #     n: number of atomic sites
    #     anisos: array of dicionary contain anisotropic coefficients
    #     Out: output results, column 0: isotropic, column 1: anisotropic
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