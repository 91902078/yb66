import numpy

from xoppylib.crystals.tools import TemperFactor

import os
# import scipy.constants as codata

# from orangecontrib.xoppy.util.xoppy_xraylib_util import f0_xop
# from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_metrictensor

# from dabax_util import calculate_f0_from_f0coeff, f0_with_fractional_charge
# from dabax_util import Crystal_GetCrystal
# from dabax_util import atomic_symbols_dabax # __symbol_to_from_atomic_number
# from dabax_util import f0_with_fractional_charge
# from dabax_util import CompoundParser
# from xoppy_xraylib_util import load_bragg_preprocessor_file
# to be removed...  TODO: move the f1 f2 routines from xraylib to dabax.
# import xraylib


#
# test routines
#
def check_temperature_factor(material_constants_library):
    #
    # anisotropy
    #

    cryst = material_constants_library.Crystal_GetCrystal('YB66')

    if 'Aniso' in cryst.keys() and cryst['Aniso'][0]['start']>0:    #most crystals have no Anisotropic input
        print(">>> is Anisotropic", len(cryst["Aniso"]))

    # Consider anisotropic temperature factor
    # X.J. Yu, slsyxj@nus.edu.sg
    # A dummy dictionary Aniso with start =0 if no aniso temperature factor input
    # start
    hh = 4
    kk = 0
    ll = 0
    dspacing = material_constants_library.Crystal_dSpacing(cryst, hh, kk, ll) #bragg_metrictensor(cryst['a'], cryst['b'], cryst['c'], cryst['alpha'], cryst['beta'], cryst['gamma'], HKL=[hh, kk, ll])
    dspacing *= 1e-8 # in cm


    if 'Aniso' in cryst.keys() and cryst['Aniso'][0]['start']>0:    #most crystals have no Anisotropic input
        print(">>> is Anisotropic")
        TFac = TemperFactor( 1.0/(2.0*dspacing*1e8),cryst['Aniso'],Miller={'h':hh,'k':kk,'l':ll}, \
            cell={'a':cryst['a'],'b':cryst['b'],'c':cryst['c']},
                             n = cryst["n_atom"]
                             )
        B_TFac = 1
    else:
        B_TFac = 0


    return TFac, cryst

def check_structure_factor(descriptor="Si", hh=1, kk=1, ll=1, energy=8000,
                           do_assert=True, models=[1,1,1],
                           material_constants_library=None):

    from xoppylib.crystals.tools import bragg_calc2, crystal_fh
    from xoppylib.crystals.tools import bragg_calc
    # from xoppy_xraylib_util import bragg_calc


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
                            fileout="xcrystal.bra",
                            material_constants_library=material_constants_library)
        os.system("cp xcrystal.bra xcrystal_1.bra")
        dic1b = crystal_fh(dic1a, energy)

    #
    # Xiaojiang
    #
    if descriptor == "YB66":
        ANISO_SEL = 2   #  0: old Temper 1:Anisotropic  2: isotropic
        do_not_prototype = 0
    else:
        ANISO_SEL = 0
        do_not_prototype = 0

    if models[2]:
        # from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_calc as bragg_calc2
        # from orangecontrib.xoppy.util.xoppy_xraylib_util import crystal_fh
        dic2a = bragg_calc2(descriptor=descriptor, hh=hh, kk=kk, ll=ll, temper=1.0,
                            emin=energy-100, emax=energy+100, estep=5.0, ANISO_SEL=ANISO_SEL,fileout="xcrystal.bra",
                            do_not_prototype=do_not_prototype,
                            verbose=False, material_constants_library=material_constants_library)
        os.system("cp xcrystal.bra xcrystal_2.bra")

        if False:
            # to test reader...
            dic2aLoaded = load_bragg_preprocessor_file("xcrystal_2.bra")
            dic2a = dic2aLoaded #############################


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
        if models[0] and models[1]:
            print ('0,1: STRUCT', dic0b['STRUCT'], dic1b['STRUCT'] )
            print ('0,1: FH    ', dic0b['FH'],     dic1b['FH']     )
            print ('0,1: FH_BAR', dic0b['FH_BAR'], dic1b['FH_BAR'] )
            print ('0,1: F_0   ', dic0b['F_0'],    dic1b['F_0']    )
            assert (numpy.abs(dic1b['STRUCT'] - dic0b['STRUCT']) < 1e-3)
            assert (numpy.abs(dic1b['FH']     - dic0b['FH'])     < 1e-3)
            assert (numpy.abs(dic1b['FH_BAR'] - dic0b['FH_BAR']) < 1e-3)
            assert (numpy.abs(dic1b['F_0']    - dic0b['F_0'])    < 1)
        if models[1] and models[2]:
            print ('1,2: STRUCT', dic1b['STRUCT'], dic2b['STRUCT'] )
            print ('1,2: FH    ', dic1b['FH'],     dic2b['FH']     )
            print ('1,2: FH_BAR', dic1b['FH_BAR'], dic2b['FH_BAR'] )
            print ('1,2: F_0   ', dic1b['F_0'],    dic2b['F_0']    )
            assert (numpy.abs(dic2b['STRUCT'] - dic1b['STRUCT']) < 1e-3)
            assert (numpy.abs(dic2b['FH']     - dic1b['FH'])     < 1e-3)
            assert (numpy.abs(dic2b['FH_BAR'] - dic1b['FH_BAR']) < 1e-3)
            assert (numpy.abs(dic2b['F_0']    - dic1b['F_0'])    < 1)
        if descriptor == "YB66":
            values_from = 1 # 0=Xiaojiang, 1=srio
            print ('STRUCT', dic2b['STRUCT'])
            print ('FH    ', dic2b['FH'])
            print ('FH_BAR', dic2b['FH_BAR'])
            print ('F_0   ', dic2b['F_0'])
            if values_from == 1:
                # assert (numpy.abs(dic2b['STRUCT'] -  (565.7225232608029 + 35.9668881704435j))  < 1e-2)
                # assert (numpy.abs(dic2b['FH'] -     (565.7225232608029 + 35.966888170443404j)) < 1e-2)
                # assert (numpy.abs(dic2b['FH_BAR'] - (565.7225232608029 + 35.96688817044359j))  < 1e-2)
                # assert (numpy.abs(dic2b['F_0'] -    (8846.406209552279 + 56.12593721027547j))  < 0.3)
                if ANISO_SEL == 0:
                    assert (numpy.abs(dic2b['STRUCT'] - (570.0726764188605+36.24657824291629j))   < 1e-2)
                    assert (numpy.abs(dic2b['FH'] -     (570.0726764188606+36.2465782429162j))    < 1e-2)
                    assert (numpy.abs(dic2b['FH_BAR'] - (570.0726764188604+36.2465782429164j))    < 1e-2)
                    assert (numpy.abs(dic2b['F_0'] -    (8848.638071350848+56.12049122626639j))   < 0.3)
                elif ANISO_SEL == 1:
                    assert (numpy.abs(dic2b['STRUCT'] - (565.7226407626008+35.963615210235865j))   < 1e-2)
                    assert (numpy.abs(dic2b['FH'] -     (565.7226407626005+35.96361521023578j))    < 1e-2)
                    assert (numpy.abs(dic2b['FH_BAR'] - (565.7226407626013+35.96361521023595j))    < 1e-2)
                    assert (numpy.abs(dic2b['F_0'] -    (8848.638071350848+56.12049122626639j))    < 0.3)
                elif ANISO_SEL == 2:
                    assert (numpy.abs(dic2b['STRUCT'] - (565.5391037232481+35.9521062287469j))    < 1e-2)
                    assert (numpy.abs(dic2b['FH'] -     (565.5391037232482+35.9521062287468j))    < 1e-2)
                    assert (numpy.abs(dic2b['FH_BAR'] - (565.5391037232481+35.952106228747j))     < 1e-2)
                    assert (numpy.abs(dic2b['F_0'] -    (8848.638071350848+56.12049122626639j))   < 0.3)
            else:
                use_Atomic_name = True
                if not use_Atomic_name:
                    assert (numpy.abs(dic2b['STRUCT'] -  (563.4529619470779+35.8256810337139j))  < 1e-2)
                    assert (numpy.abs(dic2b['FH'] -      (563.4529619470779+35.82568103371383j)) < 1e-2)
                    assert (numpy.abs(dic2b['FH_BAR'] -  (563.4529619470779+35.82568103371397j))  < 1e-2)
                    assert (numpy.abs(dic2b['F_0'] -     (8848.638071350899+56.12049122626621j))  < 0.3)
                else:
                    assert (numpy.abs(dic2b['STRUCT'] -  (565.6124450891418+35.96361291284668j))  < 1e-2)
                    assert (numpy.abs(dic2b['FH'] -      (565.612444779105+35.96362149427959j)) < 1e-2)
                    assert (numpy.abs(dic2b['FH_BAR'] -  (565.6124453991785+35.96360433141376j))  < 1e-2)
                    assert (numpy.abs(dic2b['F_0'] -     (8842.035225507192+56.120491194441975j))  < 0.3)

    return dic2b['STRUCT']

if __name__ == "__main__":
    import os
    import xraylib
    from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_calc as bragg_calc_old
    #
    # crystal
    #

    from dabax.dabax_xraylib import DabaxXraylib
    dx = DabaxXraylib()

    #
    # test temperature
    #
    if True:
        TFac, cryst = check_temperature_factor(dx)

        print("TFac: ", TFac, cryst["n_atom"], len(TFac[0]), len(TFac[1]))

        len(cryst["Aniso"])
        START = []
        for i,ele in enumerate(cryst["Aniso"]):
            START.append(ele['start'])
            print(i,ele, TFac[:,START[-1]])

        print( ">>>>>different iso 0 values = ", len ( numpy.unique(TFac[0,:] , return_index=True)[1] ))
        print( ">>>>>different iso 1 values = ", len ( numpy.unique(TFac[1,:] , return_index=True)[1] ))
        print( ">>>>>different iso 0 values = ", len ( numpy.unique(TFac[0,START] , return_index=True)[1] ), numpy.unique(TFac[0,START] , return_index=True)[1])
        print( ">>>>>different iso 1 values = ", len ( numpy.unique(TFac[1,START] , return_index=True)[1] ), numpy.unique(TFac[1,START] , return_index=True)[1])

        atom = cryst['atom']
        number_of_atoms = len(atom)
        list_Zatom = [atom[i]['Zatom'] for i in range(len(atom))]
        list_Zatom = numpy.array(list_Zatom)
        list_Zatom = list_Zatom / list_Zatom.max()
        from srxraylib.plot.gol import plot
        plot(numpy.arange(TFac.shape[1]), TFac[0,:],
             numpy.arange(TFac.shape[1]), TFac[1,:],
             # numpy.arange(TFac.shape[1]), list_Zatom,
             # numpy.arange(TFac.shape[1]), TFac[2,:] / TFac[2,:].max(),
             xtitle='atom index', ytitle='temperature factor',
             legend=['isotropic','anisosotropic',], #'Z/max(Z)','start/max(start)']
             )
        print( ">>>>>different iso values = ", len ( numpy.unique(TFac[0,:] , return_index=True)[1] ))


    #
    # test crystal
    #
    if True:
        # test Si
        print("Testing Si...")
        check_structure_factor(descriptor="Si", hh=1, kk=1, ll=1, energy=8000,
                               material_constants_library=xraylib)

        # test Muscovite
        print("Testing Muscovite...")
        check_structure_factor(descriptor="Muscovite", hh=1, kk=1, ll=1, energy=8000, do_assert=1, models=[1,1,1],
                               material_constants_library=xraylib)

        # Test YB66
        print("Testing YB66...")
        check_structure_factor(descriptor="YB66", hh=4, kk=0, ll=0, energy=8040.0, do_assert=1, models=[0,0,1],
                               material_constants_library=dx)

    #
    # f0
    #

    if True:
        from xoppylib.decorators.xraylib_decorated import XraylibDecorated
        from xoppylib.decorators.dabax_decorated import DabaxDecorated

        xrl = XraylibDecorated()
        dx = DabaxDecorated()

        Si_xrl =  xrl.f0_calc(0, "Si", 0, 6, 100)
        Si_dbx =   dx.f0_calc(0, "Si", 0, 6, 100)
        H2O_xrl = xrl.f0_calc(1, "H2O", 0, 6, 100)
        H2O_dbx =  dx.f0_calc(1, "H2O", 0, 6, 100)
        H2O_xrl = xrl.f0_calc(2, "Water, Liquid", 0, 6, 100)
        # H2O_dbx =  dx.f0_calc(2, "Water, Liquid", 0, 6, 100) # not yet implemented


        from srxraylib.plot.gol import plot
        plot(Si_xrl["data"][0,:],Si_xrl["data"][1,:],
             Si_dbx["data"][0,:],Si_dbx["data"][1,:],
             H2O_xrl["data"][0, :], H2O_xrl["data"][1, :],
             H2O_dbx["data"][0, :], H2O_dbx["data"][1, :],
             linestyle=[None,'',None,''],
             marker=[None,'+',None,'+'],
             color=['r','r','b','b'],
             legend=['Si xraylib','Si dabax','H2O xraylib','H2O dabax'])


    # crystal tests
    if True:
        from dabax.common_tools import bragg_metrictensor
        cryst = dx.Crystal_GetCrystal(filename='Crystals.dat', entry_name='YB66')
        print(cryst['a'], cryst['b'], cryst['c'], cryst['alpha'], cryst['beta'], cryst['gamma'])
        mt = bragg_metrictensor(cryst['a'], cryst['b'], cryst['c'], cryst['alpha'], cryst['beta'], cryst['gamma'],
                                RETURN_REAL_SPACE=0,RETURN_VOLUME=0, HKL=None)
        print(mt, mt[0,0])