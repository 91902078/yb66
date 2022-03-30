import numpy
import os

from srxraylib.plot.gol import plot, set_qt
from xoppylib.decorators.xraylib_decorated import XraylibDecorated
from xoppylib.decorators.dabax_decorated import DabaxDecorated


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
        from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_calc as bragg_calc_old
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

    set_qt()

    xrl = XraylibDecorated()
    dx = DabaxDecorated()

    print(dx.info())

    do_plot = 1


#
# table 3
#


    # test Muscovite
    F200 = numpy.abs( check_structure_factor(descriptor="YB66", hh=2, kk=0, ll=0, energy=8040.0, do_assert=0, models=[0,0,1], material_constants_library=dx) )
    F400 = numpy.abs( check_structure_factor(descriptor="YB66", hh=4, kk=0, ll=0, energy=8040.0, do_assert=0, models=[0,0,1], material_constants_library=dx) )
    F600 = numpy.abs( check_structure_factor(descriptor="YB66", hh=6, kk=0, ll=0, energy=8040.0, do_assert=0, models=[0,0,1], material_constants_library=dx) )
    F800 = numpy.abs( check_structure_factor(descriptor="YB66", hh=8, kk=0, ll=0, energy=8040.0, do_assert=0, models=[0,0,1], material_constants_library=dx) )

    F200_calc_ref = 8.8 * 4
    F400_calc_ref = 137.9 * 4
    F600_calc_ref = 11.2 * 4
    F800_calc_ref = 4.8 * 4

    F200_obs_ref = 6.9 * 4
    F400_obs_ref = 150.4 * 4
    F600_obs_ref = 11.9 * 4
    F800_obs_ref = 7.2 * 4

    txt = ""
    txt += "F_200 this work: F=%g, calculated(ref): Fc=%g, observed(ref): Fo=%g, |F/Fc|=%g\n" % (F200, F200_calc_ref, F200_obs_ref, F200/F200_calc_ref)
    txt += "F_400 this work: F=%g, calculated(ref): Fc=%g, observed(ref): Fo=%g, |F/Fc|=%g\n" % (F400, F400_calc_ref, F400_obs_ref, F400/F400_calc_ref)
    txt += "F_600 this work: F=%g, calculated(ref): Fc=%g, observed(ref): Fo=%g, |F/Fc|=%g\n" % (F600, F600_calc_ref, F600_obs_ref, F600/F600_calc_ref)
    txt += "F_800 this work: F=%g, calculated(ref): Fc=%g, observed(ref): Fo=%g, |F/Fc|=%g\n" % (F800, F800_calc_ref, F800_obs_ref, F800/F800_calc_ref)

    print(txt)
    f = open("table3.txt", 'w')
    f.write(txt)
    f.close()
    print("File table3.txt written to disk.")
    #
    # no prototype:
    # F_200 this work: F=39.2492, calculated(ref): Fc=35.2, observed(ref): Fo=27.6, |F/Fc|=1.11503
    # F_400 this work: F=563.74, calculated(ref): Fc=551.6, observed(ref): Fo=601.6, |F/Fc|=1.02201
    # F_600 this work: F=38.893, calculated(ref): Fc=44.8, observed(ref): Fo=47.6, |F/Fc|=0.868148
    # F_800 this work: F=24.8411, calculated(ref): Fc=19.2, observed(ref): Fo=28.8, |F/Fc|=1.29381
    #
    # with prototype ** DEFAULT **
    # F_200 this work: F=39.2639, calculated(ref): Fc=35.2, observed(ref): Fo=27.6, |F/Fc|=1.11545
    # F_400 this work: F=564.591, calculated(ref): Fc=551.6, observed(ref): Fo=601.6, |F/Fc|=1.02355
    # F_600 this work: F=39.0252, calculated(ref): Fc=44.8, observed(ref): Fo=47.6, |F/Fc|=0.871098
    # F_800 this work: F=24.9914, calculated(ref): Fc=19.2, observed(ref): Fo=28.8, |F/Fc|=1.30163
    #
    # no temperature factor (identical results useing no prototypical or with prototypical)
    # F_200 this work: F=39.3788, calculated(ref): Fc=35.2, observed(ref): Fo=27.6, |F/Fc|=1.11871
    # F_400 this work: F=571.224, calculated(ref): Fc=551.6, observed(ref): Fo=601.6, |F/Fc|=1.03558
    # F_600 this work: F=40.0644, calculated(ref): Fc=44.8, observed(ref): Fo=47.6, |F/Fc|=0.894294
    # F_800 this work: F=26.1867, calculated(ref): Fc=19.2, observed(ref): Fo=28.8, |F/Fc|=1.36389
    #
    # isotropic temperature factor
    # F_200 this work: F=39.2584, calculated(ref): Fc=35.2, observed(ref): Fo=27.6, |F/Fc|=1.1153
    # F_400 this work: F=564.274, calculated(ref): Fc=551.6, observed(ref): Fo=601.6, |F/Fc|=1.02298
    # F_600 this work: F=38.976, calculated(ref): Fc=44.8, observed(ref): Fo=47.6, |F/Fc|=0.87
    # F_800 this work: F=24.9354, calculated(ref): Fc=19.2, observed(ref): Fo=28.8, |F/Fc|=1.29872
    #
    #
    #
    #
