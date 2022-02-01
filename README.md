yb66
====

On going working for upgrading XOPPY/CRYSTAL and SHADOW with any crystal, in particular YB66


XOPPY (current)
---------------

XOPPY/CRYSTAL and also XOPPY/FH read a preprocessor file xcrystal.bra

This file was prepared with 

```
    from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_calc
    from orangecontrib.xoppy.util.xoppy_xraylib_util import crystal_fh
```

The file (or dictionary) created by bragg_calc() is used:
 
- in python crystal_fh() to calculate structure factors in XOPPY/FH
- in fortran diff_patt to calculate diffraction profiles in XOPPY/CRYSTAL

bragg_calc() and crystal_fh() have already been updated to run correctly crystals like Muscovite, that have sites with the same type of atom but different occupancies. 


XOPPY (new)
-----------

XOPPY/CRYSTAL and also XOPPY/FH read a preprocessor file xcrystal.bra


The new tools for creating the preprocessor file necessary for YB66 work using dabax instead of xraylib. They are integrated in a new package xoppylib. 
They should be installed: 

```
pip install dabax
pip install xoppylib
```


For the moment, the new bragg_calc is called bragg_calc2:

```
    from xoppylib.crystals.tools import bragg_calc2, crystal_fh
```

a typical run of bragg_calc2:


```
    from dabax.dabax_xraylib import DabaxXraylib # imports dabax decorated with interface like xraylib
    
    ENERGY = 8040.0
    SCANFROM = 0  # in microradiants
    SCANTO = 100  # in microradiants
    SCANPOINTS = 200
    bragg_dictionary = bragg_calc2(descriptor=descr"YB66",
                                       hh=4, kk=0, ll=0,
                                       temper=1.0,
                                       emin=ENERGY - 100.0, emax=ENERGY + 100.0,
                                       estep=(SCANTO - SCANFROM) / SCANPOINTS, fileout="xcrystal.bra",
                                       material_constants_library=DabaxXraylib())
```


In Oasys, one can use this dictionary to calculate the structure factors (using crystal_fh) or the diffraction pattern (using run_diff_pat):


```

from xoppylib.crystals.tools import run_diff_pat, crystal_fh

    run_diff_pat(
        bragg_dictionary,
        preprocessor_file="xcrystal.bra",
        MOSAIC=0,
        GEOMETRY=0,
        SCAN=2,
        UNIT=1,
        SCANFROM=SCANFROM,
        SCANTO=SCANTO,
        SCANPOINTS=SCANPOINTS,
        ENERGY=ENERGY,
        ASYMMETRY_ANGLE=0.0,
        THICKNESS=0.7)

    from srxraylib.plot.gol import plot

    data = numpy.loadtxt("diff_pat.dat", skiprows=5)
    plot(data[:, 0], data[:, -1])
    
    dic0b = crystal_fh(bragg_dictionary, 8000.0)
    print(dict0['info'])
    
        
```

Examples are in the directory https://github.com/91902078/yb66/yb66/new

SHADOW
------

the python tools in XOPPY allows to create the preprocessor file for shadow (the one created with bragg_calc2). 

This new preprocessor file required important changes in the SHADOW fortran kernel. It acceptrs now both preprocessors file vesrions (v1 for the old one, v2 for the new one)

To build shadow3 for yb66 use the branch https://github.com/oasys-kit/shadow3/tree/devel-gfortran-yb66

(ShadowOui does have xoppylib as dependence! YB66 is added to the list of crystals in the bragg preprocessor.)


To build OASYS/ShadowOui use the branch: https://github.com/oasys-kit/shadowOui/tree/yb66

To build OASYS/XOPPY use the branch https://github.com/oasys-kit/xoppy/tree/xoppylib
(New crystal YB66 is added to the list of crystals. It uses xraylib wwhen possible, otherwise dabax is used, like for YB66)




