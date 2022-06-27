yb66
====

On going working for upgrading XOPPY/CRYSTAL and SHADOW with any crystal, in particular YB66


XOPPY
-----

XOPPY/CRYSTAL and also XOPPY/FH read a preprocessor file xcrystal.bra

This file was prepared with (now obsolete)

```
    from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_calc
    from orangecontrib.xoppy.util.xoppy_xraylib_util import crystal_fh
```

The file (or dictionary) created by bragg_calc() is used:
 
- in python crystal_fh() to calculate structure factors in XOPPY/FH
- in fortran diff_patt to calculate diffraction profiles in XOPPY/CRYSTAL


The new tools for creating the preprocessor file necessary for YB66 work using dabax instead of xraylib. They are integrated in a new package xoppylib. 
They can be installed: 

```
pip install dabax
pip install xoppylib
```

The new bragg_calc is called bragg_calc2:

```
    from xoppylib.crystals.tools import bragg_calc2, crystal_fh
```

a typical run of bragg_calc2:


```
    from dabax.dabax_xraylib import DabaxXraylib # imports dabax decorated with interface like xraylib
    import xraylib
    
    ENERGY = 8040.0
    SCANFROM = 0  # in microradiants
    SCANTO = 100  # in microradiants
    SCANPOINTS = 200
    # if using DABAX for YB66:
    bragg_dictionary = bragg_calc2(descriptor=descr"YB66",
                                       hh=4, kk=0, ll=0,
                                       temper=1.0,
                                       emin=ENERGY - 100.0, emax=ENERGY + 100.0,
                                       estep=(SCANTO - SCANFROM) / SCANPOINTS, fileout="xcrystal.bra",
                                       material_constants_library=DabaxXraylib())
                                       
    # Example using xraylib (for Si)
    bragg_dictionary = bragg_calc2(descriptor=descr"Si",hh=1, kk=1, ll=1,temper=1.0,
                                       emin=ENERGY - 100.0, emax=ENERGY + 100.0,
                                       estep=(SCANTO - SCANFROM) / SCANPOINTS, fileout="xcrystal.bra",
                                       material_constants_library=xraylib)
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


TODO (may be): 
- Include the site-dependent temperature factor in the widgets 


SHADOW
------

the python tools in XOPPY allows to create the preprocessor file for shadow (the one created with bragg_calc2). 

This new preprocessor file required important changes in the SHADOW fortran kernel. It accepts now both preprocessors file vesrions (v1 for the old one, v2 for the new one). The new code added for the v2 preprocessor is in the branch https://github.com/oasys-kit/shadow3/tree/devel-gfortran-yb66

This branch is created on top of the https://github.com/oasys-kit/shadow3/tree/devel-gfortran branch that solved some old compilation problems ()
This is not yet functional. 
file has been adapted from the diff_pat code https://github.com/srio/crystal

The new code is added to the "develop-gfortran" branch that solves the issue: https://github.com/oasys-kit/shadow3/issues/35

ShadowOui
---------

At the moment, to install OASYS/ShadowOui use the branch: https://github.com/oasys-kit/shadowOui/tree/yb66 This includes Beryl and Muscovite in the "bragg" "default" or traditional preprocessor (v1), and also implements YB66. In addition, it creates preprocessor file v2 with a crystal from two new lists: DABAX and XRayServer. It uses xoppylib as dependence, as it reuses the bragg_calc2 in xoppy.

There are a OASYS workspaces with YB66 examples and calculations for the paper in this repository (https://github.com/91902078/yb66/tree/main/workspaces/), in files: 
- https://raw.githubusercontent.com/91902078/yb66/main/workspaces/yb66_paper.ows
- https://raw.githubusercontent.com/91902078/yb66/main/workspaces/JUMBO_YB66.ows
- https://raw.githubusercontent.com/91902078/yb66/main/workspaces/YB66_VS_GRATING.ows


FIGURES AND TABLES
------------------

Figures and Tables for the paper used these scripts and worspaces:

- Table 1 use  https://github.com/91902078/yb66/paper_figs_and_tables/table1_explore.py https://github.com/91902078/yb66/paper_figs_and_tables/table1.py and *txt
- Table 3 use  https://github.com/91902078/yb66/paper_figs_and_tables/table3.py 
- Fig. 1 use https://github.com/91902078/yb66/paper_figs_and_tables/figure1.py
- Fig. 4,5  use the Oasys workspace https://github.com/91902078/yb66/workspaces/yb66_paper.ows
- Fig. 6 use https://github.com/91902078/yb66/paper_figs_and_tables/figure6.py
- Fig. 7,8  use the Oasys workspace https://github.com/91902078/yb66/workspaces/yb66_paper.ows and https://raw.githubusercontent.com/91902078/yb66/main/workspaces/JUMBO_YB66.ows


