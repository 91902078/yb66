yb66
====

On going working for upgrading XOPPY/CRYSTAL and SHADOW with any crystal, in particular YB66


XOPPY
-----

XOPPY/CRYSTAL and also XOPPY/FH read a preprocessor file xcrystal.bra

This file is prepared with 

```
    from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_calc
    from orangecontrib.xoppy.util.xoppy_xraylib_util import crystal_fh
```

The file (or dictionary) created by bragg_calc() is used:
 
- in python crystal_fh() to calculate structure factors in XOPPY/FH
- in fortran diff_patt to calculate diffraction profiles in XOPPY/CRYSTAL

bragg_calc() and crystal_fh() have already been updated to run correctly crystals like Muscovite, that have sites with the same type of atom but different occupancies. 

The local version (slightly updated) of these routines are in a local copy of xoppy_xraylib_util.py

The bragg_calc routine has been upgraded for YB66 with new crystal data and f0 data. We use dabax files for that (http://ftp.esrf.fr/pub/scisoft/DabaxFiles/).
For the moment is called bragg_calc2 and will replace bragg_calc.  It can be found in xoppy_dabax_util.py.

The crystal_fh does not need an update.  

We will replace bragg_calc by bragg_calc2. 

SHADOW
------

with the python tools in XOPPY we should create the preprocessor file for shadow. This new preprocessor file requires also important changes in the SHADOW fortran kernel.

A new independent python code must be prepared to create the new preprocessor file (ShadowOui does not have XOPPY as dependence) and integrate it in the "bragg" preprocessor.


The SHADOW kernel will accept this new preprocessor file. Eventually, the fortran kernel would accept:
- The existing preprocessor file (not valid for new and complex crystals as YB66)
- The new preprocessor created with bragg_calc2 and also used in XOPPY/CRYSTAL and FH [STILL TO DO]

Then, we will add yb66 in the list of crystals of the bragg preprocessor. For that, a "quick-and dirty" way is to add it in the crystal list, and process it with the new python routines.  It is not very clean but it will make the job until a complete restructuration of the SHADOW crystal code in shadow4. A similar "patch" is already in use for Graphite.
