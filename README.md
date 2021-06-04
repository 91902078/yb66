yb66
====

On going working for upgrading XOPPY/CRYSTAL and SHADOW with any crystal, in particular YB66


XOPPY
-----

XOPPY/CRYSTAL and also XOPPY/FH read a preprocessor file xbrystal.bra

This file is prepared with 

```
    from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_calc
    from orangecontrib.xoppy.util.xoppy_xraylib_util import crystal_fh
```

The file (or dictionary) created by bragg_calc() is used:
 
- in python crystal_fh() to calculate structure factors in XOPPY/FH
- in fortran diff_patt to calculate diffraction profiles in XOPPY/CRYSTAL

bragg_calc() and crystal_fh() have already been updated to run correctly crystals like Muscovite, that have sites with the same type of atom but different occupancies. 

The local version of these routines are in crystal_util_debug_xoppy.py

On-going work is being doing to upgrade these routines for yb66, but need to access new crystal constants from dabax files. 

They are in crystal_util.py than contains bragg_calc2() and crystal_fh2()

The final goal is to merge the "*2()" routines into the original ones.

For that, calc_bragg2() should make the job using dabax and be compatible with bragg_calc.
Eventually, calc_bragg should be an interface to calc_bragg_xraylib and calc_bragg_dabax

crystal_fh2 should merge in crystal_fh


SHADOW
------

with the python tools in XOPPY we should create the preprocessor file for shadow.
However, a new independent code must be prepared, as ShadowOui does not have XOPPY as dependence.


The SHADOW kernel will accept this new preprocessor file. Eventually, the fortran kernel would accept:
- The current preprocessor file
- The XOPPY/CRYSTAL preprocessor file xcrystal.bra
- The new preprocessor

Then, the bragg preprocessor should accept yb66. For that, a "quick-and dirty" way is to add it in the crystal list, ans process it with the new routines. A similar "patch" is already in use for Graphite. It is not clean but it will bame the job until a complete restructuration of the SHADOW crystal code in shadow4.