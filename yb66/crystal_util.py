import sys
import xraylib
import numpy
import os
import scipy.constants as codata
# X.J. Yu, slsyxj@nus.edu.sg
# from temperature_anisotropy import TemperFactor
from orangecontrib.xoppy.util.xoppy_xraylib_util import f0_xop
from dabax_util import crystal_parser
from dabax_util import get_f0_coeffs_from_dabax_file, calculate_f0_from_f0coeff
from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_metrictensor
import re
from scipy.optimize import curve_fit

#-------------------------------------------------------------------------
toangstroms = codata.h * codata.c / codata.e * 1e10



#
# __* are auxiliar routines not to be exported outside
#

#
# tools
#

def __symbol_to_from_atomic_number(ATOM):
    atoms = [ [1 ,"H"] ,[2 ,"He"] ,[3 ,"Li"] ,[4 ,"Be"] ,[5 ,"B"] ,[6 ,"C"] ,[7 ,"N"] ,[8 ,"O"] ,[9 ,"F"] ,[10 ,"Ne"], \
                 [11 ,"Na"] ,[12 ,"Mg"] ,[13 ,"Al"] ,[14 ,"Si"] ,[15 ,"P"] ,[16 ,"S"] ,[17 ,"Cl"] ,[18 ,"Ar"] ,[19 ,"K"]
             ,[20 ,"Ca"], \
                 [21 ,"Sc"] ,[22 ,"Ti"] ,[23 ,"V"] ,[24 ,"Cr"] ,[25 ,"Mn"] ,[26 ,"Fe"] ,[27 ,"Co"] ,[28 ,"Ni"] ,[29 ,"Cu"]
             ,[30 ,"Zn"], \
                 [31 ,"Ga"] ,[32 ,"Ge"] ,[33 ,"As"] ,[34 ,"Se"] ,[35 ,"Br"] ,[36 ,"Kr"] ,[37 ,"Rb"] ,[38 ,"Sr"] ,[39 ,"Y"]
             ,[40 ,"Zr"], \
                 [41 ,"Nb"] ,[42 ,"Mo"] ,[43 ,"Tc"] ,[44 ,"Ru"] ,[45 ,"Rh"] ,[46 ,"Pd"] ,[47 ,"Ag"] ,[48 ,"Cd"] ,[49 ,"In"]
             ,[50 ,"Sn"], \
                 [51 ,"Sb"] ,[52 ,"Te"] ,[53 ,"I"] ,[54 ,"Xe"] ,[55 ,"Cs"] ,[56 ,"Ba"] ,[57 ,"La"] ,[58 ,"Ce"] ,[59 ,"Pr"]
             ,[60 ,"Nd"], \
                 [61 ,"Pm"] ,[62 ,"Sm"] ,[63 ,"Eu"] ,[64 ,"Gd"] ,[65 ,"Tb"] ,[66 ,"Dy"] ,[67 ,"Ho"] ,[68 ,"Er"] ,[69 ,"Tm"]
             ,[70 ,"Yb"], \
                 [71 ,"Lu"] ,[72 ,"Hf"] ,[73 ,"Ta"] ,[74 ,"W"] ,[75 ,"Re"] ,[76 ,"Os"] ,[77 ,"Ir"] ,[78 ,"Pt"] ,[79 ,"Au"]
             ,[80 ,"Hg"], \
                 [81 ,"Tl"] ,[82 ,"Pb"] ,[83 ,"Bi"] ,[84 ,"Po"] ,[85 ,"At"] ,[86 ,"Rn"] ,[87 ,"Fr"] ,[88 ,"Ra"] ,[89 ,"Ac"]
             ,[90 ,"Th"], \
                 [91 ,"Pa"] ,[92 ,"U"] ,[93 ,"Np"] ,[94 ,"Pu"] ,[95 ,"Am"] ,[96 ,"Cm"] ,[97 ,"Bk"] ,[98 ,"Cf"] ,[99 ,"Es"]
             ,[100 ,"Fm"], \
                 [101 ,"Md"] ,[102 ,"No"] ,[103 ,"Lr"] ,[104 ,"Rf"] ,[105 ,"Db"] ,[106 ,"Sg"] ,[107 ,"Bh"] 	]

    if isinstance(ATOM ,int):
        for a in atoms:
            if a[0] == ATOM:
                return a[1]
    for a in atoms:
        if a[1] == ATOM:
            return int(a[0])

    raise Exception("Why are you here?")

def f0_with_fractional_charge(Z,AtomicName=None):
    """
    Returns the f0 value as stired in the f0_xop DABAX file:

    Parametrization of f0 (the non-dispersive part of the atomic scttering factor) vs sin(theta/lambda).
    This file contains a Waasmaier&Kirfel-like parametrization
    for the f0 data evaluated by Kissel using modified relativistic form
    factor [1]. For the fitting, model and the data error as in [2] are used.
    The parametrization was calculating by fitting the
    original data in the RTAB database [1] (also in the DABAX file
    f0_mf_Kissel.dat) using the MPFIT IDL procedure of Markwardt
    (http://cow.physics.wisc.edu/~craigm/idl/idl.html) usind the
    function [2]:
       f0[k] = c + [SUM a_i*EXP(-b_i*(k^2)) ]
                   i=1,5
     where k = sin(theta) / lambda and c, a_i and b_i
     are the coefficients tabulated in this file (in columns:
     a1  a2  a3  a4  a5  c  b1  b2  b3  b4  b5

    References
    [1] L. Kissel, Radiation physics and chemistry 59 (2000) 185-200
        available at:
        http://www-phys.llnl.gov/Research/scattering/RTAB.html
    [2] D. Waasmaier & A. Kirfel (Acta Cryst. (1995). A51, 416-413)


     Remarks:
     1) This file contains data for only neutral atoms. For ionic states
        other DABAX files (like f0_WaasKirf.dat) could be used.
     2) The coefficients in this file are prepared for xop2.1 and are
        unpublished.
     3) We created this file for being used as a default for XOP2.1. The
        reasons were two:
        i) It is more practical to use a parametrization rather than to use
           the direct data in file f0_mf.dat because of spped and storage
           reasons.
        ii) The previous defaulty file for XOP 2.0 was f0_WaasKirf.dat, which
           contained a published parametrization [2], but the original data
           is from the International Tables for X-ray Crystallography, vol C,
           (1995) Edited by AJC Wilson, Kluwer Academic Press, table 6.1.1.3
           which contains data from multiple origins (and references). We
           prefer to use evaluated data in the same conditions for all elements
           to keep some consistency in their origin.
           Additionally, Kissel data for f1f2 and CrossSec are also used in XOP.

    Column description: a1  a2  a3  a4  a5  c  b1  b2  b3 b4  b5

    :param Z: the atomic number
    :return: the 11 coefficients for buiding f0
    """
    tmp = {
    '1':[ 0.30426303,     0.45440815,     0.07977026,     0.15268426,     0.00871624,     0.00006937,     9.61768510,    23.52152246,     3.33237545,    53.49579285,     0.80301704],
    '2':[ 0.79971010,     0.63926998,     0.25291139,     0.27260951,     0.03475295,     0.00044196,    10.35354478,     3.93677741,     1.38554540,    25.50380187,     0.34363759],
    '3':[ 0.93188799,     0.17034671,     0.71602201,     0.35839485,     0.82020547,     0.00260460,     4.17494209,     0.33818181,    85.52950781,   175.66836887,     1.37984909],
    '4':[ 1.35769775,     0.71952631,     0.76849159,     0.14970777,     0.99985829,     0.00347440,    37.89541686,     0.66643295,    92.11060309,     0.17800503,     1.90321664],
    '5':[ 2.11021585,     0.94826030,     1.03175074,     0.17991800,     0.72039282,     0.00538888,    21.36228681,     1.17425000,    65.42872639,     0.12888999,     0.44259026],
    '6':[ 2.75491071,     1.06083859,     1.39000885,   -61.38969569,     0.72824752,    61.44921510,    13.55853062,     0.76518509,    43.40152845,    -0.00002249,     0.23756415],
    '7':[ 0.59731711,     2.28313328,     0.63580984,     2.26859990,     1.18920099,     0.02241857,     0.13665445,    17.86124963,    47.32019642,     7.67047709,     0.49262938],
    '8':[ 2.85010001,     2.54410664,     0.64919585,     0.79765894,     1.12558170,     0.02853916,    12.89518279,     5.45538842,     0.11296751,    35.32893278,     0.38721225],
    '9':[ 3.40542650,     2.84112029,     0.71156027,     0.94857450,     1.05115146,     0.03584797,     9.84065316,     4.07801707,     0.09701872,    27.57669889,     0.31753486],
    '10':[ 3.93418627,     3.17292749,     0.78988890,     1.08878497,     0.96122888,     0.04476065,     7.82568355,     3.17536738,     0.08620371,    22.24803408,     0.27139624],
    '11':[ 4.67888021,     3.31246434,     1.20302749,     1.11030781,     0.60080617,     0.07903453,     2.97160025,     8.01401365,     0.09917076,   119.87596938,     0.36875620],
    '12':[ 4.85834741,     1.50402245,     1.53555848,     0.95416075,     3.00404753,     0.13060425,     4.35307904,    95.26343316,     0.11044872,    31.99246407,     1.75670443],
    '13':[ 4.68669359,     2.38879059,     1.52287056,     1.07903978,     3.16354871,     0.14361992,     3.40015051,    36.83643520,     0.09557493,   114.27524168,     1.47492947],
    '14':[ 4.98816795,     3.35710271,     1.50292204,     1.22172882,     2.76143663,     0.15142442,     2.53600438,    29.97580504,     0.08254945,    88.73513838,     1.16712390],
    '15':[ 2.10509260,     4.36716206,     1.48028761,     1.33409382,     5.53780604,     0.15604602,     0.89460753,    24.41881968,     0.07145192,    69.85860169,     1.90671824],
    '16':[ 6.20239600,     5.39538819,     1.45887725,     1.42806579,     1.33894393,     0.15420384,     1.46090299,    20.14653917,     0.06164271,    56.73801104,     0.64755460],
    '17':[ 1.42704046,     6.77857375,     6.41364632,     1.52938606,     0.68790919,     0.13835865,     0.05185139,     1.14744676,    16.83234374,    47.08267330,     0.39585031],
    '18':[ 7.12496467,     6.96362577,     0.40535492,     1.69961663,     1.50025574,     0.27745601,     0.91908882,    14.15475105,    14.15474574,    39.07709506,     0.06417165],
    '19':[ 8.37197122,     7.07989232,     1.01347196,     0.77223727,     1.46667316,     0.26287737,    11.89949237,     0.78056461,   198.44072308,    41.56210501,     0.05445260],
    '20':[ 8.64392994,     1.46357896,     1.56066861,     1.05545005,     7.04009562,     0.19999674,     9.69070226,     0.04305143,    72.04776233,   166.55421240,     0.66621783],
    '21':[ 1.68383652,     1.47480433,     1.37047425,     9.24659290,     7.02207943,     0.16099853,    49.92976638,     0.03653814,   134.31889317,     8.45519757,     0.58503649],
    '22':[ 9.88746586,     1.50747326,     1.90635535,     1.54310395,     7.00376164,     0.10558296,     7.44889667,     0.03084452,    38.16657035,   115.68708856,     0.51848158],
    '23':[10.51693874,     1.57393530,     2.16699802,     1.69193955,     6.98396214,     0.01430187,     6.58802517,     0.02534021,    30.05889048,   101.44514004,     0.46207981],
    '24':[11.07627441,     1.58619669,     2.93288168,     1.37491086,     6.96664643,     0.00459082,     5.94437112,     0.02328655,    20.59115302,    96.88722975,     0.41908215],
    '25':[11.73712145,     1.82995875,     2.81141301,     1.90065220,     6.94066322,    -0.28355091,     5.23367741,     0.01636240,    20.53719746,    82.01631881,     0.37401385],
    '26':[12.30434401,     2.32437385,     3.18807287,     1.99400470,     6.92612145,    -0.80737786,     4.69095987,     0.01074597,    17.39159075,    74.61601310,     0.33856083],
    '27':[12.87035966,     4.43881362,     3.58294665,     2.06698965,     6.91586723,    -2.95253425,     4.22745007,     0.00461408,    15.08330800,    68.44820515,     0.30788472],
    '28':[13.40782101,     6.90042218,     3.99539727,     2.14745733,   283.90205482,  -282.43818444,     3.82465894,     0.28165040,    13.17017325,    62.81493664,     0.00006082],
    '29':[13.75399998,   297.69542224,     5.21908965,     1.58299921,     6.87528047,  -296.22007511,     3.49108135,     0.00005604,    11.26181897,    62.07734534,     0.26153431],
    '30':[14.43242183,     6.83363443,     4.86041106,     2.31023153,   518.74776544,  -517.28462315,     3.16729215,     0.24008925,    10.34604629,    53.61746648,     0.00002943],
    '31':[15.57306569,     6.80777698,     4.48158313,     2.55416620,   237.18610925,  -235.71998143,     2.94104939,     0.22275367,    11.54501025,    63.56972773,     0.00006123],
    '32':[16.29213242,   271.89501093,     4.17901049,     3.17804713,     6.76929988,  -270.43761631,     2.69207459,     0.00004947,    12.10034542,    55.59088173,     0.20581071],
    '33':[16.79861850,   336.43256975,     4.18039017,     3.72954296,     6.72535523,  -334.99650264,     2.44950901,     0.00003586,    12.83743210,    47.67503824,     0.18952317],
    '34':[17.16712180,   336.96736021,     4.71747652,     3.89487231,     6.68380463,  -335.56546318,     2.22463287,     0.00003080,    14.03902313,    42.86632703,     0.17412989],
    '35':[17.40231392,     5.93891344,   354.54400524,     3.52445295,     6.64817985,  -353.19765239,     2.01532010,    15.16569556,     0.00002303,    40.69766838,     0.15923958],
    '36':[17.53267157,     7.44816522,   267.89934293,     2.98742575,     6.61999042,  -266.63403399,     1.82191497,    15.41761348,     0.00002029,    39.30642110,     0.14476941],
    '37':[ 8.87288109,     1.56926290,     6.66073618,     1.08275119,    17.57340057,     1.08665457,    14.34390518,    35.97531308,     0.12903646,   209.51615220,     1.64279492],
    '38':[17.59670810,     9.94862630,     6.04462400,     2.61652635,     0.54472165,     1.07329335,     1.49577834,    13.12002382,     0.11972149,   119.44032856,     0.11972149],
    '39':[17.65536751,    10.43931521,     5.65885161,     3.13610892,     0.86203672,     1.06017183,     1.36826402,    12.12449612,     0.11133571,    95.69676914,     0.11133576],
    '40':[17.70896284,    11.09168541,     5.76424243,     3.49865793,     0.68990337,     1.04561396,     1.25477542,    11.37544137,     0.10372418,    80.89514696,     0.10372414],
    '41':[17.82410279,    12.34652807,     4.93750714,     3.19055176,     1.46016577,     1.02630932,     1.15717637,    11.35742646,     0.09662928,    69.39554321,     0.09662930],
    '42':[ 6.17787753,    17.84044504,    13.24473559,     3.35857771,     0.15154195,     0.99960433,     0.08963202,     1.06333880,    10.67451656,    61.58100361,     0.08963181],
    '43':[17.74278282,     3.79497777,     0.79808174,    13.24774635,     6.24221603,     0.95794137,     0.97113299,    42.86390771,   132.10191226,     9.26979654,     0.08236923],
    '44':[ 6.18155275,    17.80169136,    14.63083035,     3.59439161,     0.64428903,     0.91877951,     0.07624173,     0.89717697,     9.12737711,    36.30062507,   130.83719947],
    '45':[ 6.12780040,    17.79870560,     3.70551025,     0.61262679,    15.64710027,     0.86745289,     0.07028489,     0.82885382,    33.92215621,   126.42588650,     8.60698534],
    '46':[ 6.09085570,     4.68870690,    17.15641530,     3.84300712,    13.15123080,     0.80254284,     0.06457257,     0.76863619,     8.38997512,    34.11521928,     0.76863617],
    '47':[ 6.05724103,    17.76734177,     3.85910747,     0.54976811,    17.79132902,     0.70865971,     0.05857511,     7.66424959,    30.20483695,   121.19639558,     0.70979037],
    '48':[ 6.04036427,    18.38160827,     4.20316493,     0.75471621,    17.74206099,     0.59841397,     0.05297111,     7.01093020,    30.72014689,   106.80764634,     0.65631719],
    '49':[ 6.06154651,    19.00060821,     4.45915941,     1.04270193,    17.69494384,     0.44547768,     0.04732772,     6.43011640,    31.38897338,   129.84310238,     0.60802939],
    '50':[19.44406366,     6.14105459,     5.04240100,     1.19460095,    17.63952092,     0.22942056,     5.87081521,     0.04153040,    32.37864163,   114.94916342,     0.56396218],
    '51':[ 6.01179028,     6.25268840,    19.69743565,     1.16859113,    17.56401079,    -0.01748380,    32.17182048,     0.03634623,     5.33340196,   103.27699871,     0.52334570],
    '52':[ 7.27011998,     6.55126218,    19.86580814,     0.95229359,    17.52010619,    -0.49670130,    31.27680875,     0.02990172,     4.84670590,   100.88905997,     0.48406725],
    '53':[19.93807829,     7.09598483,     8.55254548,     0.81322939,    17.46003963,    -1.21189020,     4.40199492,     0.02378768,    29.37558678,    96.37033204,     0.44869824],
    '54':[19.95850428,     8.36420076,     9.90402346,     0.65245953,    17.40931309,    -2.65593783,     4.00466253,     0.01705338,    27.35376466,    96.62465241,     0.41685323],
    '55':[17.37569453,    12.13586009,    10.55273124,     1.29778936,    19.87680486,    -6.62800954,     0.38664179,     0.00966981,    23.44864968,   220.56960899,     3.64363491],
    '56':[19.82277247,    17.37211521,    10.64912854,     2.42340529,    51.54683902,   -46.21974320,     3.34809910,     0.36080812,    20.04789088,   161.05858343,     0.00185019],
    '57':[19.96684466,  2616.77020357,    11.08605976,     2.97061929,    17.30250838, -2611.52103692,     3.11105960,     0.00003320,    18.89301922,   126.40106401,     0.33961309],
    '58':[17.29595781,  4345.36340933,    20.94755728,     2.49690482,    11.57180940, -4340.11668462,     0.32551631,     0.00001928,     3.04729063,   146.75357391,    17.76655273],
    '59':[21.56386679,    17.23892858,    11.96705606,     2.55058606,  4645.78355064, -4640.56446888,     2.91250166,     0.30985101,    16.75073376,   139.14334527,     0.00001725],
    '60':[17.18849734,  5256.12632558,    12.33127623,     2.58243269,    22.21683582, -5250.92543157,     0.29580553,     0.00001461,    15.88877581,   133.41685838,     2.79253798],
    '61':[17.14790054,  5467.09133967,    12.67590249,     2.61445255,    22.89314617, -5461.92264966,     0.28225651,     0.00001335,    15.09888908,   128.17517281,     2.67966474],
    '62':[23.57776892,  4336.95334631,    13.01186004,     2.65112019,    17.09246464, -4331.80653530,     2.57205344,     0.00001606,    14.35914195,   123.21157749,     0.26967441],
    '63':[17.07135563,  5126.49069186,    13.33466886,     2.65426567,    24.27237487, -5121.36316596,     0.25889588,     0.00001297,    13.73081795,   120.24251547,     2.47764087],
    '64':[24.64625088,    16.99747848,    13.36466048,     3.34514581,  4543.57338365, -4538.49039355,     2.34262325,     0.24679775,    12.90194424,    95.39774603,     0.00001379],
    '65':[25.66612063,  4192.75632504,    13.95736847,     2.72815093,    17.00044486, -4187.69180909,     2.29102865,     0.00001424,    12.50615875,   111.05763849,     0.23787894],
    '66':[26.38323252,  4398.24494424,    14.26736529,     2.76532287,    16.97289395, -4393.23947589,     2.20308248,     0.00001260,    11.95434176,   106.89618428,     0.22734974],
    '67':[27.09405383,    16.93749917,    14.58634079,     2.80145206,  3698.14032223, -3693.18846082,     2.11867358,     0.21749694,    11.43043354,   102.92181598,     0.00001396],
    '68':[27.81053249,  4104.64353871,    14.88325935,     2.83483528,    16.89331430, -4099.71721767,     2.04181579,     0.00001196,    10.95422544,    99.42542607,     0.20889623],
    '69':[28.52850085,  2445.57145486,    15.22240776,     2.86116529,    16.92020688, -2440.77861007,     1.96383816,     0.00001716,    10.48926663,    96.46598501,     0.19884908],
    '70':[29.24441972,  2212.27369164,    15.49510957,     2.92480251,    16.89099012, -2207.52991878,     1.89336193,     0.00001759,    10.05627636,    91.69307079,     0.19084871],
    '71':[29.69788948,    15.47177518,  2097.86718749,     3.60694158,    16.82362572, -2093.19263162,     1.80307602,     9.47929223,     0.00001698,    77.64749351,     0.18222067],
    '72':[30.11062139,    15.40407022,  2274.03155954,     4.38097094,    16.78436412, -2269.46340922,     1.71565360,     8.96253008,     0.00001361,    63.95107044,     0.17352112],
    '73':[30.57052118,    15.57376402,  1353.88450796,     4.87958160,    16.74532612, -1349.42935153,     1.63430749,     8.68104102,     0.00001949,    56.83069217,     0.16513675],
    '74':[30.98784842,    15.84365824,  1085.37020439,     5.34491722,    16.75914099, -1081.10872457,     1.55377182,     8.43294944,     0.00001745,    50.66121806,     0.15616127],
    '75':[31.36764074,    16.20662904,  1140.89802522,     5.70714498,    16.65069139, -1136.66620413,     1.48081504,     8.26474305,     0.00001603,    45.63908603,     0.14966164],
    '76':[31.75193454,    16.78095382,   410.94969397,     5.94974979,    16.75723964,  -407.05071307,     1.40570152,     8.12076847,     0.00001747,    42.41784435,     0.13977093],
    '77':[31.76739790,     1.41742145,    16.74659192,    16.21969187,     6.34032022,     3.65295339,     1.32218506,    81.15672344,     0.13119994,     7.35439190,    26.57488749],
    '78':[30.98992203,    18.15501717,    16.67920729,     6.30866056,     1.30977395,     3.63949590,     1.27222945,     7.83558995,     0.12654931,    36.86734925,     1.27222941],
    '79':[16.88254232,    19.37581635,    32.62146048,     6.01941144,   432.23434409,  -429.07568590,     0.11690541,     7.87648778,     1.20474478,    32.97062267,    -0.00003099],
    '80':[16.43979396,    19.91863890,    27.76463255,     6.41022928,     4.90481781,     3.58984976,     0.11673200,     7.55331808,     1.15426441,    33.30365931,     1.15426437],
    '81':[16.27949672,    19.61124539,    32.52808287,     1.21893192,     6.82864702,     3.56321128,     0.11193237,     6.83979837,     1.09072759,   118.61695437,    25.09178834],
    '82':[16.15573594,    32.56471040,     6.74799686,     1.62318006,    20.37867502,     3.53108830,     0.10731472,     1.03828673,    25.84002437,   104.03280242,     6.51727715],
    '83':[16.05513439,    32.57219633,     7.10079270,     1.82264739,    20.92330655,     3.49908832,     0.10310727,     0.99041508,    26.79675111,    92.39516651,     6.16859916],
    '84':[15.92966393,    32.53463621,    21.26597071,     1.85348130,     7.90481807,     3.45521535,     0.09871917,     0.94233879,     5.77485734,    84.25758041,    27.14259799],
    '85':[15.82334418,    32.46748196,     9.01250570,     1.72434540,    21.46356994,     3.42272786,     0.09486320,     0.89843768,    27.00385753,    79.63992570,     5.39445822],
    '86':[15.71714768,    32.36554794,    21.55254018,     1.54577949,    10.31957180,     3.38282129,     0.09111052,     0.85619893,     5.01523482,    76.65511543,    26.25222808],
    '87':[15.61356578,    32.29005100,    21.58808798,     1.52079712,    11.49239752,     3.33918799,     0.08744016,     0.81690795,     4.68491604,   187.13075688,    24.42129865],
    '88':[32.26137613,    21.47411433,    11.54875240,     2.70070705,    15.53353668,     3.29342205,     0.78242755,     4.39338862,    21.26815521,   145.78453108,     0.08389021],
    '89':[15.44129352,    32.19461688,    21.67018557,    11.62375491,     3.60839956,     3.23867141,     0.08037991,     0.74781467,     4.14033299,    19.84028498,   113.68244239],
    '90':[15.35041360,    32.15264867,    21.98489988,     4.20304734,    11.87984118,     3.17150887,     0.07685215,     0.71518502,     3.92354148,    96.69516744,    19.25785705],
    '91':[32.41734383,    22.05085313,    13.10871855,     3.72910194,    15.28123528,     3.11865061,     0.68922416,     3.86872813,    17.74443730,   101.90118122,     0.07385957],
    '92':[15.20314929,    32.53208188,    13.81186567,     3.76642798,    22.31042147,     3.04559348,     0.07066435,     0.66171756,    16.79213936,    97.72588956,     3.75457984],
    '93':[32.64327679,    22.60430101,    14.43752068,     3.83052598,    15.12931302,     2.98794434,     0.63635924,     3.65636186,    15.88803466,    93.50959958,     0.06783444],
    '94':[32.89114822,    23.09839499,    15.57640146,     3.06849154,    15.06334425,     2.89811885,     0.61290312,     3.63716663,    15.44356030,   103.61518316,     0.06473308],
    '95':[33.02310917,    23.59414755,    16.08710719,     3.05160429,    15.00738866,     2.79583742,     0.58946182,     3.56398355,    14.80962925,   101.26355106,     0.06167575],
    '96':[14.92952838,    33.03254055,    24.00228529,     3.79349958,    16.07933216,     2.68191187,     0.05858286,     0.56470268,     3.42619535,    86.81878918,    13.93081918],
    '97':[14.87484815,    33.25429616,    24.75369621,     3.05890997,    16.98364667,     2.55564102,     0.05561473,     0.54346670,     3.41941648,    95.53275618,    13.62900417],
    '98':[33.35395184,    25.38419399,    17.37894160,     3.08843097,    14.82268135,     2.41071058,     0.52171734,     3.34704782,    13.07655187,    91.79465535,     0.05263332],
    'B-.':[1.493,         1.0472,         0.7776,          0.64929,       1.0233,          0.050981,       21.37,          65.436,        0.36215,        21.354,          1.1387],
    'Y3+':[17.92680,   9.153100,   1.767950,   -33.10800,   40.26020,  1.354170,   11.21450,   22.65990,  -1.3190000E-02 ]
    }
#    tmp = numpy.array(tmp)
#    return tmp[Z-1].copy()
    if Z > 0: return tmp[str(Z)].copy()
    return tmp[AtomicName].copy()  #should contain a atomic string name


def __func(q, a1, a2, a3, a4, a5, a6, a7, a8, a9):
    return calculate_f0_from_f0coeff([a1, a2, a3, a4, a5, a6, a7, a8, a9], q)

def __get_f0_coeffs(atoms, list_Zatom):
    """
    Return a Dict {"B-0.0455": [f0 coefficients], ..., "Y+3":[f0 coefficients],...}
    """
    AtomicChargeList = {}
    #first row is atomic number, it is integer
    UniqueAtomicNumber = list(sorted(set(list_Zatom)))
    charge = [ atoms[i]['charge'] for i in range(len(atoms))]
    for x in  UniqueAtomicNumber:
        AtomicChargeList[str(x)]= []
    for i,x in enumerate(list_Zatom):
        if charge[i] not in AtomicChargeList[str(int(x))]:
            AtomicChargeList[str(int(x))].append(charge[i])      #Charge value
    return __crystal_get_f0_coeffs(AtomicChargeList.items())

def __crystal_get_f0_coeffs(AtomicList):
    """
    Input: AtomicList, a list of tuple {(5,[-0.0455,]), (39,[3,])}, same atom allows with different charge
    Out:   A Dict {"B-0.0455": [f0 coefficients], ..., "Y+3":[f0 coefficients],...}
    """
    f0coeffs = {}
    searchChargeNameNeg = ['1-','2-','3-','4-','5-','6-','7-']
    searchChargeNamePos = ['1+','2+','3+','4+','5+','6+','7+']
    qq = numpy.linspace(0,2,1000)       #q = 0 to 2
    for x in AtomicList:
        n = int(x[0])       #atomic number
        sym = __symbol_to_from_atomic_number(n)
        f0 = get_f0_coeffs_from_dabax_file(entry_name=sym)
        if len(f0) == 0:
                raise("cannot find f0 coefficients for '" + sym + "'")
        for charge in x[1]: #may have multiple valences for same atom, B1+, B2+, etc
            k = int(charge)
            f01 = []
            if charge < 0:
                if k == charge:  #integer charge
                    f01 = get_f0_coeffs_from_dabax_file(entry_name=sym + searchChargeNameNeg[abs(k)-1])
                if len(f01) == 0:
                    ff = []
                    for i,s in enumerate(searchChargeNameNeg):
                        f01 = get_f0_coeffs_from_dabax_file(entry_name=sym + s)
                        if len(f01) > 0:
                            ff.append((-i-1,f01))
                            if (i+1) > abs(k):      #already find one with valence higher than atom carried charge
                                break
                    if len(ff) > 0:
                        f01 = ff[-1]

            if len(f01) == 0 and 0 != charge:  #not get a f0 in negative charge direction
                ff = []
                for i,s in enumerate(searchChargeNamePos):  #try to find one with positive charge
                    f01 = get_f0_coeffs_from_dabax_file(entry_name=sym + s)
                    if len(f01) > 0:
                        ff.append((i+1,f01))
                        if (i+1) > abs(k) or charge < 0:
                            break
                if len(ff) > 0:
                    f01 = ff[-1]

            if charge == 0: #always no fit for neutral atom
                f0coeffs[sym] = f0
                continue
            #following for charged atom
            if len(f01) == 0:
                raise("No 2nd atom found for linear fit f0 coefficients")
            if charge == f01[0]: #if charged atom already listed, just get it, no fit
                f0coeffs[sym+f'%+g'%charge] = f01[1]
                continue

            #do fitting here
            f0_1 = calculate_f0_from_f0coeff(f0, qq)
            f0_2 = calculate_f0_from_f0coeff(f01[1], qq)
            f00  = f0_1 + charge / f01[0] * (f0_2 - f0_1)
            p0 = f0         #neutral f0 for p0
            #if 2nd atom with valence closer to charge, use it instead of neutral atom
            if abs(charge-f01[0]) < abs(charge):
                p0 = f01[1]
            f00_fit, pcov_fit = curve_fit(__func, qq, f00, p0=p0)
            f0coeffs[sym+f'%+g'%charge] = f00_fit
    return  f0coeffs



## TODO: simplify this routine....
def __crystal_atnum(list_AtomicName, unique_AtomicName, unique_Zatom,list_fraction, f0coeffs):
    """
    To get the atom and fractional factor in diffierent sites
    list_AtomicName:  list of all atoms in the crystal
    unique_AtomicName:  list of unique atomicname in the list
    unique_Zatom:    list of unique atomic number
    list_fraction: list of unique fractial factor

    return: (num_e, fract, n_atom, n_name)
    (number of electrons, fraction, atomic sites, Unique name)
    """
    import re
    from orangecontrib.xoppy.util.xoppy_xraylib_util import f0_xop

    num_e = []
    fract = []
    n_atom = []
    n_ATUM = []
    n_name = []
    print(">>>>>> unique_AtomicName", unique_AtomicName)
    for k,x in enumerate(unique_AtomicName):
        tmp1 = re.search('(^[a-zA-Z]*)',x)
        if tmp1.group(0) == x:
            #original notation,AtomicName only, without valence info (i.e., B, Y, O)
            print(">>>> original notation,AtomicName only, without valence info (i.e., B, Y, O)")
            f0 = f0_xop(unique_Zatom[k])
        else:
            print(">>>> second notation")
            #f0 = f0_xop(0,AtomicName=x)
            if x in f0coeffs:  #second notation
                f0 = f0coeffs[x]
            else:   #first notation
                f0 = f0_with_fractional_charge(0,x)  #### TODO: fails, only one argument is passed

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
            n_name.append(x)

    return num_e.copy(), fract.copy(), n_atom.copy(),n_ATUM.copy(), n_name.copy()



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

# TODO: This part must be made compatible with old code (the of block should be removed)

    # this is not longer working... changed to a new flag
    # if len(atom[0]) >= 6:  #6 column + 1 AtomicName or +1 SeqNo (xraylib)
    total_charge = 0
    for i in range(len(atom)):
        total_charge += numpy.abs(atom[i]['charge'])

    if total_charge != 0:
        list_AtomicName = []
        for i in range(len(atom)):
            # tmp = atom[i]['AtomicName']
            # s = symbol_to_from_atomic_number(int(cell_data[0,i]))
            # if cell_data[5, i] != 0:  #charged
            #     s = s + f'%+.6g'%cell_data[5, i]
            s = __symbol_to_from_atomic_number(atom[i]['Zatom'])
            s = s + f'%+.6g' % atom[i]['charge']
            list_AtomicName.append( s )
        # TODO: check this, it fails. May be because it does not know list_AtomicName ?
        f0coeffs = __get_f0_coeffs(atom, list_Zatom)

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
            unique_Zatom.append(__symbol_to_from_atomic_number(tmp.group(0)))
##  ------------ Singapore Synchrotorn Light Source ---------------------------------
    TmpCrystal = () # for diff_pat.exe
    #TODO: this has to be modified to make it working for old and new code
    if unique_AtomicName[0] !='':   #Complex crystal
        TmpCrystal = __crystal_atnum(list_AtomicName, unique_AtomicName, unique_Zatom,list_fraction,f0coeffs)
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
    print(">>>>>", emax,emin,estep)
    npoint  = int(numpy.ceil(( (emax - emin)/estep + 1 )))      
    txt += "# The number of energy points NPOINT: \n"
    txt +=  ("%i \n") % npoint
    output_dictionary["npoint"] = npoint
    txt += "# for each energy point, energy, F1(1),F2(1),...,F1(nbatom),F2(nbatom)\n"
    list_energy = []
    print(">>>> npoint: ", npoint)
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


if __name__ == "__main__":

    if False:
        #
        # old code Si
        #
        from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_calc, crystal_fh
        dic1a = bragg_calc(descriptor="Si",hh=1,kk=1,ll=1,temper=1.0,emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra")
        print(dic1a)

        dic1b = crystal_fh(dic1a,8000.0)

        #
        # New code Si
        #

        dic2a = bragg_calc2(descriptor="Si",hh=1,kk=1,ll=1,temper=1.0,emin=7900.0,emax=8100.0,estep=5.0,fileout="xcrystal.bra")
        print("KEYS: ",dic2a.keys())
        print(dic2a)

        dic2b = crystal_fh(dic2a,8000.0)
        print(dic2b["info"])
        print("KEYS: ",dic2b.keys())


        for key in dic1b.keys():
            if key != "info":
                print(">>>", key,dic1b[key],dic2b[key])


    #
    # New code YB66
    #

    if True:
        dic3a = bragg_calc2(descriptor="YB66",hh=4,kk=0,ll=0,temper=1.0,emin=5000.0,emax=15000.0,estep=100,fileout="xcrystal.bra")
        print("KEYS: ",dic3a.keys())
        print(dic3a)

        dic3b = crystal_fh2(dic3a,8040.0)

        print(dic3b['info'])
