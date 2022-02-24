import numpy
from crystal_util import crystal_fh2, bragg_calc2


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:

        import argparse
        parser = argparse.ArgumentParser(description='Calculation structure factor')
        #args = parser.parse_args()


        parser.add_argument('-n','--name',dest='descriptor', default=['YB66'],type=str, nargs=1, help='Crystal name')
        parser.add_argument('-m','--m', metavar='H K L', default=[4,0,0],type=int, nargs='+', help='Miller indic [H, K, L]')
        parser.add_argument('-e','--e', dest='EngRange', metavar='emin,emax,estep', default=[8040,8050,1],type=float, nargs='+', help='[emin,emax,estep]')

        args = parser.parse_args()
        print(">>>>", args)
        descriptor = args.descriptor[0].strip()
        HMILLER = args.m[0]
        KMILLER = args.m[1]
        LMILLER = args.m[2]
        ENERGY = args.EngRange[0]
        ENERGY_END = args.EngRange[1]
        estep = args.EngRange[2]
        NPOINTS = int((ENERGY_END-ENERGY)/estep + 1)

        print("Using crystal descriptor: ",descriptor)
        bragg_dictionary = bragg_calc2(descriptor=descriptor,hh=HMILLER,kk=KMILLER,ll=LMILLER,temper=1.0,
                                                emin=ENERGY,emax=ENERGY_END,estep=estep,fileout=None)   #50eV, replaced with estep
        energy = numpy.linspace(ENERGY,ENERGY_END,NPOINTS)
        print("\nCrystal = %s, Miller Index = (%d,%d,%d)\n" % (descriptor,HMILLER,KMILLER,LMILLER))
        for i,ienergy in enumerate(energy):
            dic2 = crystal_fh2(bragg_dictionary,ienergy)
            print("Energy=%g eV FH=(%g,%g)"%(ienergy,dic2["STRUCT"].real,dic2["STRUCT"].imag))

    else:

        emin = 8040
        emax = 8050
        estep = 1


        bragg_dictionary = bragg_calc2(descriptor="YB66",
                                       hh=4,
                                       kk=0,
                                       ll=0,
                                       temper=1.0,
                                       emin=emin,
                                       emax=emax,
                                       estep=estep, # 50eV, replaced with estep
                                       fileout=None)

        energy = numpy.linspace(emin, emax, 1 + int( (emax-emin) / (estep)))
        print("\nCrystal = %s, Miller Index = (%d,%d,%d)\n" % ("YB66",4,0,0))
        for i, ienergy in enumerate(energy):
            dic2 = crystal_fh2(bragg_dictionary, ienergy)

            print("Energy=%g eV F(0,0,0)=%s, FH=(%g,%g)" % (ienergy, repr(dic2["F_0"]), dic2["STRUCT"].real, dic2["STRUCT"].imag))

