#!/usr/bin/python

"""The module of the ConfSet class."""

import sys
import math
import numpy as np
from strset import StrSet
import pdb

class ConfSet(object):
    """The class determine a set of configurations for FCI method."""

    def __init__(self, **kw):

        N_el = kw.get('Nel', 1)
        N_orb = kw.get('Norb', 5)
        spin = kw.get('spin', 0)

        self.str_a = StrSet(Nel=0.5 * N_el + spin, Norb = N_orb)
        self.str_b = StrSet(Nel=0.5 * N_el - spin, Norb = N_orb)
        self.num_conf = self.str_a.num_conf * self.str_b.num_conf
        self.conf = [[0 for x in range(4)] for x in range(self.num_conf)]

        for j1 in xrange(self.str_a.num_conf):
            for j2 in xrange(self.str_b.num_conf):
                self.conf[self.str_a.num_conf*j2+j1][0]=self.str_a.num_conf*j2+j1
                self.conf[self.str_a.num_conf*j2+j1][1]=j1
                self.conf[self.str_a.num_conf*j2+j1][2]=j2
                self.conf[self.str_a.num_conf*j2+j1][3]=1

    # ---------------------------------------------------------------------

    def ic(self, a, b):
        return self.str_a.num_conf * b + a

    def print_conf(self):
        print('Number of configurations:' + " " + str(len(self.conf)))
        print('---------------------------')
        print("| ind" + " |" + "   configuration   " + "|")
        print('---------------------------')

        for i in self.conf:
            print("   " + str(i[0]) + "    " + str(zip(self.str_a.conf[i[1]].id, self.str_b.conf[i[2]].id)))

        print('---------------------------')

    @staticmethod
    def get_phase(a1,b1,a2,b2):

        if (np.nonzero(a1.id)>np.nonzero(b1.id)):
            phase1 = -1
        else:
            phase1 = 1

        if (np.nonzero(a2.id)>np.nonzero(b2.id)):
            phase2 = -1
        else:
            phase2 = 1


        return phase1*phase2

    @staticmethod
    def get_exci_orb(a1,b1,a2,b2):
        # orb[0] the initial state, orb [1] is the final state

        orb1 = [iii for iii, e in enumerate(a1.id - a2.id) if e != 0]
        si1 = [e for iii, e in enumerate(a1.id - a2.id) if e != 0]

        orb2 = [iii for iii, e in enumerate(b1.id - b2.id) if e != 0]
        si2 = [e for iii, e in enumerate(b1.id - b2.id) if e != 0]

        if (len(orb1)==2):
            orb = orb1
            si = si1
            par = 1
        elif (len(orb2)==2):
            orb = orb2
            si = si2
            par = 0
        else:
            sys.exit("The problem is in the function get_exci_orb")


        if (si[0] > si[1]):
            orb = orb[::-1]

        if (orb[0] < 1):
            phase1 = 1-2*par
        else:
            phase1 = math.pow(-1, np.sum(a2.id[:(orb[0]-par)])+np.sum(b2.id[:(orb[0]-1)]))

        if (orb[1] < 1):
            phase2 = 1-2*par
        else:
            phase2 = math.pow(-1, np.sum(a1.id[:(orb[1]-par)])+np.sum(b1.id[:(orb[1]-1)]))


        if ( orb[0] > orb[1]):
            phase3 = -1
        else:
            phase3 = 1

        phase3 = 1
   #     pdb.set_trace()
        phase = phase1 * phase2 * phase3
        return (orb[0], orb[1], phase)


#---------------------------------------------------------------------

    @staticmethod
    def phase(a1,b1,a2,b2):
        # orb[0] the initial state, orb [1] is the final state

        orb1 = [iii for iii, e in enumerate(a1.id - a2.id) if e != 0]
        si1 = [e for iii, e in enumerate(a1.id - a2.id) if e != 0]

        orb2 = [iii for iii, e in enumerate(b1.id - b2.id) if e != 0]
        si2 = [e for iii, e in enumerate(b1.id - b2.id) if e != 0]

        if (len(orb1)==2):
            orb = orb1
            si = si1
            par = 1
        elif (len(orb2)==2):
            orb = orb2
            si = si2
            par = 0
        else:
            sys.exit("The problem is in the function get_exci_orb")


        if (si[0] > si[1]):
            orb = orb[::-1]

        if (orb[0] < 1):
            phase1 = 1-2*par
        else:
            phase1 = math.pow(-1, np.sum(a2.id[:(orb[0]-par)])+np.sum(b2.id[:(orb[0]-1)]))

        if (orb[1] < 1):
            phase2 = 1-2*par
        else:
            phase2 = math.pow(-1, np.sum(a1.id[:(orb[1]-par)])+np.sum(b1.id[:(orb[1]-1)]))


        if ( orb[0] > orb[1]):
            phase3 = -1
        elif ( orb[0] < orb[1]):
            phase3 = 1
        else:
            phase3 = 0
            print '00000000000000000000000000000000000000000'

#        phase3 = 1
   #     pdb.set_trace()
        phase = phase1 * phase2 * phase3
        return phase


#---------------------------------------------------------------------
if __name__ == "__main__":

    cs = ConfSet(Nel=2, Norb=6, M=0)
    cs.print_conf()
    print cs.str_a.et[0]
    print cs.str_a.et[1]
    print cs.str_a.et[2]
