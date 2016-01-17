#!/usr/bin/python

"""The module of the ConfSet class."""

import sys
import numpy as np
from xstr import XStr
from next_permutation import next_permutation
import math


class StrSet(object):
    """The class determine a set of configurations for FCI method."""

    def __init__(self, **kw):
        self.Nel = kw.get('Nel', 1)
        self.Norb = kw.get('Norb', 5)

        self.conf = []
        ind = []
        a = [0] * int(self.Norb - self.Nel) + [1] * int(self.Nel)

        for j in next_permutation(a):
            self.conf.append(XStr(id=np.array(j)))
            ind.append(self.conf[-1].index)

        # sort list of configurations
        ind, self.conf = zip(*sorted(zip(ind, self.conf)))
        self.num_conf = len(self.conf)
        # make the table of excitations
        self.excitab, self.et = self.exci_table()

    def exci_table(self):
        exta = np.zeros((len(self.conf), len(self.conf)), dtype='int')
        et = []

        for j1 in self.conf:
            auxi = []
            for j2 in self.conf:
                if (j1.index == j2.index):
                    exta[j1.index, j2.index] = np.sum(np.abs(j1.id))*10
                    for jj in [jj1 for jj1, e in enumerate(j1.id) if e != 0]:
                        auxi.append([j2.index, jj, jj, 1])
                else:
                    aux = j1.id - j2.id
                    if ((np.sum(np.abs(aux)) == 2) & (-1 in aux) & (1 in aux)):
                        orb1, orb2 = StrSet.get_exci_orb(j2, j1)
                        if j1.index >= j2.index:
                            exta[j1.index, j2.index] = -1
                            auxi.append([j2.index, orb1, orb2, -1])
                        else:
                            exta[j1.index, j2.index] = 1
                            auxi.append([j2.index, orb1, orb2, 1])
            et.append(auxi)

        return (exta, et)

    @staticmethod
    def get_exci_orb(a1, a2):

        orb = [iii for iii, e in enumerate(a1.id - a2.id) if e != 0]
        si = [e for iii, e in enumerate(a1.id - a2.id) if e != 0]

        if (len(orb) != 2):
            sys.exit("The problem is in the function get_exci_orb")

        if (si[0] > si[1]):
            orb = orb[::-1]

        return (orb[0], orb[1])

    def print_conf(self):
        print('Number of configurations:' + " " + str(len(self.conf)))
        print('---------------------------')
        print("| ind" + " |" + "   configuration   " + "|")
        print('---------------------------')

        for i in self.conf:
            print("   " + str(i.index) + "    " + str(i.id))

        print('---------------------------')
        print('   Table of excitations')
        print('---------------------------')
        print(self.excitab)

    @staticmethod
    def phase(a1, a2):
        # orb[0] the initial state, orb [1] is the final state

        orb = [iii for iii, e in enumerate(a1.id - a2.id) if e != 0]
        si = [e for iii, e in enumerate(a1.id - a2.id) if e != 0]

        if (len(orb) != 2):
            sys.exit("The problem is in the function get_exci_orb")

        if (si[0] > si[1]):
            orb = orb[::-1]

        phase1 = math.pow(-1, np.sum(a2.id[0:(orb[0])]))
        phase2 = math.pow(-1, np.sum(a2.id[0:(orb[1])]))

        if (orb[0] > orb[1]):
            phase3 = -1
        elif (orb[0] < orb[1]):
            phase3 = 1
        else:
            phase3 = 1

        phase = phase1 * phase2 * phase3

        return phase

# --------------------------------------------------------------------

if __name__ == "__main__":

    cs = StrSet(Nel=1, Norb=6, M=0)
    cs.print_conf()
    print cs.et
