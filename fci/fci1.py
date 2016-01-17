#!/usr/bin/python

import numpy as np
from confset import ConfSet
from aux_fun import ind3mat, norb
import sys
import pdb

class FCI(object):

    def __init__(self, **kw):

        self.Nel = kw.get('Nel', 2)        # number of interacting electrons
        self.spin = kw.get('M', 0)         # spin of the system

        # initialize properties of parent class which is the class of
        # two-electron integrals

        try:
            with open("/data/users/mklymenko/work_py/mb_project_bf/int2e.dat"):
                self.int2e = np.loadtxt("/data/users/mklymenko/work_py/mb_project_bf/int2e.dat")
        except IOError:
                sys.exit("Can not open the file int2e.dat")

        self.Norb = norb(np.size(np.diag(self.int2e)))


        # verify whether spin is consistent with number of electron

        if (((0.5 * self.Nel + self.spin) % 1 != 0) |
                ((0.5 * self.Nel - self.spin) % 1 != 0)):
            sys.exit("Cann't get desired magnitude of the spin projection")

        # initialize all possible alpha and betha strings

        self.str_a = ConfSet(Nel=0.5 * self.Nel + self.spin, Norb=self.Norb)
        self.str_b = ConfSet(Nel=0.5 * self.Nel - self.spin, Norb=self.Norb)

        # try to open files with computed single-particle spectrum
        try:
            with open("/data/users/mklymenko/work_py/mb_project_bf/E.dat"):
                self.E = np.loadtxt("/data/users/mklymenko/work_py/mb_project_bf/E.dat")[0:self.Norb]
        except IOError:
            sys.exit("Can not open the file E.dat")

        # try to open files with computed single-particle spectrum
        try:
            with open("/data/users/mklymenko/work_py/mb_project_bf/exch.dat"):
                self.exch = np.loadtxt("/data/users/mklymenko/work_py/mb_project_bf/exch.dat")
        except IOError:
            sys.exit("Can not open the file exch.dat")

        # initialize configuration interaction matrix by zeros
        self._M = np.zeros((self.str_a.num_conf * self.str_b.num_conf,
                            self.str_a.num_conf * self.str_b.num_conf))

    # ---------------------------------------------------------------------

    def ic(self, a, b):
        return self.str_a.num_conf * b + a

    # ------------------------Kinetic energy-------------------------------

    def diagM(self, a1, b1):
        self._M[self.ic(a1.index, b1.index), self.ic(a1.index, b1.index)] = \
            a1.get_energy(self.E)+b1.get_energy(self.E)

    # -------------------------Correlations--------------------------------

    def sigma1(self, a1, b1):

        for i in [i1 for i1, e in enumerate(self.str_b.excitab[b1.index, :]) if e != 0]:
            orb1, orb2 = ConfSet.get_exci_orb(self.str_b.conf[i], b1)
            print orb1, i
            self._M[self.ic(a1.index, b1.index), self.ic(a1.index, i)] += \
                self.exch[orb1, orb2] * self.str_b.excitab[b1.index, i]
#            pdb.set_trace()

            for j in [j1 for j1, e in enumerate(self.str_b.excitab[self.str_b.conf[i].index, :]) if e != 0]:

                orb3, orb4 = ConfSet.get_exci_orb(self.str_b.conf[j], self.str_b.conf[i])

                self._M[self.ic(a1.index, b1.index), self.ic(a1.index, j)] += \
                    0.5*self.getInt([orb1, orb2, orb3, orb4]) *\
                    self.str_b.excitab[b1.index, i]*self.str_b.excitab[b1.index, j]
                   # np.copysign(1.0, (orb2-orb1))*np.copysign(1.0, (orb4-orb3))

                print a1.index, b1.index, a1.index, j
                print self.ic(a1.index, b1.index), self.ic(a1.index, j)
                print '---------'

    def sigma2(self, a1, b1):

        for i in [i1 for i1, e in enumerate(self.str_a.excitab[a1.index, :]) if e != 0]:
            orb1, orb2 = ConfSet.get_exci_orb(self.str_a.conf[i], a1)

            self._M[self.ic(a1.index, b1.index), self.ic(i, b1.index)] +=\
                self.exch[orb1, orb2] * self.str_b.excitab[a1.index, i]

            print a1.index, b1.index, i, b1.index
            print self.ic(a1.index, b1.index), self.ic(i, b1.index)
            print '---------'

            for j in [j1 for j1, e in enumerate(self.str_a.excitab[self.str_a.conf[i].index, :]) if e != 0]:

                orb3, orb4 = ConfSet.get_exci_orb(self.str_a.conf[j], self.str_a.conf[i])

                self._M[self.ic(a1.index, b1.index), self.ic(j, b1.index)] +=\
                    0.5*self.getInt([orb1, orb2, orb3, orb4]) *\
                    self.str_b.excitab[a1.index, i]*self.str_b.excitab[a1.index, j]
                    #np.copysign(1.0, (orb2-orb1))*np.copysign(1.0, (orb3-orb4))

                print a1.index, b1.index, j, b1.index
                print self.ic(a1.index, b1.index), self.ic(j, b1.index)
                print '---------'

    def sigma3(self, a1, b1):

        for i in [i1 for i1, e in enumerate(self.str_a.excitab[a1.index, :]) if e != 0]:
            for j in [j1 for j1, e in enumerate(self.str_b.excitab[b1.index, :]) if e != 0]:

                orb1, orb2 = ConfSet.get_exci_orb(self.str_a.conf[i], a1)
                orb3, orb4 = ConfSet.get_exci_orb(self.str_b.conf[j], b1)

                self._M[self.ic(a1.index, b1.index), self.ic(i, j)] +=\
                    self.getInt([orb1, orb2, orb3, orb4]) *\
                    self.str_a.excitab[a1.index, i] *\
                    self.str_a.excitab[b1.index, j]

    def diagonalizer(self):
        for j1 in self.str_a.conf:
            for j2 in self.str_b.conf:
#               self.diagM(j1, j2)
                self.sigma1(j1, j2)  # nonherm
#                self.sigma2(j1, j2)  # nonherm
#                self.sigma3(j1, j2)

        E, wf = np.linalg.eig(self._M)

        print '-----------------------------'
        print E
        print '-----------------------------'
        print '-----------------------------'
        print '-----------------------------'
#
#        np.savetxt(paths['psave']+"M.dat", c._M)
#        np.savetxt(paths['psave']+"EE.dat", E)

        self.str_a.print_conf()

    def clean_M(self):
        self._M = np.zeros(
            (self.str_a.num_conf *
             self.str_b.num_conf,
             self.str_a.num_conf *
             self.str_b.num_conf))

    def getInt(self, vec):
        ind1 = ind3mat(vec[0],vec[1],self.Norb,'sym')
        ind2 = ind3mat(vec[2],vec[3],self.Norb,'sym')
        return self.int2e[ind1,ind2]


if __name__ == "__main__":

    c = FCI(Nel=2, M=0)
    orb1, orb2 = ConfSet.get_exci_orb(c.str_a.conf[4], c.str_a.conf[0])
    print c.str_a.conf[0]._id
    print c.str_a.conf[4]._id
    print orb1
    print orb2

#    print c.Norb
#    print c.E
#    print c.str_a
#    #c.diagM(c.str_a.conf[0],c.str_b.conf[1])
    c.diagonalizer()
    print c._M[9:16,9:16]
#    print c.sigma1(c.str_a.conf[0],c.str_b.conf[0])
#    print c.sigma2(c.str_a.conf[0],c.str_b.conf[0])

