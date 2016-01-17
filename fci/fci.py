#!/usr/bin/python

import os
import numpy as np
from strset import StrSet
from confset import ConfSet
from aux_fun import ind3mat, norb, of
import sys
# import pdb


class FCI(object):

    def __init__(self, **kw):

        self.Nel = kw.get('Nel', 2)        # number of interacting electrons
        self.spin = kw.get('M', 0)         # spin of the system
        self.indi = kw.get('indi', 8)         # spin of the system

        # initialize properties of parent class which is the class of
        # two-electron integrals

        self.pth = os.path.dirname(os.getcwd())
        pth = self.pth + "/dis_scr/"
        self.int2e = of(pth + "!int2e"+str(self.indi)+".dat")
        self.Norb = norb(np.size(np.diag(self.int2e)))

        # verify whether spin is consistent with number of electron

        if (((0.5 * self.Nel + self.spin) % 1 != 0) |
                ((0.5 * self.Nel - self.spin) % 1 != 0)):
            sys.exit("Cann't get desired magnitude of the spin projection")

        self.E = of(pth + "!E"+str(self.indi)+".dat")[0:self.Norb]
        self.exch = -of(pth + "!exch"+str(self.indi)+".dat")[0:self.Norb,0:self.Norb]

        self.exch = self.exch + np.diag(self.E)

        # print np.diag(self.E)
        # self.Norb = 4

        # initialize all possible alpha and betha strings
        self.cs = ConfSet(Nel = self.Nel, Norb = self.Norb, spin = self.spin)

        # initialize configuration interaction matrix by zeros
        self._M = np.zeros((self.cs.num_conf, self.cs.num_conf))
        self._sym = np.zeros((self.cs.num_conf, self.cs.num_conf), dtype='S250')

    # ------------------------Kinetic energy------------------------------

    def diagM(self, a1, b1):

        self._M[self.cs.ic(a1.index, b1.index),
                self.cs.ic(a1.index, b1.index)] += \
            a1.get_energy(self.E) + b1.get_energy(self.E)

    # -------------------------Correlations--------------------------------

    def sigma1(self, a1, b1):

        for i in self.cs.str_b.et[b1.index]:

            orb1 = i[1]
            orb2 = i[2]

            if (b1.index == i[0]):
                phase1 = i[3]
            else:
                phase1 = StrSet.phase(self.cs.str_b.conf[i[0]], b1)


            self._M[self.cs.ic(a1.index, b1.index), self.cs.ic(a1.index, i[0])] += \
                self.exch[orb1, orb2] * phase1


#            print orb1, orb2, self.cs.ic(a1.index, b1.index), self.cs.ic(a1.index, i[0])

            for j in self.cs.str_b.et[i[0]]:
                orb3 = j[1]
                orb4 = j[2]
                if (j[0] == i[0]):
                    phase2 = j[3]
                else:
                    phase2 = StrSet.phase(self.cs.str_b.conf[j[0]],
                                          self.cs.str_b.conf[i[0]])

                self._M[self.cs.ic(a1.index, b1.index), self.cs.ic(a1.index, j[0])] += \
                    0.5*self.getInt([orb1, orb2, orb3, orb4]) * phase1 * phase2
#                print 0.5*self.getInt([orb1, orb2, orb3, orb4]) * phase1 * phase2


                # np.copysign(1.0, (orb2-orb1))*np.copysign(1.0, (orb4-orb3))
#                print (orb1, orb2, orb3, orb4)

    def sigma2(self, a1, b1):

        for i in self.cs.str_a.et[a1.index]:

            orb1 = i[1]
            orb2 = i[2]

            if (a1.index == i[0]):
                phase1 = i[3]
            else:
                phase1 = StrSet.phase(self.cs.str_a.conf[i[0]], a1)


            self._M[self.cs.ic(a1.index, b1.index), self.cs.ic(i[0], b1.index)] += \
                self.exch[orb1, orb2] * phase1

            for j in self.cs.str_a.et[i[0]]:
                orb3 = j[1]
                orb4 = j[2]
                if (j[0] == i[0]):
                    phase2 = j[3]
                else:
                    phase2 = StrSet.phase(self.cs.str_a.conf[j[0]],
                                          self.cs.str_a.conf[i[0]])


                self._M[self.cs.ic(a1.index, b1.index), self.cs.ic(j[0], b1.index)] += \
                    0.5*self.getInt([orb1, orb2, orb3, orb4]) * phase1 * phase2

                # np.copysign(1.0, (orb2-orb1))*np.copysign(1.0, (orb4-orb3))
                # print (orb1, orb2, orb3, orb4)

    def sigma3(self, a1, b1):

        for i in self.cs.str_a.et[a1.index]:
            for j in self.cs.str_b.et[b1.index]:

                orb1 = i[1]
                orb2 = i[2]
                orb3 = j[1]
                orb4 = j[2]

                if (a1.index == i[0]):
                    phase1 = i[3]
                else:
                    phase1 = StrSet.phase(self.cs.str_a.conf[i[0]], a1)


                if (b1.index == j[0]):
                    phase2 = j[3]
                else:
                    phase2 = StrSet.phase(self.cs.str_b.conf[j[0]], b1)


                self._M[self.cs.ic(a1.index, b1.index),
                        self.cs.ic(i[0], j[0])] +=\
                    self.getInt([orb1, orb2, orb3, orb4]) * phase1 * phase2

#                print (orb1, orb2, orb3, orb4)
#                print self.getInt([orb1, orb2, orb3, orb4])
#                print self.cs.ic(a1.index, b1.index), self.cs.ic(i[0], j[0])

    def diagonalizer(self):
        for j1 in self.cs.str_a.conf:
            for j2 in self.cs.str_b.conf:
                # continue
                # self.diagM(j1, j2)
                self.sigma1(j1, j2)  # nonherm
                self.sigma2(j1, j2)  # nonherm
                self.sigma3(j1, j2)

        E, wf = np.linalg.eig(self._M)
        #E=E.real
        idx = E.argsort()
        E = E[idx]
        wf = wf[:, idx]
        print '-----------------------------'
        print E
        print '-----------------------------'
        print np.sum(np.multiply(wf[:,0],wf[:,0]))
        print '-----------------------------'
        np.savetxt(self.pth + "/dis_scr/"+"Cfci"+str(self.indi)+".dat", wf)
        np.savetxt(self.pth + "/dis_scr/"+"E_fci"+str(self.indi)+".dat", E)
        #self.cs.print_conf()

    def clean_M(self):
        self._M = np.zeros((self.cs.num_conf, self.cs.num_conf))

    def getInt(self, vec):
        ind1 = ind3mat(vec[0], vec[1], self.Norb, 'sym')
        ind2 = ind3mat(vec[2], vec[3], self.Norb, 'sym')
        return self.int2e[ind1, ind2]

if __name__ == "__main__":

    for j in xrange(6,20):
        c = FCI(Nel=2, M=0, indi=j)
        c.diagonalizer()

#    print c.Norb
#    print c.E
#    j1=c.cs.conf
# c.sigma1(c.cs.str_a.conf[c.cs.conf[0][1]], c.cs.str_b.conf[c.cs.conf[0][2]])
# c.sigma2(c.cs.str_a.conf[c.cs.conf[0][1]], c.cs.str_b.conf[c.cs.conf[0][2]])
# c.sigma3(c.cs.str_a.conf[c.cs.conf[0][1]], c.cs.str_b.conf[c.cs.conf[0][2]])

    print c._M
    print c.exch
    print c.int2e
    print c.E

#    print c.Norb
#    print c.getInt([0, 0, 0, 0])
#    print c._sym
#    print c.cs.str_a.et
#    print c.cs.str_b.et
#    print c._M[9:16,9:16]
#    print np.char.decode(c._sym, 'UTF8')
#    print c.cs.str_b.conf[c.cs.conf[0][2]].id
#    print c.getInt([1,1,1,1])
#    print c.cs.str_b.conf[c.cs.conf[15][2]].id
#    print c.cs.str_b.conf[c.cs.conf[9][2]].id
#    print c.sigma1(c.str_a.conf[0],c.str_b.conf[0])
#    print c.sigma2(c.str_a.conf[0],c.str_b.conf[0])
