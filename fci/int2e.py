#!/usr/bin/python

"""The module of the CoulombInt class."""

import numpy as np
import crys
from itertools import combinations_with_replacement as cr
from g_fit3d import GFit
import const
import sys
import contextlib
import os


class CoulombInt(object):

    """Computes Coulomb integrals."""

    def __init__(self, N=2, **kw):

        self.N = N                              # number of basis functions
        self.comb = self.specialCombinations()  # four factor products

        self._load = kw.get('pload', '/data/users/mklymenko/work_py/mb_project/')
        self._save = kw.get('psave', '/data/users/mklymenko/work_py/mb_project/')

        self.cint = []
        if os.path.isfile('/data/users/mklymenko/work_py/mb_project/coint2.dat'):
            os.remove('/data/users/mklymenko/work_py/mb_project/coint2.dat')
        try:               # try to read integrals from the disk
            with open(self._load + "coint2.dat"):
                self.cint = np.loadtxt(self._load + "coint2.dat")
        except IOError:    # there are no integrals on the disk

            wfs_x = []       # array of wave functions represented as a set of fitting parameters
            wfs_y = []       # array of wave functions represented as a set of fitting parameters
            wfs_z = []       # array of wave functions represented as a set of fitting parameters

            for j in xrange(self.N):
                wfs_x.append(GFit(qn=j, p='/data/users/mklymenko/science/H2_100/programing/dis/v0/', **kw)) # executing the fitting procedure N-times
                wfs_y.append(GFit(qn=j, p='/data/users/mklymenko/science/H2_100/programing/dis/v1/', **kw)) # executing the fitting procedure N-times
                wfs_z.append(GFit(qn=j, p='/data/users/mklymenko/science/H2_100/programing/dis/v2/', **kw)) # executing the fitting procedure N-times

            print(wfs_x[1].gf)
            # computing integrals for all possible combinations stored in self.comb
            self.cint = [self.comp_int(wfs_x[j[0]],wfs_x[j[1]],wfs_x[j[2]],wfs_x[j[3]]) for j in self.comb]
            self.save()

    #-----------------------------------
    def getInt(self, vec):
        i = self.comb.index(tuple(sorted(vec)))
        return self.cint[i]

    #-----------------------------------
    def specialCombinations(self):

        """This function computes the number of possible configurations
        of N basis functions into poducts of four functions."""

        vec = list(xrange(self.N))
        b = []
 #       for a1 in cr(vec, 2):
 #           for a2 in cr(vec, 2):
 #               b.append(a1+a2)
        for a in cr(vec, 4):
            b.append(a)
        return b

    #-----------------------------------

    @staticmethod
    def overlap_int(gf1,gf2):
        """Overlap integral for two s-gaussian functions"""

        r1=pow((gf1[0]-gf2[0]),2)+pow((gf1[1]-gf2[1]),2)+pow((gf1[2]-gf2[2]),2)
        return pow((const.pi/(gf1[3]+gf2[3])),(3/2))*np.exp(-(r1*gf1[3]*gf2[3])/(gf1[3]+gf2[3]))

    #----------------------------------

    @staticmethod
    def comp_int(wf1,wf2,wf3,wf4):

        integral=0

        for j1 in xrange(wf1._num_fu):
            for j2 in xrange(wf2._num_fu):
                for j3 in xrange(wf3._num_fu):
                    for j4 in xrange(wf4._num_fu):

#                         over0=CoulombInt.overlap_int(wf1.gf[j1],wf2.gf[j2])
#                         over1=CoulombInt.overlap_int(wf1.gf[j1],wf3.gf[j3])
#                         over2=CoulombInt.overlap_int(wf1.gf[j1],wf4.gf[j4])
#                         over3=CoulombInt.overlap_int(wf2.gf[j2],wf4.gf[j4])
#                         over4=CoulombInt.overlap_int(wf3.gf[j3],wf4.gf[j4])
#                         over5=CoulombInt.overlap_int(wf2.gf[j2],wf3.gf[j3])

#                         print(over0)
#                         print(over1)
#                         print(over2)
#                         print(over3)
#                         print(over4)
#                         print(over5)

                        #if (min([over0,over1,over2,over3,over4,over5])>=0.5):
                        integral += crys.coulomb_repulsion((wf1.gf[j1][0],wf1.gf[j1][1],wf1.gf[j1][2]),
                                                            wf1.gf[j1][4],(0,0,0),1.0/wf1.gf[j1][3],\
                                                           (wf2.gf[j2][0],wf2.gf[j2][1],wf2.gf[j2][2]),
                                                            wf2.gf[j2][4],(0,0,0),1.0/wf2.gf[j2][3],\
                                                           (wf3.gf[j3][0],wf3.gf[j3][1],wf3.gf[j3][2]),
                                                            wf3.gf[j3][4],(0,0,0),1.0/wf3.gf[j3][3],\
                                                           (wf4.gf[j4][0],wf4.gf[j4][1],wf4.gf[j4][2]),
                                                            wf4.gf[j4][4],(0,0,0),1.0/wf4.gf[j4][3])
                        #else:
                        #    integral +=0
        return integral
    #-----------------------------------

    def sum2(self, j, i):
        ans = 0
        for k in xrange(self.N):
            ans+=self.getInt([j,i,k,k])

        return ans

    def save(self):     # fixed number of functions, varied amplitudes and widths
        if (self._save != '0'):
            p = self._save + "coint2.dat"
            np.savetxt(p, self.cint)
        else:
            sys.exit("Wrong path to save")

if __name__ == "__main__":

    c=CoulombInt()
    print(np.array([c.comb,c.cint]).T)
