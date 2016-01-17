#!/usr/bin/python

import math
import numpy as np
from itertools import combinations_with_replacement as cr
from PyQuante import Ints
from g_fit3d import GFit
from invdisttree import Invdisttree
import matplotlib.pyplot as plt
from wf import WF
import copy

class SDet(object):
    """Single slater determinant class"""
    def __init__(self,**kw):

        self._alpha_id=kw.get('d_a_id',np.array([1,0,0,0]))    #id of the string
        self._betha_id=kw.get('d_b_id',np.array([1,0,0,0]))    #id of the string

        self.K=len(self._d_id)                           #number of orbitals

        self.Na=np.sum(self._alpha_id)                   #number of electrons in the string
        self.Nb=np.sum(self._betha_id)

        self._index_a=self.get_index(self._alpha_id)     #index the string
        self._index_b=self.get_index(self._betha_id)     #index the string

    @property
    def d_id(self):
        """I'm the 'x' property."""
        return self._d_id

    @d_id.setter
    def d_id(self, value):
        self._d_id = value


    def get_index(self,id):
        """this function assigns the unique index to the configuration """
        N=np.sum(id)
        auxiliary_mat=np.ones((self.K-N+1,N+1))

        for orb_ind in xrange(1,(self.K-N+1)):
            for el_ind in xrange(1,(N+1)):
                auxiliary_mat[orb_ind,el_ind]=auxiliary_mat[orb_ind-1,el_ind]+auxiliary_mat[orb_ind,el_ind-1]

        YY=np.zeros((self.K-N+1,N))
        YY[1:,:]=auxiliary_mat[:-1,1:]

        j1=0;j2=0;ind=0

        for orb_ind in xrange(self.K):
            if (self._d_id[orb_ind]==1):
                ind=ind+YY[j2,j1]
                j1+=1
            else:
                j2+=1

        return int(ind)



    def excitation(self,j,i):

        det1=copy.deepcopy(self)

        if (det1.d_id[j]==0):
            det1.d_id[j]=det1.d_id[j]+1
        else:
            det1.d_id[j]=0

        if (det1.d_id[i]==0):
            det1.d_id[i]=det1.d_id[i]+1
        else:
            det1.d_id[i]=0

        return det1

class ConfSet:
    """The class determine a set of configurations for FCI method"""

    def __init__(self,**kw):

        self.Nel=kw.get('Nel',2)
        self.Norb=kw.get('Norb',3)
        self.spin=kw.get('M',0)

        self.Nel_a=0.5*self.Nel+self.spin    #number of electrons in alpha's string
        self.Nel_b=0.5*self.Nel-self.spin    #number of electrons in betha's string

        def next_permutation(seq, pred=cmp):
            """Like C++ std::next_permutation() but implemented as
            generator. Yields copies of seq."""
            def reverse(seq, start, end):
                # seq = seq[:start] + reversed(seq[start:end]) + \
                #       seq[end:]
                end -= 1
                if end <= start:
                    return
                while True:
                    seq[start], seq[end] = seq[end], seq[start]
                    if start == end or start+1 == end:
                        return
                    start += 1
                    end -= 1
            if not seq:
                raise StopIteration
            try:
                seq[0]
            except TypeError:
                raise TypeError("seq must allow random access.")
            first = 0
            last = len(seq)
            seq = seq[:]
            # Yield input sequence as the STL version is often
            # used inside do {} while.
            yield seq
            if last == 1:
                raise StopIteration
            while True:
                next = last - 1
                while True:
                    # Step 1.
                    next1 = next
                    next -= 1
                    if pred(seq[next], seq[next1]) < 0:
                        # Step 2.
                        mid = last - 1
                        while not (pred(seq[next], seq[mid]) < 0):
                            mid -= 1
                        seq[next], seq[mid] = seq[mid], seq[next]
                        # Step 3.
                        reverse(seq, next1, last)
                        # Change to yield references to get rid of
                        # (at worst) |seq|! copy operations.
                        yield seq[:]
                        break
                    if next == first:
                        raise StopIteration
            raise StopIteration

        self.conf_a=[]
        a=[0]*int(self.Norb-self.Nel_a)+[1]*int(self.Nel_a)

        for j in next_permutation(a):
            self.conf_a.append(SDet(d_id=np.array(j)))

        self.conf_b=[]
        a=[0]*int(self.Norb-self.Nel_b)+[1]*int(self.Nel_b)

        for j in next_permutation(a):
            self.conf_b.append(SDet(d_id=np.array(j)))

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

def Coulomb_integrals(N):
    """computes Coulomb integrals"""
    #-----------------------------------
    def specialCombinations(vec):

        b=[];j=0
        for a in cr(vec,4):
            b.append(a)
            j+=1
        return b
    #-----------------------------------
    def overlap_int(gf1,gf2):
        """Overlap integral for two s-gaussian functions"""

        r1=pow((gf1[0]-gf2[0]),2)+pow((gf1[1]-gf2[1]),2)+pow((gf1[2]-gf2[2]),2)
        return pow((pi/(gf1[3]+gf2[3])),(3/2))*math.exp(-(r1*gf1[3]*gf2[3])/(gf1[3]+gf2[3]))
    #----------------------------------
    def comp_int(wf1,wf2,wf3,wf4):

        integral=0

        for j1 in xrange(wf1.N):
            for j2 in xrange(wf2.N):
               for j3 in xrange(wf3.N):
                   for j4 in xrange(wf4.N):

                        over0=overlap_int(wf1.gf[j1],wf2.gf[j2])
                        over1=overlap_int(wf1.gf[j1],wf3.gf[j3])
                        over2=overlap_int(wf1.gf[j1],wf4.gf[j4])
                        over3=overlap_int(wf2.gf[j2],wf4.gf[j4])
                        over4=overlap_int(wf3.gf[j3],wf4.gf[j4])
                        over5=overlap_int(wf2.gf[j2],wf3.gf[j3])

                        if (min([over0,over1,over2,over3,over4,over5])>=1.0):
                            integral += Ints.coulomb_repulsion((wf1.gf[j1][0],wf1.gf[j1][1],wf1.gf[j1][2]),wf1.gf[j1][4],(0,0,0),wf1.gf[j1][3],\
                                                               (wf2.gf[j2][0],wf2.gf[j2][1],wf2.gf[j2][2]),wf2.gf[j2][4],(0,0,0),wf2.gf[j2][3],\
                                                               (wf3.gf[j3][0],wf3.gf[j3][1],wf3.gf[j3][2]),wf3.gf[j3][4],(0,0,0),wf3.gf[j3][3],\
                                                               (wf4.gf[j4][0],wf4.gf[j4][1],wf4.gf[j4][2]),wf4.gf[j4][4],(0,0,0),wf4.gf[j4][3])
                        else:
                            integral +=0

        return integral
    #-----------------------------------

    N=5 #number of basis functions

    a=list(xrange(N))#this is a good place in the code to do the interpolation
    wfs=[]

    jjjj=0;p1=np.loadtxt('/data/users/mklymenko/science/H2/programing/15d2/'+'results/ff_'+str(jjjj)+'.dat')

    for j in a:
        print(j)
        wfs.append(WF(p1,Nst=j))

    b=specialCombinations(a)
    j1=0;i=[]

    i=[(j, comp_int(wfs[j[0]],wfs[j[1]],wfs[j[2]],wfs[j[3]])) for j in b]

    return i

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

if __name__ == "__main__":

    #Matrix formation:



    y = np.linspace(-3.0, 3.0, 300)
    x = np.linspace(0.0, 9.0, 300)
    xi,yi = np.meshgrid(x,y)

    x, y = xi.flatten(), yi.flatten()
    z=0.0*x
    XX=np.vstack((x,y,z));

    wf=WF(Nst=0,Nsys=0,save='/data/users/mklymenko/work_py/mb_project/',mf=2)
    wf.draw_func_grid(XX.T,wf.show_func(XX),par='2d')

#     a=ConfSet(Norb=7,Nel=8,M=0)
#     print('Number of configurations:'+" "+str(len(a.conf_a)))
#     print('---------------------------')
#     print("| ind"+" |"+"   configuration   "+"|")
#     print('---------------------------')
#
#     for i in a.conf_a:
#         print("   "+str(i.index)+"    "+str(i.d_id))
#
#     print('---------------------------')
#
#     N=5 #number of basis functions
#
#     a=list(xrange(N))#this is a good place in the code to do the interpolation
#     wfs=[]
#
#     for j in a:
#         print(j)
#         wfs.append(WF(Nst=j,Nsys=0,save='/data/users/mklymenko/work_py/mb_project/',mf=2))




    # a=[1, 0, 1, 0, 0]
    # a=sorted(a)

    # for p in next_permutation(a):
    #     print p

#    BS=ConfSet(Nel=2,Norb=3,M=0)



