#!/usr/bin/python

import numpy as np
import sys
import pdb


def el2int(kk,bas_fun):
    integ=np.trapz(np.trapz(np.trapz(bas_fun,dx=kk,axis=2),dx=kk,axis=1),dx=kk,axis=0)/(2*(np.pi**2))/1.0919
    return integ

def el2int_c(kk,bas_fun,R_tr):

    kX, kY, kZ = np.meshgrid(kk,kk,kk)

    Q=np.power(kX,2)+np.power(kY,2)+np.power(kZ,2)

    V = (1-np.cos(np.sqrt(Q)*R_tr))/Q
    V[np.where(Q==0)] = 0.5*R_tr**2

    integ=np.trapz(np.trapz(np.trapz(bas_fun*V,dx=kk,axis=2),dx=kk,axis=1),dx=kk,axis=0)/(2*np.power(np.pi,2))/1.0919

    return integ

def u_i(j1,j2,G,tab):

    s = tab.shape
    tab_i = 0

    for j in xrange(s[3]):

        if (tab[j1,j2,j,0]==-G[0])&(tab[j1,j2,j,1]==-G[1])&(tab[j1,j2,j,2]==-G[2]):
            tab_i=tab[j1,j2,j,3]+1j*tab[j1,j2,j,4]

    return tab_i

def mat3ind(ind,N):

    a=np.where(np.reshape(np.cumsum(\
                np.triu(np.ones((N,N))))-1,(N,N))*np.triu(np.ones((N,N)))==ind)
    print a
#    ind = ind+1
#    j = np.arange(N)+1
#
#    di = (2*N-j)/2*(j-1)+j
#
#    i = ind - di
#    i = len(i[np.where(i >= 0)])
#    j = ind-di[i-1]+i
    i = a[0]
    j = a[1]
    return (i.item(0), j.item(0))


def ind3mat(i,j,N,str):
    i = i+1
    j = j+1
    if (str == 'sym'):
        if (j >= i):
            ind = (2*N-i)/2*(i-1)+j
        else:
            ind = (2*N-j)/2*(j-1)+i
    else:
        if (j >= i):
            ind = (2*N-i)/2*(i-1)+j
        else:
            ind = float('nan')

    return ind-1

def ind3mat_v(j1,v1,j2,v2,Nbands,fl):
    i=j1+Nbands*v1
    j=j2+Nbands*v2
    ind = ind3mat(i, j, 6*Nbands, fl)
    return ind

def mat3ind_v(ind,Nbands):

    p1, q1 = mat3ind(ind,Nbands*6)

    val1 = (p1//Nbands)
    st_ind1 = p1 - Nbands*val1

    val2 = (q1//Nbands)
    st_ind2 = q1 - Nbands*val2
    return (st_ind1,val1,st_ind2,val2)

def k_inv(k):

    if (k==0):
        kk=1

    if (k==1):
        kk=0

    if (k==2):
        kk=3

    if (k==3):
        kk=2

    if (k==4):
        kk=5

    if (k==5):
        kk=4

    return kk


def norb(N):
    return (np.sqrt(1+8*N)-1)/2

def of(fn):
    # try to open files with computed single-particle spectrum

    try:
        with open(fn):
            E = np.loadtxt(fn)
    except IOError:
        sys.exit("Can not open the file"+fn)

    return E

if __name__ == "__main__":

#    print(ind3mat_v(0,5,0,5,1,'sym'))
  #  print(mat3ind_v(14, 1))
     #print(ind3mat_v(0,5,0,5,1,'sym'))
     print(mat3ind(0, 1))
