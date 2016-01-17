#!/usr/bin/python

import numpy as np
import sys

def el2int(kk,bas_fun):
    integ=np.trapz(np.trapz(np.trapz(bas_fun,dx=kk,axis=2),dx=kk,axis=1),dx=kk,axis=0)/(2*np.power(np.pi,2))/1.0919
    return integ


def mat3ind(ind,N):
    ind=ind+1
    j=np.arange(N)+1

    di=(2*N-j)/2*(j-1)+j

    i=(ind-di)
    i = len([np.where(i >= 0)])
    j=ind-di(i)+i

    return (i-1,j-1)


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

    print(mat3ind(1, 3))
