#!/usr/bin/python

import numpy as np
import sys


def norb(N):

    return (np.sqrt(1+8*N)-1)/2


def mat3ind(ind, N):

    ind = ind+1
    j = np.linspace(1, N, N)
    di = (2*N-j)/2*(j-1)+j
    i = ind - di
    i = len([np.where(i >= 0)])
    j = ind-di[i-1]+i

    return (i-1, j-1)


def ind3mat(i, j, N, str):
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
