#!/usr/bin/python

import numpy as np


def el2int(kk,bas_fun):

    integ=np.trapz(np.trapz(np.trapz(bas_fun,dx=kk,axis=2),dx=kk,axis=1),dx=kk,axis=0)/(2*np.power(pi,2))/1.0919
    return integ


def mat3ind(ind,N):


    j=1:N;

    di=(2*N-j)./2.*(j-1)+j;

    i=(ind-di);
    i=length(i(i>=0));
    j=ind-di(i)+i;

    return (i,j)
