#!/usr/bin/python

import numpy as np
from const import MyConst
import sys
import scipy.io


def u_tab(s):

    k0 = MyConst.k0*MyConst.ab
    kk = k0*np.array([[1,  0,  0],\
        [-1,  0,  0],\
        [0,  1,  0],\
        [0, -1,  0],\
        [0,  0,  1],\
        [0,  0, -1]])

    G = np.zeros((6,6,s,6))
    num = np.zeros((6,6))

    for j2 in xrange(6):
        for j4 in xrange(6):
            wf = np.conj(abi_read(1, 70, kk[j2,:]))*(abi_read(1, 70, kk[j4,:]))
            G[j2,j4,:,:] = per_freq(wf,s)
            j = 0
            while ((0.05 < G[j2,j4,j,5]) and (j<s)):
                num[j2, j4] = j
                j = j + 1

    return (G, num)

def per_freq(wf,num=0):

    s = wf.shape
    lc = MyConst.a_Si/MyConst.ab
    dV = np.power((lc/s(0)),3)

    ub=np.fft.fftshift(np.fft.fftn(wf))*dV

    s = ub.shape
    x = np.arange(s[0])-s[0]//2
    X, Y, Z = np.meshgrid(x,x,x)
    G = np.zeros(((s[0]*s[1]*s[2]),6))

    G[:,0] = X.ravel*2*pi/lc
    G[:,1] = Y.ravel*2*pi/lc
    G[:,2] = Z.ravel*2*pi/lc
    G[:,3] = np.real(ub.ravel)
    G[:,4] = np.imag(ub.ravel)
    G[:,5] = np.power((np.abs(ub.ravel)),2)

    IX = np.argsort(np.squeeze(G[:,5]))
    IX = IX[::-1]
    G=G[IX,:]

    if (num!=0):
        G=G[0:num,:]

    return G

def abi_read(fac,T,valley):
    # the function reads the periodic functions computed by ABINIT
    # for fac number of unit cells

    wf1 = read_wf(T,valley) # read the wave function for a single unit cells

    # if valley(find(valley))<0
    #     wf1=-wf1;
    # end;

    # compose a fac number of cells
    wf = np.zeros((fac*T,fac*T,fac*T))

    for j1 in xrange(fac):
        for j2 in xrange(fac):
            for j3 in xrange(fac):
                wf[((j1-1)*T+1):(j1*T),((j2-1)*T+1):(j2*T),((j3-1)*T+1):(j3*T)] = wf1

    return wf

def read_wf(T,k1):

    # associate valley index with abinit wave-function

    if k1[0]!=0:

        if k1[0]>0:
            indi = 5

        if k1[0]<0:
            indi = 6

    if k1[1]~=0:

        if k1[1]>0:
            indi = 2

        if k1[1]<0:
            indi = 4

    if k1[2]~=0:

        if k1[2]>0:
            indi = 1

        if k1[2]<0:
            indi = 3

    # check if the unitcell function is already stored (both real and imaginary parts)

    pth = os.path.dirname(os.getcwd())
    pth1 = pth + '/dis_scr/wfr_' + str(indi) + '_' + str(T) + '.mat'
    pth2 = pth + '/dis_scr/wfi_' + str(indi) + '_' + str(T) + '.mat'

    try:
        with open(pth):
            print('(Data is saved)...')
            wf1 = scipy.io.loadmat(pth1)
            wf2 = scipy.io.loadmat(pth2)
            wf1 = wf1['wfr'] + 1j*wf2['wfi']
    except:
        pass

       # fprintf('(Data is not saved)...\n');

       # wf=ncread(strcat(pwd,'/dis_scr/tbase3_xo_DS2_PS_WFK-etsf.nc'),...
       #     'real_space_wavefunctions',[1 1 1 1 1 5 1 1],[Inf Inf Inf Inf 1 1 Inf 1],[1 1 1 1 1 1 1 1]);

       # sc=1;

       # wfr=squeeze(wf(1,1:sc:end,1:sc:end,1:sc:end,1,1,indi));
       # wfi=squeeze(wf(2,1:sc:end,1:sc:end,1:sc:end,1,1,indi));

       # wfr=transform_to_uc1(wfr,T);
       # wfi=transform_to_uc1(wfi,T);

       # save(strcat(pwd,['/dis_scr/wfr_',num2str(indi),'_',num2str(T),'.mat']), 'wfr');
       # save(strcat(pwd,['/dis_scr/wfi_',num2str(indi),'_',num2str(T),'.mat']), 'wfi');

       # wf1=wfr+ij.*wfi;

    kk = (MyConst.a_Si/MyConst.ab)/T
    ME = np.trapz(np.trapz(np.trapz(\
        np.power(np.abs(wf1),2),dx=kk,axis=2),dx=kk,axis=1),dx=kk,axis=0)/\
        (2*np.power(np.pi,2))/1.0919

    return (1/np.sqrt(ME))*wf1
