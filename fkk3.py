#!/usr/bin/python

import scipy.io
import os
from coordsys import CoordSys
from scipy.interpolate import griddata
import numpy as np
from const import MyConst
import array as array
from aux_func1 import *
# -----------------------------------------------------------
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

@profile
def fkk3():
    indi=17
    sav_e='no'

    cx = np.arange(25)
    cx=MyConst.a_Si/MyConst.ab*cx*0.5

    pth = os.path.dirname(os.getcwd())
    EigVec1 = scipy.io.loadmat(pth+'/dis_scr/M'+str(indi)+'.mat')
    kk = EigVec1['kk']
    Nbands = np.squeeze(EigVec1['Nbands'])
    bands = EigVec1['bands']

    Nbands=1

    EEE = EigVec1['a1']/40
    EigVec1 = EigVec1['EigVec1']

    cs = CoordSys(num_cells=80,arrays_sizes=2,units='au')
    cs.set_origin_cells(41)
    x=cs.x()
    X,Y,Z = np.meshgrid(cs.x(),cs.x(),cs.x())
    Xshp=X.shape
    #X=np.vstack((X.flatten(), Y.flatten(), Z.flatten())).T
    #Y=[]
    #Z=[]

    fl = scipy.io.loadmat(pth+'/dis_scr/fl'+str(indi)+'.mat')
    fl = fl['fl']

    M1 = np.zeros((Nbands,3,cs.coord_sizes,cs.coord_sizes,cs.coord_sizes))

    for jj1 in xrange(3):

        p1 = np.loadtxt(pth+'/dis_scr/v'+str(2-jj1)+'/ff_'+str(indi)+'.dat')

        for jj3 in xrange(Nbands):

            print(jj3)

            a = []
            a = np.where(
                ((p1[:, 0] == 111) & (p1[:, 1] == 111) & (p1[:, 2] == 111)))[0]

            X1 = np.array(p1[a[jj3]+1:a[jj3+1], 0])
            Y1 = np.array(p1[a[jj3]+1:a[jj3+1], 1])
            Z1 = np.array(p1[a[jj3]+1:a[jj3+1], 2])
            F1 = np.array(p1[a[jj3]+1:a[jj3+1], 3])
            X1=np.vstack((X1.flatten(), Y1.flatten(), Z1.flatten())).T
            Y1=[]
            Z1=[]

            #invdisttree = Invdisttree(X1, F1, leafsize=10, stat=1)
            #M = invdisttree(X, nnear=30, eps=0, p=1)

            M = griddata(X1, F1, (X, Y, Z), method='linear')
            M[np.isnan(M)] = 0

            #fig = plt.figure()
            #ax = fig.add_subplot(111, projection='3d')
            #ax.plot_surface(X[:,:,81],Y[:,:,81],M[:,:,81], cmap=cm.jet, linewidth=0.2)
            #plt.hold(True)
            #plt.show()
            M1[jj3,jj1,:,:,:] = fl[0,jj1]*M

    #############################################################

    X = []
    Y = []
    Z = []
    X1 = []
    Y1 = []
    Z1 = []
    M = []

    #############################################################

    ap = 200
    Rc = 5.9

    x1 = np.fft.fftshift(np.fft.fftfreq(cs.coord_sizes+ap, x[3]-x[2]))

    k_tr = 0
    c2 = np.argmin(np.abs(x1-k_tr))
    c1 = np.argmin(np.abs(x1+k_tr))
    x1 = x1[c1:c2]

    Ff = np.zeros((Nbands,6,Nbands,6,Nbands,6,Nbands,6))
    Ff_c = np.zeros((Nbands,6,Nbands,6,Nbands,6,Nbands,6))

    for j1 in xrange(ind3mat_v(Nbands,6,Nbands,6,Nbands,'sym')):
        print(j1)
        st1, v1, st2, v2 = mat3ind_v(j1,Nbands)
        for st3 in xrange(Nbands):
            for st4 in xrange(Nbands):
                if v1==v2:
                    for v3 in xrange(6):

                        j2=ind3mat_v(st3,v3,st4,v3,Nbands,'nonsym')

                        if (not np.isnan(j2))&(Ff[st1,v1,st2,v2,st3,v3,st4,v3]==0):

                            Ff1 = np.abs((np.fft.ifftshift(np.fft.fftn(np.fft.fftshift(np.pad(\
                                    np.squeeze(M1[st1,v1//2,:,:,:])*np.squeeze(M1[st2,v2//2,:,:,:]),\
                                    ((ap//2, ap//2), (ap//2, ap//2), (ap//2, ap//2)),'constant'))))))*np.power(x[3]-x[2],3)

                            Ff2 = np.abs((np.fft.ifftshift(np.fft.fftn(np.fft.fftshift(np.pad(\
                                    np.squeeze(M1[st3,v3//2,:,:,:])*np.squeeze(M1[st4,v3//2,:,:,:]),\
                                    ((ap//2, ap//2), (ap//2, ap//2), (ap//2, ap//2)),'constant'))))))*np.power(x[3]-x[2],3)

                            Ff[st1,v1,st2,v2,st3,v3,st4,v3]=el2int(x1,Ff1[c2:-1:c1,c2:-1:c1,c2:-1:c1]*Ff2[c1:c2,c1:c2,c1:c2])

                            Ff[st3,v3,st4,v3,st1,v1,st2,v1]=Ff[st1,v1,st2,v1,st3,v3,st4,v3]
                            Ff[st2,v1,st1,v1,st4,v3,st3,v3]=Ff[st1,v1,st2,v1,st3,v3,st4,v3]
                            Ff[st4,v3,st3,v3,st2,v1,st1,v1]=Ff[st1,v1,st2,v1,st3,v3,st4,v3]
                            Ff[st2,v1,st1,v1,st3,v3,st4,v3]=Ff[st1,v1,st2,v1,st3,v3,st4,v3]
                            Ff[st4,v3,st3,v3,st1,v1,st2,v1]=Ff[st1,v1,st2,v1,st3,v3,st4,v3]
                            Ff[st1,v1,st2,v1,st4,v3,st3,v3]=Ff[st1,v1,st2,v1,st3,v3,st4,v3]
                            Ff[st3,v3,st4,v3,st2,v1,st1,v1]=Ff[st1,v1,st2,v1,st3,v3,st4,v3]

                            Ff_c[st1,v1,st2,v1,st3,v3,st4,v3]=el2int_c(x1,Ff1[c2:-1:c1,c2:-1:c1,c2:-1:c1]*Ff2[c1:c2,c1:c2,c1:c2],Rc)

                            Ff_c[st3,v3,st4,v3,st1,v1,st2,v1]=Ff_c[st1,v1,st2,v1,st3,v3,st4,v3]
                            Ff_c[st2,v1,st1,v1,st4,v3,st3,v3]=Ff_c[st1,v1,st2,v1,st3,v3,st4,v3]
                            Ff_c[st4,v3,st3,v3,st2,v1,st1,v1]=Ff_c[st1,v1,st2,v1,st3,v3,st4,v3]
                            Ff_c[st2,v1,st1,v1,st3,v3,st4,v3]=Ff_c[st1,v1,st2,v1,st3,v3,st4,v3]
                            Ff_c[st4,v3,st3,v3,st1,v1,st2,v1]=Ff_c[st1,v1,st2,v1,st3,v3,st4,v3]
                            Ff_c[st1,v1,st2,v1,st4,v3,st3,v3]=Ff_c[st1,v1,st2,v1,st3,v3,st4,v3]
                            Ff_c[st3,v3,st4,v3,st2,v1,st1,v1]=Ff_c[st1,v1,st2,v1,st3,v3,st4,v3]

                else:
                    #                 j2=ind3mat_v(st3,v2,st4,v1,Nbands,'nonsym')
                    #                 j2
                    #                 if ~isnan(j2)
                    if (Ff[st1,v1,st2,v2,st3,v2,st4,v1]==0):

                        Ff1 = np.abs((np.fft.ifftshift(np.fft.fftn(np.fft.fftshift(np.pad(\
                            np.squeeze(M1[st1,v1//2,:,:,:])*np.squeeze(M1[st2,v2//2,:,:,:]),\
                            ((ap//2, ap//2), (ap//2, ap//2), (ap//2, ap//2)),'constant'))))))*np.power(x[3]-x[2],3)


                        Ff2 = np.abs((np.fft.ifftshift(np.fft.fftn(np.fft.fftshift(np.pad(\
                            np.squeeze(M1[st3,v2//2,:,:,:])*np.squeeze(M1[st4,v1//2,:,:,:]),\
                        ((ap//2, ap//2), (ap//2, ap//2), (ap//2, ap//2)),'constant'))))))*np.power(x[3]-x[2],3)


                        Ff[st1,v1,st2,v2,st3,v2,st4,v1] = el2int(x1,Ff1[c2:-1:c1,c2:-1:c1,c2:-1:c1]*Ff2[c1:c2,c1:c2,c1:c2])

                        Ff[st3,v2,st4,v1,st1,v1,st2,v2]=Ff[st1,v1,st2,v2,st3,v2,st4,v1]
                        Ff[st2,v2,st1,v1,st4,v1,st3,v2]=Ff[st1,v1,st2,v2,st3,v2,st4,v1]
                        Ff[st4,v1,st3,v2,st2,v2,st1,v1]=Ff[st1,v1,st2,v2,st3,v2,st4,v1]
                        Ff[st2,v2,st1,v1,st3,v2,st4,v1]=Ff[st1,v1,st2,v2,st3,v2,st4,v1]
                        Ff[st4,v1,st3,v2,st1,v1,st2,v2]=Ff[st1,v1,st2,v2,st3,v2,st4,v1]
                        Ff[st1,v1,st2,v2,st4,v1,st3,v3]=Ff[st1,v1,st2,v2,st3,v2,st4,v1]
                        Ff[st3,v2,st4,v1,st2,v2,st1,v1]=Ff[st1,v1,st2,v2,st3,v2,st4,v1]

                    #-----------------------------------------------------------------------

                    if (Ff[st1,v1,st2,v2,st3,v2,st4,v1]==0):

                        Ff1 = ((np.fft.ifftshift(np.fft.fftn(np.fft.fftshift(np.pad(\
                            np.squeeze(M1[st1,v1//2,:,:,:])*np.squeeze(M1[st2,v2//2,:,:,:]),\
                        ((ap//2, ap//2), (ap//2, ap//2), (ap//2, ap//2)),'constant'))))))*np.power(x[3]-x[2],3)

                        Ff2 = ((np.fft.ifftshift(np.fft.fftn(np.fft.fftshift(np.pad(\
                            np.squeeze(M1[st3,k_inv(v1)//2,:,:,:])*np.squeeze(M1[st4,k_inv(v2)//2,:,:,:]),\
                        ((ap//2, ap//2), (ap//2, ap//2), (ap//2, ap//2)),'constant'))))))*np.power(x[3]-x[2],3)

                        Ff[st1,v1,st2,v2,st3,v2,st4,v1] = el2int(x1,Ff1[c2:-1:c1,c2:-1:c1,c2:-1:c1]*Ff2[c1:c2,c1:c2,c1:c2])

                        Ff[st3,k_inv(v1),st4,k_inv(v2),st1,v1,st2,v2]=Ff[st1,v1,st2,v2,st3,k_inv(v1),st4,k_inv(v2)]
                        Ff[st2,v2,st1,v1,st4,k_inv(v2),st3,k_inv(v1)]=Ff[st1,v1,st2,v2,st3,k_inv(v1),st4,k_inv(v2)]
                        Ff[st4,k_inv(v2),st3,k_inv(v1),st2,v2,st1,v1]=Ff[st1,v1,st2,v2,st3,k_inv(v1),st4,k_inv(v2)]
                        Ff[st2,v2,st1,v1,st3,k_inv(v1),st4,k_inv(v2)]=Ff[st1,v1,st2,v2,st3,k_inv(v1),st4,k_inv(v2)]
                        Ff[st4,k_inv(v2),st3,k_inv(v1),st1,v1,st2,v2]=Ff[st1,v1,st2,v2,st3,k_inv(v1),st4,k_inv(v2)]
                        Ff[st1,v1,st2,v2,st4,k_inv(v2),st3,k_inv(v1)]=Ff[st1,v1,st2,v2,st3,k_inv(v1),st4,k_inv(v2)]
                        Ff[st3,k_inv(v1),st4,k_inv(v2),st2,v2,st1,v1]=Ff[st1,v1,st2,v2,st3,k_inv(v1),st4,k_inv(v2)]


    Ff = np.real(Ff)

    ##

    Gg = np.zeros((6,6,6,6,np.power(Nbands,4)))
    I1 = np.zeros((6,6,6,6,np.power(Nbands,4)))
    I2 = np.zeros((6,6,6,6,np.power(Nbands,4)))
    I3 = np.zeros((6,6,6,6,np.power(Nbands,4)))
    I4 = np.zeros((6,6,6,6,np.power(Nbands,4)))
    M_lim = np.zeros((6,6,6,6))

    Gg_c = np.zeros((6,6,6,6,np.power(Nbands,4)))
    I1_c = np.zeros((6,6,6,6,np.power(Nbands,4)))
    I2_c = np.zeros((6,6,6,6,np.power(Nbands,4)))
    I3_c = np.zeros((6,6,6,6,np.power(Nbands,4)))
    I4_c = np.zeros((6,6,6,6,np.power(Nbands,4)))
    M_lim_c = np.zeros((6,6,6,6))

    for v1 in xrange(6):
        for v2 in xrange(6):
            for v3 in xrange(6):
                for v4 in xrange(6):

                    A = np.squeeze(Ff[:,v1,:,v2,:,v3,:,v4])
                    B = np.sort(np.abs(A.ravel))
                    IND = np.argsort(np.abs(A.ravel))
                    B = B[::-1]
                    IND = IND[::-1]

                    C = B[np.where(B>0.1)]

                    I1[v1,v2,v3,v4,:], I2[v1,v2,v3,v4,:],\
                        I3[v1,v2,v3,v4,:], I4[v1,v2,v3,v4,:] =\
                        np.unravel_index(IND, (Nbands, Nbands, Nbands, Nbands))

                    Gg[v1,v2,v3,v4,:] = B

                    if not C:
                        M_lim[v1,v2,v3,v4] = 0
                    else:
                        M_lim[v1,v2,v3,v4] = len(C)


                    A=np.squeeze(Ff_c[:,v1,:,v2,:,v3,:,v4])
                    B = np.sort(np.abs(A.ravel))
                    IND = np.argsort(np.abs(A.ravel))
                    B = B[::-1]
                    IND = IND[::-1]

                    C = B[np.where(B>0.1)]

                    I1_c[v1,v2,v3,v4,:], I2_c[v1,v2,v3,v4,:],\
                        I3_c[v1,v2,v3,v4,:], I4_c[v1,v2,v3,v4,:] =\
                        np.unravel_index(IND, (Nbands, Nbands, Nbands, Nbands))

                    Gg_c[v1,v2,v3,v4,:] = B

                    if not C:
                        M_lim_c[v1,v2,v3,v4] = 0
                    else:
                        M_lim_c[v1,v2,v3,v4] = len(C)

    ###

    num_bs = 24
    ind_len = ind3mat(num_bs,num_bs,num_bs,'sym')

    num_el = 40
    tab, num_el = u_tab(num_el)
    integ = np.zeros((ind_len, ind_len))

    print("Large loop starts")

    for jj1 in xrange(ind_len):
        print "---"
        print jj1
        for jj2 in xrange(ind_len):
            if (jj2>=jj1):

                j1_1, j3_1 = mat3ind(jj1,num_bs)
                j1_2, j3_2 = mat3ind(jj2,num_bs)

                for v1 in xrange(6):
                    for v2 in xrange(6):

                        if (v1==v2):

                            for v3 in xrange(6):
                                if (M_lim[v1,v2,v3,v3]!=0):
                                    Ff_1=0
                                    Ff_2=0

                                    for jjj in xrange(np.max([M_lim(v1,v2,v3,v3), M_lim_c(v1,v2,v3,v3)])):
                                        Ff_1 = Ff_1 + fl[v1//2]*fl[v2//2]*fl[v3//2]*fl[v3//2]*\
                                                np.conj(EigVec1[I1[v1,v2,v3,v3,jjj]+Nbands*v1,j1_1])*\
                                                        EigVec1[I2[v1,v2,v3,v3,jjj]+Nbands*v2,j3_1]*\
                                                np.conj(EigVec1[I3[v1,v2,v3,v3,jjj]+Nbands*v3,j1_2])*\
                                                        EigVec1[I4[v1,v2,v3,v3,jjj]+Nbands*v3,j3_2]*Gg[v1,v2,v3,v3,jjj]

                                        Ff_2 = Ff_2 + fl[v1//2]*fl[v2//2]*fl[v3//2]*fl[v3//2]*\
                                                np.conj(EigVec1[I1_c[v1,v2,v3,v3,jjj]+Nbands*v1,j1_1])*\
                                                        EigVec1[I2_c[v1,v2,v3,v3,jjj]+Nbands*v2,j3_1]*\
                                                np.conj(EigVec1[I3_c[v1,v2,v3,v3,jjj]+Nbands*v3,j1_2])*\
                                                        EigVec1[I4_c[v1,v2,v3,v3,jjj]+Nbands*v3,j3_2]*Gg_c[v1,v2,v3,v3,jjj]


                                    for j1 in xrange(num_el[v1,v2]):
                                        if np.sqrt(tab[v1,v2,j1,1]**2+tab[v1,v2,j1,2]**2+tab[v1,v2,j1,3]**2)==0:
                                            u1 = tab[v1,v2,j1,4]+1j*tab[v1,v2,j1,5]
                                            integ[jj1,jj2] = integ[jj1,jj2]+u1*u1*Ff_2
                                        else:
                                            u1 = tab[v1,v2,j1,4]+1j*tab[v1,v2,j1,5]
                                            u2 = u_i(v3,v3,(tab[v1,v2,j1,1], tab[v1,v2,j1,2], tab[v1,v2,j1,3]),tab)
                                            Q = (kk[v2,1]-kk[v1,1]+tab[v1,v2,j1,1])**2+\
                                                (kk[v2,2]-kk[v1,2]+tab[v1,v2,j1,2])**2+\
                                                (kk[v2,3]-kk[v1,3]+tab[v1,v2,j1,3])**2

                                            V = (1-np.cos(np.sqrt(Q)*Rc))/Q

                                            integ[jj1,jj2] = integ[jj1,jj2]+u1*u2*Ff_1*V

                        else:
                            Ff_1=0
                            for jjj in xrange(M_lim[v1,v2,v2,v1]):
                                Ff_1 = Ff_1 + fl[v1//2]*fl[v2//2]*fl[v2//2]*fl[v1//2]*\
                                        np.conj(EigVec1[I1[v1,v2,v2,v1,jjj]+Nbands*v1,j1_1])*\
                                                EigVec1[I2[v1,v2,v2,v1,jjj]+Nbands*v2,j3_1]*\
                                        np.conj(EigVec1[I3[v1,v2,v2,v1,jjj]+Nbands*v2,j1_2])*\
                                                EigVec1[I4[v1,v2,v2,v1,jjj]+Nbands*v1,j3_2]*Gg[v1,v2,v2,v1,jjj]

                            for j1 in xrange( num_el[v1,v2]):
                                u1 = u_i(v2,v1,(tab[v1,v2,j1,1], tab[v1,v2,j1,2], tab[v1,v2,j1,3]),tab)
                                u2 = tab[v1,v2,j1,4]+1j*tab[v1,v2,j1,5]

                                Q = (kk[v2,1]-kk[v1,1]+tab[v1,v2,j1,1])**2+\
                                    (kk[v2,2]-kk[v1,2]+tab[v1,v2,j1,2])**2+\
                                    (kk[v2,3]-kk[v1,3]+tab[v1,v2,j1,3])**2

                                V = (1-np.cos(np.sqrt(Q)*Rc))/Q

                                integ[jj1,jj2] = integ[jj1,jj2]+u1*u2*Ff_1*V

                            Ff_1=0
                            for jjj in xrange(M_lim[v1,v2,k_inv(v1),k_inv(v2)]):
                                Ff_1 = Ff_1 + fl[v1//2]*fl[v2//2]*fl[k_inv(v1)//2]*fl[k_inv(v2)//2]*\
                                        np.conj(EigVec1[I1[v1,v2,k_inv(v1),k_inv(v2),jjj]+Nbands*v1,j1_1])*\
                                                EigVec1[I2[v1,v2,k_inv(v1),k_inv(v2),jjj]+Nbands*v2,j3_1]*\
                                        np.conj(EigVec1[I3[v1,v2,k_inv(v1),k_inv(v2),jjj]+Nbands*k_inv(v1),j1_2])*\
                                                EigVec1[I4[v1,v2,k_inv(v1),k_inv(v2),jjj]+Nbands*k_inv(v2),j3_2]*\
                                        Gg[v1,v2,k_inv(v1),k_inv(v2),jjj]

                            for j1 in xrange(num_el[v1,v2]):
                                u1=u_i(k_inv(v1),k_inv(v2),(tab[v1,v2,j1,1], tab[v1,v2,j1,2], tab[v1,v2,j1,3]),tab)
                                u2 = tab[v1,v2,j1,4]+1j*tab[v1,v2,j1,5]

                                Q = (kk[v2,1]-kk[v1,1]+tab[v1,v2,j1,1])**2+\
                                    (kk[v2,2]-kk[v1,2]+tab[v1,v2,j1,2])**2+\
                                    (kk[v2,3]-kk[v1,3]+tab[v1,v2,j1,3])**2

                                V = (1-np.cos(np.sqrt(Q)*Rc))/Q
                                integ[jj1,jj2] = integ[jj1,jj2]+u1*u2*Ff_1*V

            else:
                integ[jj1, jj2] = integ[jj2, jj1]

    ###

    exch = np.zeros((num_bs,num_bs))

    for j1 in xrange(num_bs):
        for j2 in xrange(num_bs):
            exch[j1,j2] = 0
            for j3 in xrange(num_bs):
                jj1 = ind3mat(j1,j3,num_bs,'sym')
                jj2 = ind3mat(j3,j2,num_bs,'sym')
                exch[j1,j2]=exch[j1,j2]-0.5*np.abs(integ[jj1,jj2])

    if (sav_e == 'yes'):
        # save([pwd,'/dis_scr/bas_fun_ms',num2str(indi),'.mat'], 'bas_fun_ms', '-v7.3')
        np.savetxt(pth+'/dis_scr/!E'+str(indi)+'.dat', EEE, delimiter=' ')
        np.savetxt(pth+'/dis_scr/!int2e'+str(indi)+'.dat', np.abs(integ), delimiter=' ')
        np.savetxt(pth+'/dis_scr/!exch'+str(indi)+'.dat', np.abs(exch), delimiter=' ')

    return

fkk3()
