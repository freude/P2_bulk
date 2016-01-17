from ctypes import *

class Mkl_Fft:
    c_double_p = POINTER(c_double)

    def __init__(self,num_threads=8):
        self.dfti = cdll.LoadLibrary("mk2_rt.dll")
        self.dfti.MKL_Set_Num_Threads(num_threads)
        self.Create = self.dfti.DftiCreateDescriptor_d_md
        self.Commit = self.dfti.DftiCommitDescriptor
        self.ComputeForward = self.dfti.DftiComputeForward

    def fft(self,a):
        Desc_Handle = c_void_p(0)
        dims = (c_int*2)(*a.shape)
        DFTI_COMPLEX = c_int(32)
        rank = 2

        self.Create(byref(Desc_Handle), DFTI_COMPLEX, rank, dims )
        self.Commit(Desc_Handle)
        self.ComputeForward(Desc_Handle, a.ctypes.data_as(self.c_double_p) )


x = np.linspace(-10.5, 10.5, 700)
y = np.linspace(0.0, 12.9, 700)

XXX=np.vstack((x,x*0.0+6.45,x*0.0))

xi,yi = np.meshgrid(x,y)
x, y = xi.flatten(), yi.flatten()
z=x*0.0
XX=np.vstack((x,y,z))

wf=GFit(qn=0,sn=5,mf=2,num_fu=21,psave='/data/users/mklymenko/work_py/mb_project/',pload='/data/users/mklymenko/work_py/mb_project/')
wf.save()

#wf.draw_func(x,y,par='2d')
g=wf.show_func(XX)
#g1=wf.show_gf(XXX)

X,F=wf.read_from_file(5,0)
invdisttree = Invdisttree(X.T,F, leafsize=10, stat=1)
AA=invdisttree(XX.T, nnear=130, eps=0, p=1)

fig=plt.figure()

#for j in range(0,wf._num_fu):
#    plt.plot(XXX[0,:].T,g1[:,j])

#plt.plot(XXX[0,:].T,g)
#plt.plot(XXX[0,:].T,AA)

#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(xi,yi,g.reshape(xi.shape), cmap=cm.jet, linewidth=0.2)

plt.contour(xi, yi, -AA.reshape(xi.shape), colors='red')
plt.contour(xi, yi, -g.reshape(xi.shape), colors='blue')

plt.hold(True)
plt.show()