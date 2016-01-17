function ind=ind3mat_v(j1,v1,j2,v2,Nbands,fl)

i=j1+Nbands*(v1-1);
j=j2+Nbands*(v2-1);

ind=ind3mat(i,j,6*Nbands,fl);