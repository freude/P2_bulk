function [st_ind1,val1,st_ind2,val2]=mat3ind_v(ind,Nbands)

[p1,q1]=mat3ind(ind,Nbands*6);
     
    val1=fix((p1-1)/Nbands)+1;
    st_ind1=p1-Nbands*(val1-1);
    
    val2=fix((q1-1)/Nbands)+1;
    st_ind2=q1-Nbands*(val2-1);