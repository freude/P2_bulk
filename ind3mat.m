function ind=ind3mat(i,j,N,str)


if strcmp(str,'sym')
    if (j>=i)
        ind=(2*N-i)/2*(i-1)+j;
    else
        ind=(2*N-j)/2*(j-1)+i;
    end;    
else
    if (j>=i)
        ind=(2*N-i)/2*(i-1)+j;
    else
        ind=NaN;
    end;    
end;