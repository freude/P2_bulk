


def mat3ind(ind,N):

%not finished yet

    j=1:N;

    di=(2*N-j)./2.*(j-1)+j;

    i=(ind-di);
    i=length(i(i>=0));
    j=ind-di(i)+i;

    return (i,j)
