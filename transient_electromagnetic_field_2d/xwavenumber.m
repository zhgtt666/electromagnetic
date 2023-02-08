function [kx]=xwavenumber(nx,dx)
 %====================================================================
 %    [kx]=xwavenumber(nx,dx)
 %    create wavenumber kx
 %    nx: number of element along x (2nd) direction
 %    dx: interval along x (2nd) direction 
 %    return value:
 %    kx: kx array with size of 'nz x nx'   
 %====================================================================
    % kx
    deltF=1/(nx*dx);
    %==========
    wx=2*pi*deltF;
    index_max_wavenumber=fix(nx/2);
    k=zeros(1,nx);
    k(1,1:index_max_wavenumber+1)=0:index_max_wavenumber;%求0：N/2的波数
    for j=index_max_wavenumber+2:nx
       k(1,j)=-(nx-(j-1));%N/2+1:N-1的波数
    end
    kx=wx*k;
end