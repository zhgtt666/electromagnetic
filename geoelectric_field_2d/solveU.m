function [U] = solveU(N,M,h,sigma1,sigma,I1,spos1,I2,spos2)
% INPUT
% N: number of sampling points along z dimention
% M: number of sampling points along x dimention
% h: samping interval along h dimention
% sigma1: earth conductivity
% sigma: conductivity
% I1: source 1 current intensity
% spos1: Location of the source 1, a binary array of the positions of x and z
% I2: source 2 current intensity
% spos2: Location of the source 2, a binary array of the positions of x and z
% OUTPUT
% U: Electrodynamic Potential
% 参考文献：邓正栋 等，2001，《稳定地电场三维有限差分正演模拟》，石油物探，40（1）
x=((1:M)-100)*h;
z=(1:N)*h;
xs1=x(spos1(1));
zs1=z(spos1(2));
xs2=x(spos2(1));
zs2=z(spos2(2));
%% Calculate A
% --- b ------ a ----- d ------- 0 
%     |        |       |         |  
%     |        |       |         |        
% --- c ------ b ----- a ------- 0
%     |        |       |         |        
%     |        |       |         |        
% --- 0 ------ c------ b ------- a 
%     |        |       |         |       
%     |        |       |         |        
% --- 0 ------ e ----- c ------- b
 
% Initialize factor
i=2:N-1; j=2:M-1;
a=(sigma(i+1,j)-sigma(i-1,j))/(4*h^2) + sigma(i,j)/h^2;
b=-4*sigma(i,j)/h^2;
c=-(sigma(i+1,j)-sigma(i-1,j))/(4*h^2) + sigma(i,j)/h^2;
d=(sigma(i,j+1)-sigma(i,j-1))/(4*h^2) + sigma(i,j)/h^2;
e=-(sigma(i,j+1)-sigma(i,j-1))/(4*h^2) + sigma(i,j)/h^2;
% Join b
vlen=(N-2)*(M-2)+2*(N-2)*(M-3)+2*((N-2)*(M-2)-1);
iv=1;i=zeros(1,vlen);j=zeros(1,vlen);v=zeros(1,vlen);
for k=1:(N-2)*(M-2)
    i(iv)=k;
    j(iv)=k;
    v(iv)=b(k);
    iv=iv+1;
end
% Calculate d and e
for k=1:(N-2)*(M-3)
    % Join d
    id=k;
    i(iv)=id;
    j(iv)=N-2+id;
    v(iv)=d(id);
    iv=iv+1;
    % Join e
    ie=(N-2)*(M-2)+1-k;
    i(iv)=ie;
    j(iv)=ie-(N-2);
    v(iv)=e(ie);
    iv=iv+1;
end
% Calculate c
for k=2:(N-2)*(M-2)
    i(iv)=k;
    j(iv)=k-1;
    v(iv)=c(k);
    if mod(k,N-2)==1
        v(iv)=0;
    end
    iv=iv+1;
end
% Calculate a
for k=1:(N-2)*(M-2)-1
    i(iv)=k;
    j(iv)=k+1;
    v(iv)=a(k);
    if mod(k,N-2)==1
        v(iv)=v(iv)+c(k);
    elseif mod(k,N-2)==0
        v(iv)=0;
    end
    iv=iv+1;
end
A=sparse(i,j,v);
%% Calculate F
F=zeros((N-2)*(M-2),1);
F(spos1(1)*(N-2)+spos1(2)-1,1)=I1;
F(spos2(1)*(N-2)+spos2(2)-1,1)=I2;

% NOTE: There are still problems with analytical solutions as boundary
% conditions !!! The analytical solution value is too big. 
% So when this problem is fixed, uncomment the follow lines.

for k=1:(N-2)*(M-2)
    if mod(k,N-2)==0
        F(k)=F(k)-a(N-2,k/(N-2))*(u_solution(x(mod(k,N-2)+1),z(N),xs1,zs1,sigma1,I1)...
                                 +u_solution(x(mod(k,N-2)+1),z(N),xs2,zs2,sigma1,I2));
    end
    if k<=N-2
         F(k)=F(k)-e(k,1)*(u_solution(x(1),z(k+1),xs1,zs1,sigma1,I1)...
                                 +u_solution(x(1),z(k+1),xs2,zs2,sigma1,I2));
    elseif k>(N-2)*(M-3)
        F(k)=F(k)-d(k-(N-2)*(M-3),M-2)*(u_solution(x(M),z(k-(N-2)*(M-3)+1),xs1,zs1,sigma1,I1)...
                                 +u_solution(x(M),z(k-(N-2)*(M-3)+1),xs2,zs2,sigma1,I2));
    end
end

%% Solve AU=F
U=A\F;
U=reshape(U,[N-2,M-2]);
end

