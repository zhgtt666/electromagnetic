function [D] = eal2d(vp,nlayer,dx,dz,NX,NZ)
%% 2D effective absorbing boundary condition
% INPUT
% vp: velocity of P wave
% nlayer: number of absorbing boundary layers
% dx: samping interval along x dimention
% dz: samping interval along z dimention
% NX: number of sampling points along z dimention
% NZ: number of sampling points along z dimention
% OUTPUT
% D: attenuation damping coefficient
% reference
% Yao, G., Da Silva, N. V., & Wu, D. (2018). An effective absorbing layer 
% for the boundary condition in acoustic seismic wave simulation. Journal 
% of Geophysics and Engineering, 15(2), 495-511.
%%
vp_max = max(max(vp));
Lx = nlayer*dx;
Lz = nlayer*dz;
R = 10^(-(log10(nlayer)-1)/log10(2)-3);
d0x = -log(R)*(3*vp_max)/(2*Lx);
d0z = -log(R)*(3*vp_max)/(2*Lz);
%%
ddz = zeros(NZ,1);
ddx = zeros(1,NX);
for i = 1:nlayer
    ddz(i) = d0z*((nlayer-i)/nlayer)^2;
    ddz(NZ-i+1) = ddz(i);
    ddx(i) = d0x*((nlayer-i)/nlayer)^2;
    ddx(NX-i+1) = ddx(i);
end
D = sqrt(ddx.^2+ddz.^2);
end