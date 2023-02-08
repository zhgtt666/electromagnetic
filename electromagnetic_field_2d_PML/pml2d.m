function [ddz,ddx] = pml2d(v,nlayer,dx,dz,NX,NZ)
%% 2D perfectly matched layer condition
% INPUT
% v: velocity
% nlayer: number of absorbing boundary layers
% dx: samping interval along x dimention
% dz: samping interval along z dimention
% NX: number of sampling points along z dimention
% NZ: number of sampling points along z dimention
% OUTPUT
% ddz: attenuation damping coefficient along z dimention
% ddx: attenuation damping coefficient along x dimention
%%
v_max = max(max(v));
Lx = nlayer*dx;
Lz = nlayer*dz;
R = 10^(-(log10(nlayer)-1)/log10(2)-3);
d0x = -log(R)*(3*v_max)/(2*Lx);
d0z = -log(R)*(3*v_max)/(2*Lz);
%%
ddz = zeros(NZ,NX);
ddx = zeros(NZ,NX);
for i = 1:nlayer
    ddz(i,:) = d0z*((nlayer-i)/nlayer)^2;
    ddz(NZ-i+1,:) = ddz(i,:);
    ddx(:,i) = d0x*((nlayer-i)/nlayer)^2;
    ddx(:,NX-i+1) = ddx(:,i);
end
end