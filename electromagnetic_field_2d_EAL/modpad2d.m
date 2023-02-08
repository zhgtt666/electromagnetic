function [mod] = modpad2d(mod,nlayer,NZ,NX)
%% pading 2D model parameter for absorbing boundary condition
% INPUT
% mod: model parameter
% nlayer: number of absorbing boundary layers
% NX: number of sampling points along z dimention
% NZ: number of sampling points along z dimention
% OUTPUT
% mod: model parameter after pading boundary

[nz,nx] = size(mod);
modp = zeros(NZ,NX);
modp(nlayer+1:nz+nlayer,nlayer+1:nx+nlayer) = mod;
modp(1:nlayer,:) = modp(nlayer+1,:).*ones(nlayer,NX);
modp(nz+nlayer+1:NZ,:) = modp(nz+nlayer,:).*ones(nlayer,NX);
modp(:,1:nlayer) = modp(:,nlayer+1).*ones(NZ,nlayer);
modp(:,nx+nlayer+1:NX) = modp(:,nx+nlayer).*ones(NZ,nlayer);
mod = modp;

end