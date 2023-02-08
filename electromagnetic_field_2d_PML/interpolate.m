function [mod] = interpolate(mod,direction,nz,nx)
%interpolate 
% INPUT
% mod: model parameter
% direction: =1, z direction; =2, x direction 
% nz: z direction number
% nx: x direction number
switch direction
    case 1
        z=1:nz;
        zi=1:0.5:nz+1;
        for i=1:nx
            yi = interp1(z,mod(:,i),zi,'Spline');
            mod(:,i) = yi(2:2:end);
        end
    case 2
        x=1:nx;
        xi=1:0.5:nx+1;
        for i=1:nz
            yi = interp1(x,mod(i,:),xi,'Spline');
            mod(i,:) = yi(2:2:end);
        end
end
% OUTPUT
% mod: model parameter after interpolating 
end

