%% ========================================================================
% Code for 2D TE mode wave propagation with EAL absorbing boundary
% Written by: Chen Quyang 
% Date: 2023.1.1
%% ========================================================================

close; clear; clc;

%% Parameter definition
nx = 101;                                                                   
nz = 101; 
nlayer = 30;
dh = 0.06;
dt = 0.25e-9;
nt = 401;
NZ = nz+2*nlayer;
NX = nx+2*nlayer;
x = (0:nx-1)*dh;
z = (0:nz-1)*dh;
t = (0:nt-1)*dt;

% Recviver definition
rx = 1:nx;
rz = 1;
rec=zeros(nt,length(rx));

% Source definition
sz = 1;
sx = 1:nx;
t0 = 20;                                                                   
tau = 6; 
gaussian = @(it)exp(-0.5*((t0-it)/tau)^2);
fp = 200e6;                                   
t_delay = 1.5/fp;                                 
ricker = @(it)(1-2*(pi*fp*(it*dt-t_delay))^2)*exp(-(pi*fp*(it*dt-t_delay))^2);
w0=100e7;
a=0.93*w0;
pulse = @(it)(1.5e18*(it*dt)^2*exp(-a*it*dt)*sin(w0*it*dt)); 
src = ricker; % set source type

% Medium parameter
c = 2.99792458e8; % electromagnet wave velocity
sigma0 = 0; % 0 conductivity in a vacuum
mu0 = 4*pi*10^-7; % permeability of vacuum
epsilon0 = 1/(mu0*c^2); % dielectric constant of vacuum

% Model definition
sigma=(5e-4)*ones(nz,nx);% conductivity
mur = 1*ones(nz,nx); % relative permeability
epsilonr = 6*ones(nz,nx); % relative dielectric constant
mu = mur*mu0; % absolute permeability
epsilon = epsilonr*epsilon0; % absolute permeability

ci=[50,50,80];
cj=[floor(nx/2)-20,floor(nx/2)+20,floor(nx/2)];
r=3;
for k=1:3
    for i=1:nz
        for j=1:nx
            if (ci(k)-i)^2+(cj(k)-j)^2<=r^2
                sigma(i,j)=0.05;
                epsilon(i,j)=16;
            end
        end
    end
end

% Show Model
% figure;
% imagesc(x,z,sigma);
% title('Model');
% xlabel('Distance (m)');
% ylabel('Depth (m)');

[muz] = interpolate(mu,1,nz,nx);
[mux] = interpolate(mu,2,nz,nx);

[sigma] = modpad2d(sigma,nlayer,NZ,NX);
[muz] = modpad2d(muz,nlayer,NZ,NX);
[mux] = modpad2d(mux,nlayer,NZ,NX);
[epsilon] = modpad2d(epsilon,nlayer,NZ,NX);

% EAL absorbing boundary 
[D] = eal2d(c*ones(NZ,NX),nlayer,dh,dh,NX,NZ);
coef1 = (D*dt)/2;

% Function coefficient
b=(sigma*dt)./(2*epsilon);
A11=(1-b-coef1)./(1+b+coef1);
A12=(dt./epsilon)./((1+b+coef1)*dh);
A21=(1-coef1)./(1+coef1);
A22=dt./((mux*dh).*(1+coef1)); 
A31=A21;
A32=dt./((muz*dh).*(1+coef1)); 

% FD coefficient
c1 = 9/8; c2 = -1/24;

% Save snpshot
filename = 'snapshot.gif';
snpshotDelayTime = 0.1;% two snpshot interval time
slice = 5; % the number of time steps saved as a snpshot
ntrack = floor(nx/2); % the number of track  saved as a snpshot
%% Calculate 
for k = 1:nx
% Initialize wavefield
Ey = zeros(NZ,NX);
Hz = zeros(NZ,NX);
Hx = zeros(NZ,NX);
for it = 1:nt       
    % Calculate the Ey field
    for j = 3:NZ-1
        for i = 3:NX-1
            Ey(i,j) = A11(i,j)*Ey(i,j) + A12(i,j)*(c1*(Hx(i,j)-Hx(i-1,j))...
                                                +c2*(Hx(i+1,j)-Hx(i-2,j))...
                                                - (c1*(Hz(i,j)-Hz(i,j-1))...
                                                +c2*(Hz(i,j+1)-Hz(i,j-2))));
        end
    end

    % Calculate the Hz and Hx field
    for j = 2:NZ-2
        for i = 2:NX-2
            Hx(i,j) = A21(i,j)*Hx(i,j) + A22(i,j)*(c1*(Ey(i+1,j)-Ey(i,j))...
                                                 + c2*(Ey(i+2,j)-Ey(i-1,j)));
            Hz(i,j) = A31(i,j)*Hz(i,j) + A32(i,j)*(c1*(Ey(i,j)-Ey(i,j+1))...
                                                  + c2*(Ey(i,j-1)-Ey(i,j+2)));
        end
    end

    % Load source 
    Ey(nlayer+sz,nlayer+sx(k)) = Ey(nlayer+sz,nlayer+sx(k)) + src(it);
%     Ey(nlayer+sz,nlayer+sx) = Ey(nlayer+sz,nlayer+sx) + src(it);
    
    % Receive seismic record
    rec(it,sx(k))=Ey(nlayer+rz,nlayer+sx(k));
%     rec(it,rx)=Ey(nlayer+rz,nlayer+rx);

    % Show result
    if(k==ntrack)
        pcolor(x,z,Ey(nlayer+1:nlayer+nz,nlayer+1:nlayer+nx)); 
        shading interp;axis ij
        title(['Ey at ',num2str(dt*(it-1)*10^9),' ns']);
        xlabel('Distance (m)');
        ylabel('Depth (m)');
        colormap redblue;
        caxis([-0.01 0.01]);
        pause(0.01);

    
        frame = getframe(gcf);
        imind = frame2im(frame); 
        [imind, cmap] = rgb2ind(imind, 256); 
        if it == 1
            imwrite(imind,cmap, filename,'gif', 'Loopcount',inf,'DelayTime',snpshotDelayTime);
        else 
            if(mod(it,slice)==0)
            imwrite(imind,cmap, filename,'gif','WriteMode','append','DelayTime',snpshotDelayTime);
            end
        end
    end
end
disp(k);
end

%% Show seismic record 
figure;
imagesc(x,t*10^9,rec); 
title('receiver record');
xlabel('Distance (m)');
ylabel('time (ns)');
colormap gray;
caxis([-0.01 0.01])