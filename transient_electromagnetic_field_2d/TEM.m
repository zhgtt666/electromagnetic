% Description: simulation of transient electromagnetic field
% Version: 1.0
% Autor: Chen Quyang
% Date: 2023-1-7
% LastEditors: Yao Gang
% LastEditTime: 2023-01-15
% Copyright (c) 2023 WaveTomo. All rights reserved. 
%%
close;clear all;clc;

%%% Parameter define
dz=1; 
dx=1;
nx=1101; 
nz=601;  
x=((1:nx)-501)*dx;
z=((1:nz)-1)*dz;
xlen=length(x);
zlen=length(z);

mu=4*pi*1e-7;
sigma0=1/300; % background conductivity
sigma=sigma0*ones(nz,nx); 
% Anomalous body
lefttop_iz=100; % index of upper left corner of abnormal body in z direction
lefttop_ix=800;% index of upper left corner of abnormal body in x direction
width=20; height=300;
sigma(lefttop_iz:lefttop_iz+height,lefttop_ix:lefttop_ix+width)=(1e3)*1/300;
% imagesc(x,z,sigma0);
% rectangle('Position',[x(lefttop_ix),z(lefttop_iz),width,height],'Linewidth',1,'LineStyle','-','EdgeColor','black')
%% Time parameter 
% t0=1.13*mu*mean(sigma(:))*dx^2;
t0=1e-7;

nt1=1000;
nt2=10000;
nt3=9000;
dt1=1e-8;
dt2=1e-7;
dt3=5e-7;

if dt1 > sqrt(mu*min(sigma(:)))*min([dx dz])/2
    errordlg('dt1 is inappropriate!');
end

if dt2 > sqrt(mu*min(sigma(:)))*min([dx dz])/2
    errordlg('dt2 is inappropriate!');
end

if dt3 > sqrt(mu*min(sigma(:)))*min([dx dz])/2
    errordlg('dt3 is inappropriate!');
end
% Source define
I1=1;
I2=-1;
xs1=0; %The units is in meters
xs2=-400; 
zs1=0;
zs2=0;

% Initial condition 
E=zeros(zlen,xlen);
E0=zeros(zlen,xlen);
E1=zeros(zlen,xlen);
E2=zeros(zlen,xlen);
for ix=1:xlen
    for iz=2:zlen
        xpos1=x(1,ix)-xs1;
        zpos1=z(1,iz)-zs1;
        xpos2=x(1,ix)-xs2;
        zpos2=z(1,iz)-zs2;
        E0(iz,ix)=TEM_half_space_response(xpos1,zpos1,t0,sigma0,I1)...
                 +TEM_half_space_response(xpos2,zpos2,t0,sigma0,I2);
        E1(iz,ix)=TEM_half_space_response(xpos1,zpos1,t0+dt1,sigma0,I1)...
                 +TEM_half_space_response(xpos2,zpos2,t0+dt1,sigma0,I2);
    end
end

% Top boundary computation: method 1 - field continuation in wavenumber domain
xArrayTmp1=zeros(1,xlen+800);
kx=xwavenumber(xlen+800,dx);

xArrayTmp1(1,401:400+xlen)=E0(2,:);
xArrayTmp2=ifft(exp(-dz*abs(kx)).*fft(xArrayTmp1));
E0(1,:) =xArrayTmp2(1,401:400+xlen);

xArrayTmp1(1,401:400+xlen)=E1(2,:);
xArrayTmp2=ifft(exp(-dz*abs(kx)).*fft(xArrayTmp1));
E1(1,:) =xArrayTmp2(1,401:400+xlen);

% Top boundary computation: method 2 - field continuation using integeration
%  Note: method1 is more accurate than method 2.
%  E0(1,:) = 0;
%  for ix=1:xlen
%       for jx=1:xlen
%             E0(1,ix)= E0(1,ix) + (dz/pi)*E0(2,jx)*dx/((x(ix)-x(jx))^2+dz^2);
%        end
%  end
%      
%  E1(1,:) = 0;
%  for ix=1:xlen
%       for jx=1:xlen
%             E1(1,ix)= E1(1,ix) + (dz/pi)*E1(2,jx)*dx/((x(ix)-x(jx))^2+dz^2);
%        end
%  end
     
% E1(1,:) = E1(2,:); %% test Riemann boundary conditions
figure;
imagesc(E1);

% Save as GIF
filename = 'TEM.gif';
snpshotDelayTime = 0.2;
sliceInterval = 100;
%% Simulation part 1
tic;
% Function coefficient
A1=2*dt1./(mu*sigma);
A2=1./(1+A1/dz^2+A1/dx^2);

t_start=t0;
dt=dt1;
for it=1:nt1
    t = t_start+it*dt;
    % Calculate Es
    for ix=1:xlen
            xpos1=x(ix)-xs1; zpos1=z(zlen)-zs1;
            xpos2=x(ix)-xs2; zpos2=z(zlen)-zs2;
            E(zlen,ix) = TEM_half_space_response(xpos1,zpos1,t,sigma0,I1)+...
                       TEM_half_space_response(xpos2,zpos2,t,sigma0,I2);
    end
    for iz=1:zlen
            xpos1=x(1)-xs1; zpos1=z(iz)-zs1;
            xpos2=x(1)-xs2; zpos2=z(iz)-zs2;
            E(iz,1) = TEM_half_space_response(xpos1,zpos1,t,sigma0,I1)+...
                       TEM_half_space_response(xpos2,zpos2,t,sigma0,I2);
            xpos1=x(xlen)-xs1;zpos1=z(iz)-zs1;
            xpos2=x(xlen)-xs2;zpos2=z(iz)-zs2;
            E(iz,xlen) = TEM_half_space_response(xpos1,zpos1,t,sigma0,I1)+...
                       TEM_half_space_response(xpos2,zpos2,t,sigma0,I2);
    end
    
    % Top boundary computation: method 1 - field continuation in wavenumber domain
    xArrayTmp1(1,401:400+xlen)=E1(2,:);
    xArrayTmp2=ifft(exp(-dz*abs(kx)).*fft(xArrayTmp1));
    E1(1,:) =xArrayTmp2(1,401:400+xlen);

    % Bottom, left, right boundary 
    E1(zlen,:)=E(zlen,:);
    E1(:,1)=E(:,1);
    E1(:,xlen)=E(:,xlen);
    
    % Finite difference
    ix=2:xlen-1;
    iz=2:zlen-1;
    E2(iz,ix)= A2(iz,ix).*(A1(iz,ix).*((E1(iz+1,ix)-E0(iz,ix)+E1(iz-1,ix))/dz^2 ...
                      +(E1(iz,ix+1)-E0(iz,ix)+E1(iz,ix-1))/dx^2)+ E0(iz,ix));                 
  
    % Show result
    if mod(it,sliceInterval)==0 || it==1
        temp=E1*1e6;
        E_log=sign(temp).*log(abs(temp)+1);
        imagesc(x,z,E_log);
        rectangle('Position',[x(lefttop_ix),z(lefttop_iz),width,height],'Linewidth',1,'LineStyle','-','EdgeColor','black')
        title(['Ey at ',num2str(t*1e3),' ms']);
        % caxis([-1e5,1e5]);
        colormap jet; colorbar;
        set(get(colorbar,'title'),'string','log(E) (\muV/m)');
        set(gca,'YDir','reverse');
        xlabel('X / m');
        ylabel('Depth / m');
        max_value=max(max(log(abs(temp)+1)));
%         max_value=13;
        caxis([-max_value,max_value]);
        pause(0.001);
        if it==1
            savegif(it,filename,sliceInterval,snpshotDelayTime,1);
        else
            savegif(it,filename,sliceInterval,snpshotDelayTime,0);
        end
    end
    disp(it);
    if it==nt1-round(dt2/dt1)
        E0_pre=E1;
    end
    E0=E1;
    E1=E2;

end 

%%  Simulation part 2
sliceInterval = 200;
% Function coefficient
A1=2*dt2./(mu*sigma);
A2=1./(1+A1/dz^2+A1/dx^2);

E1=E0;
E0=E0_pre;
t_start=t0+nt1*dt1;
dt=dt2;
for it=1:nt2
    t = t_start+(it-1)*dt;
    % Calculate Es
    for ix=1:xlen
            xpos1=x(ix)-xs1; zpos1=z(zlen)-zs1;
            xpos2=x(ix)-xs2; zpos2=z(zlen)-zs2;
            E(zlen,ix) = TEM_half_space_response(xpos1,zpos1,t,sigma0,I1)+...
                       TEM_half_space_response(xpos2,zpos2,t,sigma0,I2);
    end
    for iz=1:zlen
            xpos1=x(1)-xs1; zpos1=z(iz)-zs1;
            xpos2=x(1)-xs2; zpos2=z(iz)-zs2;
            E(iz,1) = TEM_half_space_response(xpos1,zpos1,t,sigma0,I1)+...
                       TEM_half_space_response(xpos2,zpos2,t,sigma0,I2);
            xpos1=x(xlen)-xs1;zpos1=z(iz)-zs1;
            xpos2=x(xlen)-xs2;zpos2=z(iz)-zs2;
            E(iz,xlen) = TEM_half_space_response(xpos1,zpos1,t,sigma0,I1)+...
                       TEM_half_space_response(xpos2,zpos2,t,sigma0,I2);
    end
    
    % Top boundary computation: method 1 - field continuation in wavenumber domain
    xArrayTmp1(1,401:400+xlen)=E1(2,:);
    xArrayTmp2=ifft(exp(-dz*abs(kx)).*fft(xArrayTmp1));
    E1(1,:) =xArrayTmp2(1,401:400+xlen);

    % Bottom, left, right boundary 
    E1(zlen,:)=E(zlen,:);
    E1(:,1)=E(:,1);
    E1(:,xlen)=E(:,xlen);
    
    % Finite difference
    ix=2:xlen-1;
    iz=2:zlen-1;
    E2(iz,ix)= A2(iz,ix).*(A1(iz,ix).*((E1(iz+1,ix)-E0(iz,ix)+E1(iz-1,ix))/dz^2 ...
                      +(E1(iz,ix+1)-E0(iz,ix)+E1(iz,ix-1))/dx^2)+ E0(iz,ix));                 
  
    % Show result
    if mod(it,sliceInterval)==0
        temp=E1*1e6;
        E_log=sign(temp).*log(abs(temp)+1);
        imagesc(x,z,E_log);
        rectangle('Position',[x(lefttop_ix),z(lefttop_iz),width,height],'Linewidth',1,'LineStyle','-','EdgeColor','black')
        title(['Ey at ',num2str(t*1e3),' ms']);
        % caxis([-1e5,1e5]);
        colormap jet; colorbar;
        set(get(colorbar,'title'),'string','log(E) (\muV/m)');
        set(gca,'YDir','reverse');
        xlabel('X / m');
        ylabel('Depth / m');
        max_value=max(max(log(abs(temp)+1)));
%         max_value=13;
        caxis([-max_value,max_value]);

        pause(0.001);
        savegif(it,filename,sliceInterval,snpshotDelayTime,0);
    end
    disp(it);
    if it==nt2-round(dt3/dt2)
        E0_pre=E1;
    end
    E0=E1;
    E1=E2;
end 

%%  Simulation part 3
sliceInterval = 400;
% Function coefficient
A1=2*dt3./(mu*sigma);
A2=1./(1+A1/dz^2+A1/dx^2);

E1=E0;
E0=E0_pre;
t_start=t0+nt1*dt1+(nt2-1)*dt2;
dt=dt3;
for it=1:nt3
    t = t_start+(it-1)*dt;
    % Calculate Es
    for ix=1:xlen
            xpos1=x(ix)-xs1; zpos1=z(zlen)-zs1;
            xpos2=x(ix)-xs2; zpos2=z(zlen)-zs2;
            E(zlen,ix) = TEM_half_space_response(xpos1,zpos1,t,sigma0,I1)+...
                       TEM_half_space_response(xpos2,zpos2,t,sigma0,I2);
    end
    for iz=1:zlen
            xpos1=x(1)-xs1; zpos1=z(iz)-zs1;
            xpos2=x(1)-xs2; zpos2=z(iz)-zs2;
            E(iz,1) = TEM_half_space_response(xpos1,zpos1,t,sigma0,I1)+...
                       TEM_half_space_response(xpos2,zpos2,t,sigma0,I2);
            xpos1=x(xlen)-xs1;zpos1=z(iz)-zs1;
            xpos2=x(xlen)-xs2;zpos2=z(iz)-zs2;
            E(iz,xlen) = TEM_half_space_response(xpos1,zpos1,t,sigma0,I1)+...
                       TEM_half_space_response(xpos2,zpos2,t,sigma0,I2);
    end
    
    % Top boundary computation: method 1 - field continuation in wavenumber domain
    xArrayTmp1(1,401:400+xlen)=E1(2,:);
    xArrayTmp2=ifft(exp(-dz*abs(kx)).*fft(xArrayTmp1));
    E1(1,:) =xArrayTmp2(1,401:400+xlen);
     
    % Bottom, left, right boundary 
    E1(zlen,:)=E(zlen,:);
    E1(:,1)=E(:,1);
    E1(:,xlen)=E(:,xlen);
    
    % Finite difference
    ix=2:xlen-1;
    iz=2:zlen-1;
    E2(iz,ix)= A2(iz,ix).*(A1(iz,ix).*((E1(iz+1,ix)-E0(iz,ix)+E1(iz-1,ix))/dz^2 ...
                      +(E1(iz,ix+1)-E0(iz,ix)+E1(iz,ix-1))/dx^2)+ E0(iz,ix));                 
  
    % Show result
    if mod(it,sliceInterval)==0
        temp=E1*1e6;
        E_log=sign(temp).*log(abs(temp)+1);
        imagesc(x,z,E_log);
        rectangle('Position',[x(lefttop_ix),z(lefttop_iz),width,height],'Linewidth',1,'LineStyle','-','EdgeColor','black')
        title(['Ey at ',num2str(t*1e3),' ms']);
        % caxis([-1e5,1e5]);
        colormap jet; colorbar;
        set(get(colorbar,'title'),'string','log(E) (\muV/m)');
        set(gca,'YDir','reverse');
        xlabel('X / m');
        ylabel('Depth / m');
        max_value=max(max(log(abs(temp)+1)));
%         max_value=13;
        caxis([-max_value,max_value]);
        pause(0.001);
        savegif(it,filename,sliceInterval,snpshotDelayTime,0);
    end
    disp(it);
    if it==nt3-10
        E0_pre=E1;
    end
    E0=E1;
    E1=E2;
end 
toc;