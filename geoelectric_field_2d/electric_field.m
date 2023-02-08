% Description: simulation of electrostatic field
% Version: 1.0
% Autor: Chen Quyang
% Date: 2023-1-19
% LastEditors: Chen Quyang
% LastEditTime: 2023-01-31
% Copyright (c) 2023 WaveTomo. All rights reserved. 

clear;clc;close all;

% Define parameter
M=802; % nx
N=202; % nz
h=0.01; % dx=dz=h 
sigma1=1/100;
sigma=sigma1*ones(N,M);
% Define source
I1=1;
I2=-1;
ix1=round((M-2)/4);
iz1=2;
ix2=round((M-2)/4*3);
iz2=2;
U=solveU(N,M,h,sigma1,sigma,I1,[ix1 iz1],I2,[ix2 iz2]);
%% Define low conductivity anomalous body
iax=round(M/2)+100;
iaz=60;
r=50;
for iz =1:N
    for ix = 1:M
        if((iz-iaz)^2+(ix-iax)^2<=r^2)
            sigma(iz,ix)=1/1000;
        end
    end
end
U1 = solveU(N,M,h,sigma1,sigma,I1,[ix1 iz1],I2,[ix2 iz2]);
% Define high conductivity anomalous body
for iz =1:N
    for ix = 1:M
        if((iz-iaz)^2+(ix-iax)^2<=r^2)
            sigma(iz,ix)=1/10;
        end
    end
end
U2 = solveU(N,M,h,sigma1,sigma,I1,[ix1 iz1],I2,[ix2 iz2]);
%% Show result
u(:,:,1)=U;
u(:,:,2)=U1;
u(:,:,3)=U2;
strlist=["Homogeneous","Low conductivity anomaly","High conductivity anomaly"];
fig=figure(1);set(fig,'position',[100 100 1000 600]);
x=((1:M-2)-round((M-2)/2))*h;
z=(1:N-2)*h;
pos=[x(iax-r),z(iaz-r),2*r*h,2*r*h];
for i=1:3
    subplot(3,1,i);
    contourf(x,z,u(:,:,i),100);
    if i~=1
        rectangle('Position',pos,'Curvature',[1 1],'EdgeColor','r','LineWidth',3)
    end
    title(strlist(i));
    colormap jet; colorbar;
    set(get(colorbar,'title'),'string','V');
    set(gca,'YDir','reverse');
    xlabel('X / m');
    ylabel('Depth / m');
end
saveas( 1, 'electric_field.jpg');
%% Plot Surface Electrodynamic Potential
figure(2);
plot(x,U(1,:));
hold on;
plot(x,U1(1,:));
hold on;
plot(x,U2(1,:));
title("Surface voltage measurement");
legend(strlist,'Location','NorthWest');
xlabel('X / m');
ylabel('Voltage / V');
grid on;
saveas( 2, 'surface_voltage.jpg');