% QUAD-ROTOR AERSP 304
% By : Nicholas Giampetro, Craig Stenstrom, Payton Glynn

clc;
clear;
close all;

% const values
g = 9.81 ; % m/s
m = 0.450 ; % kg
l = 0.225 ; % 
k = 2.98*10^-7 ;
b = 1.14*10^-6 ;
I = [ 4.85*10^-3 , 0, 0 ;
     0 , 4.85*10^-3 , 0 ;
     0 , 0 , 8.80*10^-3 ] ;

% Part A

OmegaH = sqrt((m*g)/(4*k));
r0 = 0 ;
t = linspace(0,1,1000);
options = odeset('reltol',1e-12,'abstol',1e-12);
[t,vz,z] = ode45(@(t,x) (4*k*(OmegaH+70*sin(0.5*pi*t)))/m, t, r0, options);





    

