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
options = odeset('reltol',1e-12,'abstol',1e-12);
[t,z] = integrator(@(t,z) ((4/m*k)*(OmegaH+70*sin(0.5*pi*t))), [0,0], 0, 1, 1000);

 



    

function [t,r] = integrator(func, init, tStart, tEnd, pnts)
    t = linspace(tStart,tEnd,pnts);
    options = odeset('reltol',1e-12,'abstol',1e-12);
    [t,r] = ode45(@(t,r) odeFun( @(t,r)func, t, r), t, init, options);
 end

function    rDot = odeFun(func,t,r)
    r1 = r(1);
    r2 = r(2);
    
    r1dot = r2;
    r2dot = @(t,r)func;
  
    rDot = [r1dot; r2dot];
end