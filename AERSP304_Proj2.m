% QUAD-ROTOR AERSP 304
% By : Nicholas Giampetro, Craig Stenstrom, Payton Glynn

clc;
clear;
close all;

% const values
setGlobalG(9.81) ; % m/s
setGlobalM(0.450) ; % kg
setGlobalL(0.225) ; % 
setGlobalK(2.98*10^-7) ;
setGlobalB(1.14*10^-6) ;

I = [ 4.85*10^-3 , 0, 0 ;
     0 , 4.85*10^-3 , 0 ;
     0 , 0 , 8.80*10^-3 ] ;

% Part A

t = linspace(0,2,1000) ;

options = odeset('reltol',1e-12,'abstol',1e-12);
[t,z] = ode45(@(t,z) Q1afun(t,z), t , [0,0] , options);

 


function    rDot = Q1afun(t,z)
    g = getG ;
    m = getM ;
    k = getK ;
    
    OmegaH = sqrt((m*g)/(4*k)) ;

    z1 = z(1);
    z2 = z(2);
    
    z1dot = z2;
    if t < 1 
        z2dot = (4 * k)/m * (OmegaH + 70 * sin(0.5 * pi * t))^2;
    elseif t >= 1 
        z2dot = (4 * k)/m * (OmegaH - 77 * sin(0.5 * pi * t))^2;
    end

    rDot = [z1dot; z2dot];
end



% GLOBAL VARIABLES

function setGlobalM(val)
global m
m = val;
end

function r = getM
global m
r = m;
end

function setGlobalK(val)
global k
k = val;
end

function r = getK
global k
r = k;
end

function setGlobalG(val)
global g
g = val;
end

function r = getG
global g
r = g;
end

function setGlobalL(val)
global l
l = val;
end

function r = getL
global l
r = l;
end

function setGlobalB(val)
global b
b = val;
end

function r = getB
global b
r = b;
end