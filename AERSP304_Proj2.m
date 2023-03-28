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

t = linspace(0,2,2000) ;
z = zeros(2,1);
options = odeset('reltol',1e-12,'abstol',1e-12);
[t,z] = ode45(@(t,z) Q1afun(t,z), t , [0,0] , options);

 



function    rDot = Q1afun(t,r)
    r1 = r(1);
    r2 = r(2);
    
    r1dot = r2;
    if t < 1 
        r2dot = (4*k)/m*((sqrt((m*g)/(4*k)))+70*sin(0.5*pi*t))^2;
    elseif t >= 1 
        r2dot = (4*k)/m*((sqrt((m*g)/(4*k)))-77*sin(0.5*pi*t))^2;
    end

    rDot = [r1dot; r2dot];
end