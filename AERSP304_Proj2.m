% QUAD-ROTOR AERSP 304
% By : Nicholas Giampetro, Craig Stenstrom, Payton Glynn

clc;
clear;
close all;

% setting constants values
setGlobalG(9.81) ; % m/s
setGlobalM(0.450) ; % kg
setGlobalL(0.225) ; % 
setGlobalK(2.98*10^-7) ;
setGlobalB(1.14*10^-6) ;

I = [ 4.85*10^-3 , 0, 0 ;
     0 , 4.85*10^-3 , 0 ;
     0 , 0 , 8.80*10^-3 ] ;

% setting ODE45 options
options = odeset('reltol',1e-12,'abstol',1e-12);

% Question 1 Part A
t = linspace(0,2,2000) ;
[t,r] = ode45(@(t,r) Q1fun(t,r,I), t , [0,0,0,0,0,0,0,0,0,0,0,0] , options);

figure
plot(t,r(:,5))
title('z(t) vs. t');
xlabel('t');
ylabel('z');
ax = gca ;
exportgraphics(ax,'zQ1A.jpg')

figure
plot(t,r(:,6))
title('Vz(t) vs. t');
xlabel('t');
ylabel('Vz');
ax = gca ;
exportgraphics(ax,'vzQ1A.jpg')



% acceleration equation for Q1
function    rDot = Q1fun(t,d,I)
    % pulling consts into function
    g = getG ;
    m = getM ;
    k = getK ;
    
    % calculating hover omega
    OmegaH = sqrt((m*g)/(4*k)) ;

    % setting state space
    x1 = d(1);
    x2 = d(2);
    y1 = d(3);
    y2 = d(4);
    z1 = d(5);
    z2 = d(6);
    phi1dot = d(7);
    phi2dot = d(8);
    theta1dot = d(9);
    theta2dot = d(10);
    psi1dot = d(11);
    psi2dot = d(12);

    % setting default accelerations
    x1dot = 0 ;
    x2dot = 0 ;
    y1dot = 0 ;
    y2dot = 0 ;
    z1dot = 0 ;
    z2dot = 0 ;
    phi1dot = 0 ;
    phi2dot = 0 ;
    theta1dot = 0 ;
    theta2dot = 0 ;
    psi1dot = 0 ;
    psi2dot = 0 ;

    
    % part A function
    z1dot = z2;
    if t < 1 
        T1 = k*(OmegaH + 70 * sin(0.5 * pi * t))^2 ;
        T2 = k*(OmegaH + 70 * sin(0.5 * pi * t))^2 ;
        T3 = k*(OmegaH + 70 * sin(0.5 * pi * t))^2 ;
        T4 = k*(OmegaH + 70 * sin(0.5 * pi * t))^2 ;
        z2dot = (T1+T2+T3+T4)/m - g ;
    elseif t >= 1 && t <= 2
        T1 = k*(OmegaH - 77 * sin(0.5 * pi * t))^2 ;
        T2 = k*(OmegaH - 77 * sin(0.5 * pi * t))^2 ;
        T3 = k*(OmegaH - 77 * sin(0.5 * pi * t))^2 ;
        T4 = k*(OmegaH - 77 * sin(0.5 * pi * t))^2 ;
        z2dot = (T1+T2+T3+T4)/m - g ;
    end
    
    % part B function

    % return
    rDot = [x1dot; x2dot; y1dot; y2dot; z1dot; z2dot; phi1dot; phi2dot; theta1dot; theta2dot; psi1dot; psi2dot];
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