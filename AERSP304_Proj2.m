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
t = linspace(0,4,1000) ;
[t,d] = ode45(@(t,d) Q1fun(t,d,I), t , [0,0,0,0,0,0,0,0,0,0,0,0] , options);

figure
plot(t,d(:,7))
title('phi(t) vs. t');
xlabel('t');
ylabel('phi');
ax = gca ;
exportgraphics(ax,'phi.jpg')

figure
plot(t,d(:,3))
title('y(t) vs. t');
xlabel('t');
ylabel('y');
ax = gca ;
exportgraphics(ax,'y.jpg')

figure
plot(t,d(:,4))
title('Vy(t) vs. t');
xlabel('t');
ylabel('Vy');
ax = gca ;
exportgraphics(ax,'vy.jpg')

figure
plot(t,d(:,5))
title('z(t) vs. t');
xlabel('t');
ylabel('z');
ax = gca ;
exportgraphics(ax,'z.jpg')

figure
plot(t,d(:,6))
title('Vz(t) vs. t');
xlabel('t');
ylabel('Vz');
ax = gca ;
exportgraphics(ax,'vz.jpg')



% acceleration equation for Q1
function    rDot = Q1fun(t,d,I)
    % pulling consts into function
    g = getG ;
    m = getM ;
    k = getK ;
    l = getL ;
    
    % calculating hover omega
    OmegaH = sqrt((m*g)/(4*k)) ;

    % setting state space
    x1 = d(1);
    x2 = d(2);
    y1 = d(3);
    y2 = d(4);
    z1 = d(5);
    z2 = d(6);
    phi1 = d(7);
    phi2 = d(8);
    theta1 = d(9);
    theta2 = d(10);
    psi1 = d(11);
    psi2 = d(12);

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
    T1 = k*OmegaH^2 ;
    T3 = k*OmegaH^2 ;

    if t > 2 && t < 3
        T2 = k*(OmegaH^2 - 70^2 * sin(0.5*pi*(t-2))) ;
        T4 = k*(OmegaH^2 + 70^2 * sin(0.5*pi*(t-2))) ;
        L = l*(-T2+T4) ;
        phi1dot = phi2 ;
        phi2dot = L/I(1,1) ;
        z1dot = z2 ;
        z2dot = ((T1+T2+T3+T4)/m) * cos(phi1) - g ;
        y1dot = y2 ;
        y2dot = ((T1+T2+T3+T4)/m) * sin(phi1) ;
    elseif t >=3 && t < 4
        T2 = k*(OmegaH^2 + 70^2 * sin(0.5*pi*(t-2))) ;
        T4 = k*(OmegaH^2 - 70^2 * sin(0.5*pi*(t-2))) ;
        L = l*(-T2+T4) ;
        phi1dot = phi2 ;
        phi2dot = L/I(1,1) ;
        z1dot = z2 ;
        z2dot = ((T1+T2+T3+T4)/m) * cos(phi1) - g ;
        y1dot = y2 ;
        y2dot = ((T1+T2+T3+T4)/m) * sin(phi1) ;
    end

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