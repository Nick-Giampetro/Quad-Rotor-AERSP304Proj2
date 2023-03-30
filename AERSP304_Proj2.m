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

% Question 1
init = [0,0,0,0,0,0,0,0,0,0,0,0] ;
t = linspace(0,6,1200) ;
[t,d] = ode45(@(t,d) Q1fun(t,d,I), t , init , options);

%Q1 Plots
figure
plot(d(:,1),d(:,3))
title('x(t) vs. y(t)');
xlabel('x');
ylabel('y');
ax = gca ;
exportgraphics(ax,'XvY.jpg')

figure
plot(t,d(:,7))
title('phi(t) vs. t');
xlabel('t');
ylabel('phi');
ax = gca ;
exportgraphics(ax,'phi.jpg')

figure
plot(t,d(:,9))
title('theta(t) vs. t');
xlabel('t');
ylabel('theta');
ax = gca ;
exportgraphics(ax,'theta.jpg')

figure
plot(t,d(:,11))
title('psi(t) vs. t');
xlabel('t');
ylabel('psi');
ax = gca ;
exportgraphics(ax,'psi.jpg')

figure
plot(t,d(:,1))
title('x(t) vs. t');
xlabel('t');
ylabel('x');
ax = gca ;
exportgraphics(ax,'x.jpg')
figure
plot(t,d(:,2))
title('Vx(t) vs. t');
xlabel('t');
ylabel('Vx');
ax = gca ;
exportgraphics(ax,'vx.jpg')

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
    b = getB ;
    
    % calculating hover omega
    OmegaH = sqrt((m*g)/(4*k)) ;

    % setting state space
    x1 = d(1) ;
    x2 = d(2) ;
    y1 = d(3) ;
    y2 = d(4) ;
    z1 = d(5) ;
    z2 = d(6) ;
    phi1 = d(7) ;
    phi2 = d(8) ;
    theta1 = d(9) ;
    theta2 = d(10) ;
    psi1 = d(11) ;
    psi2 = d(12) ;

    phi1dot = phi2 ;
    theta1dot = theta2 ;
    psi1dot = psi2 ;
    x1dot = x2 ;
    y1dot = y2 ;
    z1dot = z2 ;
   
    % part A function
    if t <= 2
        if t < 1
            T1 = k*(OmegaH + 70 * sin(0.5 * pi * t))^2 ;
            T2 = k*(OmegaH + 70 * sin(0.5 * pi * t))^2 ;
            T3 = k*(OmegaH + 70 * sin(0.5 * pi * t))^2 ;
            T4 = k*(OmegaH + 70 * sin(0.5 * pi * t))^2 ;
        elseif t >= 1 && t <= 2
            T1 = k*(OmegaH - 77 * sin(0.5 * pi * t))^2 ;
            T2 = k*(OmegaH - 77 * sin(0.5 * pi * t))^2 ;
            T3 = k*(OmegaH - 77 * sin(0.5 * pi * t))^2 ;
            T4 = k*(OmegaH - 77 * sin(0.5 * pi * t))^2 ;
        end
    
    % part B function
    elseif t > 2 && t < 4
        T1 = k*OmegaH^2 ;
        T3 = k*OmegaH^2 ;
        if t > 2 && t < 3
            T2 = k*(OmegaH^2 - 70^2 * sin(0.5*pi*(t-2))) ;
            T4 = k*(OmegaH^2 + 70^2 * sin(0.5*pi*(t-2))) ;
        elseif t >= 3 && t < 4
            T2 = k*(OmegaH^2 + 70^2 * sin(0.5*pi*(t-2))) ;
            T4 = k*(OmegaH^2 - 70^2 * sin(0.5*pi*(t-2))) ;
        end
        
    % part C function
    elseif t >= 4 && t <= 6
        T2 = k*OmegaH^2 ;
        T4 = k*OmegaH^2 ;
        if t >=4 && t < 5
            T1 = k*(OmegaH^2 - 70^2 * sin(0.5*pi*(t-4))) ;
            T3 = k*(OmegaH^2 + 70^2 * sin(0.5*pi*(t-4))) ;
        elseif t >=5 && t <= 6
            T1 = k*(OmegaH^2 + 70^2 * sin(0.5*pi*(t-4))) ;
            T3 = k*(OmegaH^2 - 70^2 * sin(0.5*pi*(t-4))) ;
        end
    end
    
    % same fundamental equations use for parts A,B,C
    L = l*(-T2+T4) ; 
    M = l*(-T1+T3) ;
    N = (b/k)*(T1-T2+T3-T4) ;
    phi2dot = L/I(1,1) ;
    theta2dot = M/I(2,2) ;
    psi2dot = N/I(3,3) ;
    x2dot = ((T1+T2+T3+T4)/m) * (cos(psi1)*sin(theta1)*cos(phi1)+sin(psi1)*sin(phi1)) ;
    y2dot = ((T1+T2+T3+T4)/m) * (sin(phi1)*sin(theta1)*cos(phi1)-cos(psi1)*sin(phi1)) ;
    z2dot = ((T1+T2+T3+T4)/m) * (cos(phi1)*cos(theta1)) - g ;

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