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

%{
% Question 2
init = [0,0,0,0,0,0,0,0,0,0,0,0] ;
t = linspace(0,6,1200) ;
[t,a] = ode45(@(t,d) Q2fun(t,d,I,[10,0,0,0,0,0,0,0]), t , init , options);
%}

%Q1 Plots
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

figure
plot(t,d(:,8))
title('Pvel(t) vs. t');
xlabel('t');
ylabel('Pvel');
ax = gca ;
exportgraphics(ax,'pvel.jpg')

figure
plot(t,d(:,10))
title('Qvel(t) vs. t');
xlabel('t');
ylabel('Qvel');
ax = gca ;
exportgraphics(ax,'qvel.jpg')

figure                      % questionable about how R vs t looks
plot(t,d(:,12))
title('Rvel(t) vs. t');
xlabel('t');
ylabel('Rvel');
ax = gca ;
exportgraphics(ax,'rvel.jpg')

% acceleration equation for Q1
function    dDot = Q1fun(t,d,I)
    % pulling constants into function
    g = getG ;
    m = getM ;
    k = getK ;
    l = getL ;
    b = getB ;
    
    % calculating hover omega
    OmegaH = sqrt((m*g)/(4*k)) ;

    % setting state space
    x1 = d(1) ;         % x position
    x2 = d(2) ;         % x velocity
    y1 = d(3) ;         % y position
    y2 = d(4) ;         % y velocity
    z1 = d(5) ;         % z position
    z2 = d(6) ;         % z velocity
    phi = d(7) ;        % phi
    p = d(8) ;          % p
    theta = d(9) ;      % theta
    q = d(10) ;         % q
    psi = d(11) ;       % psi
    r = d(12) ;         % r

    x1dot = x2 ;
    y1dot = y2 ;
    z1dot = z2 ;
   
    % part A function
    if t <= 2
        if t < 1
            Omega1 = OmegaH + 70 * sin(0.5 * pi * t) ;
            Omega2 = OmegaH + 70 * sin(0.5 * pi * t) ;
            Omega3 = OmegaH + 70 * sin(0.5 * pi * t) ;
            Omega4 = OmegaH + 70 * sin(0.5 * pi * t) ;
        elseif t >= 1 && t <= 2
            Omega1 = OmegaH - 77 * sin(0.5 * pi * t) ;
            Omega2 = OmegaH - 77 * sin(0.5 * pi * t) ;
            Omega3 = OmegaH - 77 * sin(0.5 * pi * t) ;
            Omega4 = OmegaH - 77 * sin(0.5 * pi * t) ;
        end
    
    % part B function
    elseif t > 2 && t < 4
        Omega1 = OmegaH ;
        Omega3 = OmegaH ;
        if t > 2 && t < 3
            Omega2 = sqrt((OmegaH^2 - 70^2 * sin(0.5*pi*(t-2)))) ;
            Omega4 = sqrt((OmegaH^2 + 70^2 * sin(0.5*pi*(t-2)))) ;
        elseif t >= 3 && t < 4
            Omega2 = sqrt((OmegaH^2 + 70^2 * sin(0.5*pi*(t-2)))) ;
            Omega4 = sqrt((OmegaH^2 - 70^2 * sin(0.5*pi*(t-2)))) ;
        end
        
    % part C function
    elseif t >= 4 && t <= 6
        Omega2 = OmegaH ;
        Omega4 = OmegaH ;
        if t >=4 && t < 5
            Omega1 = sqrt((OmegaH^2 - 70^2 * sin(0.5*pi*(t-4)))) ;
            Omega3 = sqrt((OmegaH^2 + 70^2 * sin(0.5*pi*(t-4)))) ;
        elseif t >=5 && t <= 6
            Omega1 = sqrt((OmegaH^2 + 70^2 * sin(0.5*pi*(t-4)))) ;
            Omega3 = sqrt((OmegaH^2 - 70^2 * sin(0.5*pi*(t-4)))) ;
        end
    end
    
    A = [ 1 , 0, -sin(theta) ;
          0 , cos(phi) , cos(theta)*sin(phi) ;
          0 , -sin(phi) , cos(theta)*cos(phi) ] ;
    
    inA = inv(A);

    % same fundamental equations use for parts A,B,C
    L = l*k*(-Omega2^2+Omega4^2) ;                       % Roll Moment
    M = l*k*(-Omega1^2+Omega3^2) ;                       % Pitch Moment
    N = b*(Omega1^2-Omega2^2+Omega3^2-Omega4^2) ;        % Yaw Moment

    pdot = ((q*r*(I(3,3)-I(2,2))) + L)/I(1,1) ;         % I was getting something different for pdot, qdot, and rdot, but I still need to check it out
    qdot = ((p*r*(I(1,1)-I(3,3))) + M)/I(2,2) ;
    rdot = ((p*q*(I(2,2)-I(1,1))) + N)/I(3,3) ;

    phidot = inA(1,:) * [p,q,r]' ;
    thetadot = inA(2,:) * [p,q,r]' ;
    psidot = inA(3,:) * [p,q,r]' ;

    x2dot = (k*(Omega1^2+Omega2^2+Omega3^2+Omega4^2)/m) * (cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi)) ;
    y2dot = (k*(Omega1^2+Omega2^2+Omega3^2+Omega4^2)/m) * (sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi)) ;
    z2dot =((k*(Omega1^2+Omega2^2+Omega3^2+Omega4^2)/m) * (cos(phi)*cos(theta))) - g ;

    % return
    dDot = [x1dot; x2dot; y1dot; y2dot; z1dot; z2dot; phidot; pdot; thetadot; qdot; psidot; rdot];
end

%{
function    dDot = Q2fun(t,d,I,rS)
    % pulling consts into function
    g = getG ;
    m = getM ;
    k = getK ;
    l = getL ;
    b = getB ;
    
    % calculating hover omega
    OmegaH = sqrt((m*g)/(4*k)) ;

    % setting state space
    x1 = d(1) ;             % x position
    x2 = d(2) ;             % x velocity
    y1 = d(3) ;             % y position
    y2 = d(4) ;             % y velocity
    z1 = d(5) ;             % z position
    z2 = d(6) ;             % z velocity
    phi = d(7) ;            % phi
    p = d(8) ;              % p
    theta = d(9) ;          % theta
    q = d(10) ;             % q
    psi = d(11) ;           % psi
    r = d(12) ;             % r

    x1dot = x2 ;
    y1dot = y2 ;
    z1dot = z2 ;
    
  
    

    % governing physics

    A = [ 1 , 0 , -sin(theta) ;
          0 , cos(phi) , cos(theta)*sin(phi) ;
          0 , -sin(phi) , cos(theta)*cos(phi) ] ;
    
    inA = inv(A);

    T = (g + (rS(2)-z2) + (rS(1)-z1)) * m/(cos(phi)*cos(theta)) ;
    L = I(1,1) * ((rS(3)-phidot) + (rS(4)-phi)) ; 
    M = I(2,2) * ((rS(5)-thetadot) + (rS(6)-theta)) ;
    N = I(3,3) * ((rS(7)-psidot) + (rS(8)-psi)) ;

    pdot = ((q*r*(I(3,3)-I(2,2))) + L)/I(1,1) ;
    qdot = ((p*r*(I(1,1)-I(3,3))) + M)/I(2,2) ;
    rdot = ((p*q*(I(2,2)-I(1,1))) + N)/I(3,3) ;

    phidot = inA(1,:) * [p,q,r]' ;
    thetadot = inA(2,:) * [p,q,r]' ;
    psidot = inA(3,:) * [p,q,r]' ;

    x2dot = (T/m) * (cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi)) ;
    y2dot = (T/m) * (sin(phi)*sin(theta)*cos(phi)-cos(psi)*sin(phi)) ;
    z2dot = (T/m) * (cos(phi)*cos(theta)) - g ;

    % return
    dDot = [x1dot; x2dot; y1dot; y2dot; z1dot; z2dot; phidot; pdot; thetadot; qdot; psidot; rdot];
end
%}

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