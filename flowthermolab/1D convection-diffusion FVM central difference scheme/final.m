% Convection–Diffusion in 1D using Finite Volume Method
% Central differencing scheme
% User specifies number of nodes N

clc;
clear;
close all;

% --- Provided data ---
L = 1.0;          % domain length [m]
rho = 1.0;        % density [kg/m^3]
Gamma = 0.1;      % diffusion coefficient [kg/m·s]

% --- Input parameters ---
u = input('Enter velocity u [m/s]: ');
N = input('Enter number of control volumes: ');

% --- Derived quantities ---
dx = L / N; % grid size
F = rho * u;      % advection strength
D = Gamma / dx;   % diffusion strength

% --- Boundary conditions ---
phi0 = 1;         % at x = 0, left boundary
phiL = 0;         % at x = L, right boundary

% ---Matrix Form---
A = zeros(N,N);
b = zeros(N,1);

for i=1:N
  if i==1
    % left boundary
    aW = 0;
    aE = D - F/2;
    Sp = -(2*D + F);
    Su = (2*D + F)*phi0;
  elseif i==N
    % right boundary
    aW = D + F/2;
    aE = 0;
    Sp = -(2*D - F);
    Su = (2*D - F)*phiL;
  else
    % interior cells
    aW = D + F/2;
    aE = D - F/2;
    Sp = 0;
    Su = 0;
  end

  aP = aW + aE - Sp;

  %fill co-eff matrix
  A(i,i) = aP;

  if i>1
    A(i,i-1) = -aW;
  end

  if i<N
    A(i,i+1) = -aE;
  end

  b(i) = Su;

end

% --- Solve ---
phi = A\b;
disp(phi);
Pe = F / D;
disp("Peclet Number: ");
disp(Pe);
