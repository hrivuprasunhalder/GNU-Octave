% Convection–Diffusion in 1D using Finite Volume Method
% Central differencing scheme
% User specifies number of nodes N

clear; clc;

% --- Input parameters ---
L = 1.0;          % domain length [m]
rho = 1.0;        % density [kg/m^3]
Gamma = 0.1;      % diffusion coefficient [kg/m·s]
u = input('Enter velocity u [m/s]: ');
N = input('Enter number of nodes: ');

% --- Derived quantities ---
dx = L / (N - 1); % grid size
F = rho * u;      % advection strength
D = Gamma / dx;   % diffusion strength

% --- Coefficient arrays ---
aW = zeros(N,1);
aE = zeros(N,1);
aP = zeros(N,1);
Su = zeros(N,1);
SP = zeros(N,1);

% --- Boundary conditions ---
phi0 = 1;         % at x = 0
phiL = 0;         % at x = L

% --- Coefficient for interior nodes ---
for i = 2:N-1
   aW(i) = D + F/2;
   aE(i) = D - F/2;
   SP(i) = 0;
   Su(i) = 0;
   aP(i) = aW(i) + aE(i) - SP(i);
endfor

% --- Boundary nodes ---
% Node 1 (x =0)
aW(1) = 0;
aE(1) = D - F/2;
SP(1) = -2*(D + F);
Su(1) = (2*D + F)*phi0;
aP(1) = aW(1) + aE(1) - SP(1);

% Node N (x = L)
aW(N) = D + F/2;
aE(N) = 0;
SP(N) = -2*(D - F);
Su(N) = (2*D - F)*phiL;
aP(N) = aW(N) + aE(N) - SP(N);


% --- Assemble coefficient matrix and Right Hand Side constant matrix ---
A = zeros(N,N);
b = zeros(N,1);

for i = 1:N
  A(i,i) = aP(i);
  if i > 1
    A(i, i-1) = -aW(i);
  endif
  if i < N
    A(i, i+1) = -aE(i);
  endif
  b(i) = Su(i);
endfor

% --- Solve System ---
phi = A\b;
disp(phi);

