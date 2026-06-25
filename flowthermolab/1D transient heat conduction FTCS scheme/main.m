% 1d Heat heat conduction problem using FTCS scheme for transient problem
clc; clear all; close all;


% Geometry

xmin = 0;
xmax = 0.2; % rod length(meter)
N = 20     %number of nodes
dx = (xmax - xmin) / (N - 1) %grid space
x = xmin:dx:xmax;
dt = 1e-3
tmax = 0.25; % total simulation time
t = 0:dt:tmax;
alpha = 0.05; %diffusion coefficient

% Initial and Boundary Conditions

Tcurrent = ones(1,N)* 30; %initial temperature across the rod
Tbase = 300;
Ttip = 50;

% Solution

d = (alpha * dt) / (dx^2)  %for parabolic equation as such

for j=2:length(t) % loop for time step

   T = Tcurrent;

   for i=1:N      % Space stepping

     if i == 1
       T(i) = Tbase;
     elseif i == N
       T(i) = Ttip;
     else
       T(i) = T(i) + (d * (T(i + 1) - 2 * T(i) + T(i - 1)));
     end

   endfor

   Tcurrent = T;
   time = j*dt

   plot(x, Tcurrent)
   xlabel("Length of the rod(m)")
   ylabel("Temperature (^{o}C)")
   shg
   pause(0.01)

end
