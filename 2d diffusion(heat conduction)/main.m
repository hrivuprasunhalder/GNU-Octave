clear all;
close all;
clc;

% Geometry
L = 10e-2; %10 cm
H = 10e-2; %10 cm

% Meshing
nx = 10;
ny = 10;   % number of nodes

x = linspace(0, L, nx);
y = linspace(0, H, ny);
dx = L / (nx-1);
dy = H / (ny-1);

% Initialization
T = zeros(nx, ny);

% Boundary conditions
T(:, 1) = 80;
T(1, :) = 100;
T(end, :) = 50;
T(:, end) = 20;
T_old = T;

beta = (dx/dy)^2;
err = 200;
tol = 1e-3;
k = 0;

while err > tol
  k = k + 1;
  for i=2:nx-1
    for j=2:ny-1
      T(i,j) = (1 / (2 * (1 + (beta)^2))) * (T(i+1, j) + T(i-1, j) + ((beta^2) * T(i, j+1)) + ((beta^2) * T(i, j-1)));
    endfor
  endfor

  err = abs(max(max(T - T_old)));
  err_p(k) = err;
  T_old = T;

  shg;
  pause(0.001);

  figure(1);
  plot(err_p);
  xlabel('Number of ierations');
  ylabel('Error/Residual');
  title('Residual Graph');

  figure(2);
  contourf(x, y, T, 'ShowText', 'on'); %In MATLAB and GNU Octave, contour draws contour lines (isolines) only, while contourf fills the regions between those lines with colors based on the colormap.
  xlabel('Length');
  ylabel('Height');
  title(' 2D steady state heat conduction');

end






