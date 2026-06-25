clc; clear all; close all;
l = 0:.01:(10/100 - .01); % different lengths from the base(m)
n = 10;    % number of nodes

T = zeros(n,1); %initial guess values

Tb = 1000; % base temperature
Ttip = -273.16;    % tip temperature


% solution

k = 200; %number of iterations

%case 1
for j=1:1:k

  T(1,1) = Tb; %first node

  for i = 2:1:n-1
    T(i,1) = (T(i+1,1) + T(i-1,1)) / 2; % temp inside the domain
  end

  T(n,1) = Ttip; %last node


  figure(1);
  plot(l, T); hold on;
  title("Temperature distribution between 1000 degree C base and absolute zero temperature tip of a straight rod")
  xlabel("Length from the base (m)");
  ylabel("Temperature (Degree C)");
  grid on;

end

hold off;

%case 2

for j=1:1:k

  T(1,1) = T(2,1); %first node

  for i = 2:1:n-1
    T(i,1) = (T(i+1,1) + T(i-1,1)) / 2; % temp inside the domain
  end

  T(n,1) = Ttip; %last node

  figure(2)
  plot(l, T); hold on;
  title("Temperature distribution between adiabatic base(insulated) and absolute zero temperature tip of a straight rod")
  xlabel("Length from the base (m)");
  ylabel("Temperature (Degree C)");
  grid on;

end

hold off;






