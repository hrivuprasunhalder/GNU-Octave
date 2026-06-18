% Gauss Elimination Method
% Co-eff matrix A
% constants b
% x_gauss = solution vector

function x_gauss = GauEli(A,b)

% Augmented matrix [A|b]
Ab = [A b];

%number of equations
n = size(A, 1);

for k = 1:n-1
    %pivoting
    [~, idx] = max(abs(Ab(k:n,k)));
    idx = idx + k - 1;
    if idx ~= k
        Ab([k idx],:) = Ab([idx k],:);
    end
    %elimination
    for i = k+1:n
        factor = Ab(i,k) / Ab(k,k);
        Ab(i,:) = Ab(i, :) - factor * Ab(k, :);
    end
end

% Back substitution
x_gauss = zeros(n,1);
for i = n:-1:1
    x_gauss(i) = (Ab(i,end) - Ab(i,1:n)*x_gauss) / Ab(i,i);
end


