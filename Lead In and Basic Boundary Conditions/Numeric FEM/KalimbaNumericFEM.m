clc; clear all;
addpath('Functions\');
%% KalimbaNumericFEM.m
%  Joseph Anthony
%
% Created:         7/16/25
% Last Modified:   7/16/25
%
% Description: Numerically solves the Euler-Bernoulli beam equation 
%   for a kalimba using numeric FEM and with piecewise functions.

n = 10000;      % Mesh size, number of interior points
shapes = 100;   % Number of shape functions used
L = 1;          % Beam length
a = 0.1;        % Position of the simple support
modeCount = 5;  % Number of displayed modes

%% Generate Piecewise Basis
% Basis function: x²(x-a)(1/h⁴)

xvals = linspace(0,L,n+2);
deltax = mean(diff(xvals));
h = L/(shapes+1);

basis = zeros(shapes,n+2);
Dbasis = zeros(shapes, n+1);
D2basis = zeros(shapes,n);

% for i = 1:shapes
%     x0 = (i-1)*L/(shapes+1);
%     for j = 1:length(xvals)
%         x = xvals(j);
%         if (i-1)*h-L/(n+1) < x
%             basis(i,j) = (x-x0)^2/h^4;
%         end
%     end
%     basis(i,:) = basis(i,:) .* (xvals - a*ones(1,n+2));
% end


% Generate derivatives
for i = 1:shapes
    Dbasis(i,:) = diff(basis(i,:))/deltax;
    D2basis(i,:) = diff(Dbasis(i,:))/deltax;
end  

%% Build Stiffness Matrix
K = zeros(shapes);
for row = 1:shapes
    for col = 1:row
        K(row,col) = trapz(xvals(2:end-1), D2basis(row,:).*D2basis(col,:));
        
    end
end

K = K + K' - diag(diag(K));

% Determine sorted eigenbasis
[evecs, evals] = eig(K);
[~, index] = sort(diag(abs(evals))); 
evals = evals(index,index);
evecs = evecs(:, index);

% Display mode shapes
modeshapes = basis'*evecs;
modeshapes = modeshapes';

for i = 1:modeCount
    modeshapes(i,:) = modeshapes(i,:)/max(abs(modeshapes(i,:)));
end

figure();
hold on
for i = 1:modeCount
    plot(xvals, modeshapes(i,:));
end
xlabel("Kalimba Position [m]");
ylabel("Relative Displacement");
title("Mode Shapes of a Kalimba");

% Plot eigenvalue spectrum
figure();
espectrum = abs(diag(evals))/abs(evals(1,1));

scaled_espectrum = espectrum.^(1/4);

hold on
plot(scaled_espectrum, "LineStyle", "none", "Marker", ".");
title("Adjusted Kalimba Eigenvalue Magnitudes");
ylabel("Fourth Root of λ");
