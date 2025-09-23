clc; clear all;
addpath('Functions\');
%% DNStringNumericFEM.m
%  Joseph Anthony
%
% Created:         9/22/25
% Last Modified:   9/22/25
%
% Description: Numerically solves the wave equation for a string with 
%   Dirichlet-Neumann boundary conditions using FEM.
%   (x(0) = 0, x'(L) = 0).

n = 1000;       % Mesh size, number of interior points
shapes = 100;   % Number of shape functions, typically equal to n
L = pi;          % Beam length
modeCount = 5;  % Number of displayed modes

xvals = linspace(0,L,n+2);
basis = zeros(shapes, n+2); 
Dbasis = zeros(shapes,n+1); 

%% TODO: FIX DERIVATIVES OF SHAPE FUNCITONS AT x = L

deltax = xvals(2)-xvals(1);

% Generate shape functions and their respective derivatives
for i = 1:shapes
    basis(i,:) = sin((2*i-1)/2*xvals);
    Dbasis(i,:) = diff(basis(i,:))/deltax;
end

% Normalize basis functions
for i = 1:shapes
    basis(i,:) = basis(i,:) / (2i-1)*2;
end

% Generate symmetric stiffness matrices
K = zeros(shapes);

% Create a different mesh of the xvals
Dxvals = zeros(1, n+1);
for i = 1:length(xvals) - 1
    Dxvals(i) = (xvals(i) + xvals(i+1))/2;
end

for row = 1:shapes
    for col = 1:row
        K(row,col) = trapz(Dxvals, Dbasis(row,:).*Dbasis(col,:));
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

figure();
hold on
for i = 1:modeCount
    plot(xvals, modeshapes(i,:));
end
xlabel("String Position [m]");
ylabel("Relative Displacement");
title("Mode Shapes of a Dirichlet-Neumann String");

% Plot eigenvalue spectrum
espectrum = abs(diag(evals))/abs(evals(1,1));
rootvals = sqrt(espectrum);
figure();
plot(rootvals, "LineStyle", "none", "Marker", ".");
title("Squareroots of Eigenvalue Magnitudes");
