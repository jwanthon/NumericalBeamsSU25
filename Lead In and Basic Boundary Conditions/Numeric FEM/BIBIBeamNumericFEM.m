clc; clear all;
addpath('Functions\');
%% BIBIBeamNumericFEM.m
%  Joseph Anthony
%
% Created:         6/24/25
% Last Modified:   6/24/25
%
% Description: Numerically solves the Euler-Bernoulli beam equation 
%   for a built in-built in beam using numeric FEM and determines the solution 
%   structure.
%

n = 1000;       % Mesh size, number of interior points
shapes = 500;   % Number of shape functions
L = pi;         % Beam length
modeCount = 5;  % Number of displayed modes

xvals = linspace(0,L,n+2);
deltax = mean(diff(xvals));
beambasis = zeros(shapes, n+2);
Dbeambasis = zeros(shapes, n+1);
D2beambasis = zeros(shapes, n);

% Generate shape functions and their respective derivatives
for i = 1:shapes
    beambasis(i,:) = xvals.*(xvals-L*ones(1,n+2)).*sin(i*xvals);
    Dbeambasis(i,:) = diff(beambasis(i,:))/deltax;
    D2beambasis(i,:) = diff(Dbeambasis(i,:))/deltax;
end

% Generate symmetric stiffness matrices
K = zeros(shapes);

for row = 1:shapes
    for col = 1:row
        K(row,col) = trapz(xvals(2:end-1), D2beambasis(row,:).*D2beambasis(col,:));
    end
end
K = K + K' - diag(diag(K));

% Determine sorted eigenbasis
[evecs, evals] = eig(K);
[~, index] = sort(diag(abs(evals))); 
evals = evals(index,index);
evecs = evecs(:, index);

% Display mode shapes
modeshapes = beambasis'*evecs;
modeshapes = modeshapes';

figure();
hold on
for i = 1:modeCount
    plot(xvals, modeshapes(i,:));
end
xlabel("Beam Position [m]");
ylabel("Relative Displacement");
title("Mode Shapes of a Built In-Built In Beam");


% Plot eigenvalue spectrum
figure();
espectrum = abs(diag(evals))/abs(evals(1,1));

scaled_espectrum = espectrum.^(1/4);

hold on
plot(scaled_espectrum, "LineStyle", "none", "Marker", ".");
title("Adjusted BI-BI Eigenvalue Magnitudes");
ylabel("Fourth Root of λ");
