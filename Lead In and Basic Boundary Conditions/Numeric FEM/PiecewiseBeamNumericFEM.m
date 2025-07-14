clc; clear all;
addpath('Functions\');
%% PiecewiseBeamNumericFEM.m
%  Joseph Anthony
%
% Created:         7/14/25
% Last Modified:   7/14/25
%
% Description: Numerically solves the Euler-Bernoulli beam equation 
%   for a built in-built in beam using numeric FEM and with different types
%   of piecewise functions.

n = 1000;       % Mesh size, number of interior points
types = 2;      % Number of types of shape functions; normally set to 2
L = pi;         % Beam length
modeCount = 5;  % Number of displayed modes
piecemesh = 100; % Piecewise mesh size, number of points spanned

%% Generate # of types of piecewise functions
xvals = linspace(0,L,n+2);
deltax = mean(diff(xvals));
h = piecemesh/2

% Type 1, Hump: x²(x-2h)²(1/h⁴)
hump = zeros(1, piecemesh);
for i = 1:piecemesh
    hump(i) = xvals(i)^2*(xvals(i)-xvals(2*h))^2/(xvals(h)^4);
end

% Type 2, Wiggle:  x²(x-2h)²(x/h⁴-1/h³)
wiggle = zeros(1, piecemesh);
for i = 1:piecemesh
    wiggle(i) = xvals(i)^2*(xvals(i)-xvals(2*h))^2*(xvals(i)/xvals(h)^4-1/xvals(h)^3);
end

%% Populate basis functions (TODO: MAKE VAR NAMING MORE DESCRIPTIVE)
shapes = types*(n+2-piecemesh);
beambasis = zeros(shapes,n+2);

for i = 1:shapes/types+1
    beambasis(i, i:i+piecemesh-1) = hump;
end

for i = 1:shapes/types+1
    beambasis(i+shapes/types, i:i+piecemesh-1) = wiggle;
end

Dbeambasis = zeros(shapes,n+1);
D2beambasis = zeros(shapes, n);
for i = 1:shapes
    Dbeambasis(i,:) = diff(beambasis(i,:))/deltax;
    D2beambasis(i,:) = diff(Dbeambasis(i,:))/deltax;
end


% Generate symmetric stiffness matrices
K = zeros(length(beambasis));

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

for i = 1:modeCount
    modeshapes(i,:) = modeshapes(i,:)/max(abs(modeshapes(i,:)));
end

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
