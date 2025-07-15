clc; clear all;
addpath('Functions\');
%% PiecewiseBeamNumericFEM.m
%  Joseph Anthony
%
% Created:         7/14/25
% Last Modified:   7/15/25
%
% Description: Numerically solves the Euler-Bernoulli beam equation 
%   for a built in-built in beam using numeric FEM and with different types
%   of piecewise functions.

n = 10000;       % Mesh size, number of interior points
shapes = 100;     % Number of shape functions used
L = 1;         % Beam length
modeCount = 5;  % Number of displayed modes

%% Generate Different Bases
xvals = linspace(0,L,n+2);
deltax = mean(diff(xvals));
h = L/(shapes+1);

% Type 1, Hump: x²(x-2h)²(1/h⁴)
hump_basis = zeros(shapes,n+2);
hump_Dbasis = zeros(shapes, n+1);
hump_D2basis = zeros(shapes,n);

for i = 1:shapes
    x0 = (i-1)*L/(shapes+1);
    for j = 1:length(xvals)
        x = xvals(j);
        if (i-1)*h < x && x < (i+1)*h
            hump_basis(i,j) = (x-x0)^2*(x-x0-2*h)^2/h^4;
        end
    end
end

% Type 2, Wiggle:  x²(x-2h)²(x/h⁴-1/h³)
wigg_basis = zeros(shapes,n+2);
wigg_Dbasis = zeros(shapes, n+1);
wigg_D2basis = zeros(shapes,n);

for i = 1:shapes
    x0 = (i-1)*L/(shapes+1);
    for j = 1:length(xvals)
        x = xvals(j);
        if (i-1)*h < x && x < (i+1)*h
            wigg_basis(i,j) = (x-x0)^2*(x-x0-2*h)^2*(1/h^4*(x-x0)-1/h^3);
        end
    end
end
hold on
% Generate derivatives
for i = 1:shapes
    hump_Dbasis(i,:) = diff(hump_basis(i,:))/deltax;
    hump_D2basis(i,:) = diff(hump_Dbasis(i,:))/deltax;
    wigg_Dbasis(i,:) = diff(wigg_basis(i,:))/deltax;
    wigg_D2basis(i,:) = diff(wigg_Dbasis(i,:))/deltax;
end  

beambasis = [hump_basis ; wigg_basis];

%% Build Stiffness Matrix

% Hump submatrix
hump_K = zeros(shapes);
for row = 1:shapes
    for col = 1:row
        hump_K(row,col) = trapz(xvals(2:end-1), hump_D2basis(row,:).*hump_D2basis(col,:));
        
    end
end

% Wiggle submatrix
wigg_K = zeros(shapes);
for row = 1:shapes
    for col = 1:row
        wigg_K(row,col) = trapz(xvals(2:end-1), hump_D2basis(row,:).*hump_D2basis(col,:));
        
    end
end

% % Build the full matrix
% K = [hump_K,        zeros(shapes);
%      zeros(shapes), wigg_K];

K = wigg_K;

K = K + K' - diag(diag(K));

% Determine sorted eigenbasis
[evecs, evals] = eig(K);
[~, index] = sort(diag(abs(evals))); 
evals = evals(index,index);
evecs = evecs(:, index);

% % Display mode shapes
% modeshapes = beambasis'*evecs;
% modeshapes = modeshapes';

modeshapes = wigg_basis'*evecs
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
