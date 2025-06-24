clc; clear all;
addpath('Functions\');
%% FixedFreeString.m
%  Joseph Anthony
%
% Created:         6/24/25
% Last Modified:   6/24/25
%
% Description: Numerically solves the Euler-Bernoulli beam equation 
%   for a string with a fixed end and with a free end, and compares its
%   eigenvalue structure with the eigenstructure of a fixed-fixed string.
%

n = 1000;       % Mesh size, number of interior points
shapes = 200;   % Number of shape functions, typically equal to n
L = pi;          % Beam length
modeCount = 5;  % Number of displayed modes

xvals = linspace(0,L,n+2);
fixed_basis = zeros(shapes, n+2); % Fixed-fixed string
fixed_Dbasis = zeros(shapes,n+1); 
free_basis = zeros(shapes, n+2); % Fixed-free string
free_Dbasis = zeros(shapes,n+1);
deltax = xvals(2)-xvals(1);

% Generate shape functions and their respective derivatives
for i = 1:shapes
    free_basis(i,:) = xvals.*cos(i*xvals);
    fixed_basis(i,:) = sin(i*xvals);
    free_Dbasis(i,:) = diff(free_basis(i,:))/deltax;
    fixed_Dbasis(i,:) = diff(fixed_basis(i,:))/deltax;
end

% Generate symmetric stiffness matrices
fixed_K = zeros(shapes);
free_K = zeros(shapes);

% Create a different mesh of the xvals
Dxvals = zeros(1, n+1);
for i = 1:length(xvals) - 1
    Dxvals(i) = (xvals(i) + xvals(i+1))/2;
end

for row = 1:shapes
    for col = 1:row
        fixed_K(row,col) = trapz(Dxvals, fixed_Dbasis(row,:).*fixed_Dbasis(col,:));
        free_K(row,col) = trapz(Dxvals, free_Dbasis(row,:).*free_Dbasis(col,:));
    end
end
fixed_K = fixed_K + fixed_K' - diag(diag(fixed_K));
free_K = free_K + free_K' - diag(diag(free_K));

% Generate nonsymmetric components
% Assumption: for a sufficiently large mesh and a smooth basis function,
%   the derivative at x = L is similar to the derivative at x = L - Δx
free_natK = zeros(shapes);
for row = 1:shapes
    for col = 1:shapes
        free_natK(row,col) = free_basis(col, end)*free_Dbasis(row, end);
    end
end

free_K = free_natK - free_K;

% Determine sorted eigenbasis
[free_evecs, free_evals] = eig(free_K);
[~, index] = sort(diag(abs(free_evals))); 
free_evals = free_evals(index,index);
free_evecs = free_evecs(:, index);

[fixed_evecs, fixed_evals] = eig(fixed_K);
[~, index] = sort(diag(abs(fixed_evals))); 
fixed_evals = fixed_evals(index,index);
fixed_evecs = fixed_evecs(:, index);

% Display mode shapes
fixed_modes = fixed_basis'*fixed_evecs;
fixed_modes = fixed_modes';
free_modes = free_basis'*free_evecs;
free_modes = free_modes';

figure();
fig1 = tiledlayout(1,2);
nexttile;
hold on
for i = 1:modeCount
    plot(xvals, fixed_modes(i,:));
end
xlabel("String Position [m]");
ylabel("Relative Displacement");
title("Mode Shapes of a Fixed-Fixed String");
nexttile;
hold on
for i = 1:modeCount
    plot(xvals, free_modes(i,:));
end
xlabel("String Position [m]");
ylabel("Relative Displacement");
title("Mode Shapes of a Fixed-Free String");

figure();

% Plot eigenvalue spectrum
fixed_espectrum = abs(diag(fixed_evals))/abs(fixed_evals(1,1));
free_espectrum = abs(diag(free_evals))/abs(free_evals(1,1));

fixed_rootvals = sqrt(fixed_espectrum);
free_rootvals = sqrt(sqrt(free_espectrum));

hold on
plot(fixed_rootvals, "LineStyle", "none", "Marker", ".");
plot(free_rootvals, "LineStyle", "none", "Marker", ".");
title("Adjusted Eigenvalue Magnitudes");
legend("Fixed-Fixed String (square root)", "Fixed-Free String (quarter root)");