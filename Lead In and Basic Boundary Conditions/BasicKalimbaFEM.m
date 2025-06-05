clc; clear all;
addpath('Functions\');
%% BIBIBeamFEM.m
%  Joseph Anthony
%
% Created:         5/20/25
% Last Modified:   5/27/25
%
% Description: Solves the Euler-Bernoulli beam equation for a simple model
%   of a kalimba using FEM techniques.
%
% Required:
%   Symbolic Math Toolbox               
%% BCs Building Polynomial Test Functions 
% ρ * utt = Γ * uxxx                                           (E-B)
%
% Clamped at x = 0,     u(0)   = u'(0)   = 0
% free-free at x = L,   u''(L) = u'''(L) = 0
% no motion at x = b,   u(b)   = 0                             (BCs)

rho         = 7.85;       % Mass per unit length [kg/m]
gamma       = 5.643;        % Flexural rigidity [N·m²]
L           = 1;        % Beam length [m]
b           = 0.25;      % Bridge point [m]
n           = 20;       % Mesh size
modeCount   = 5;        % Number of displayed modes

% Initial function: 
syms phi0(x);
phi0 = x^2*(x-b)*((x-L)^4-x);

% Create vector of test function polynomials and their derivatives
syms phi_i(x)
phi = sym(1:n);
D2phi = sym(1:n);
for i = 1:n
    fprintf('Creating test function %d\n', i);
    phi_i = phi0;
    for j = 1:i
        if j ~= 1
            phi_i = times(phi_i, (x - L*j/i));      % Multiply in equally-spaced zeros
        end
    end
    phi(i) = phi_i;
    D2phi(i) = diff(diff(phi_i, x),x);                   % Find each function's derivative
end

%% Building the Mass and Stiffness Matrices
% Kij = ∫ Φi * Φj, Mij = ∫ Φi' * Φj'.

% Computing upper triangles of M and K
M = zeros(n, n);
K = M;

for col = 1:n
    for row = 1:col
        M(row, col) = double(int(phi(row)  * phi(col),  [0, L]));
        K(row, col) = double(int(D2phi(row) * D2phi(col), [0, L]));
    end
    fprintf('Column %d of M and K matrices calculated\n', col);
end

% Reflect upper triangles across the main diagonal and correct duplicate
%   entries on the main diagonal
M = M + M' - diag(diag(M));
K = K + K' - diag(diag(K));

%% Working with FEM Matrix
FEM = K;
disp('FEM matrix built.');

% Find and sort eigenvectors
[eigVecs, eigVals] = eig(FEM);
[d, index] = sort(diag(eigVals));
eigVals = eigVals(index,index);
eigVecs = eigVecs(:, index);
eigVals = flip(flip(eigVals, 1), 2);
eigVecs = flip(eigVecs, 1);

% Turn basis functions into usable form for MATLAB
ht = zeros(1,n);
for i = 1:n
    ht = matlabFunction(phi);
end
deltax = L/(n+1);

% Determine the effect of each basis function on the beam
beamBasis = zeros(n,n);
for i = 1:n
    beamBasis(i,:) = ht(deltax*i);
end
disp('Basis function string effects calculated.');

% Normalize the basis functions
for i = 1:n
    beamBasis(:,i) = beamBasis(:,i)/max(abs(beamBasis(:,i)));
end

% Create the mode shapes by scaling the basis functions by each eigenvector
modeShapes = beamBasis*eigVecs;
disp('Mode shapes calculated.');

% Plot modecount eigenvectors
hold on;
for i = 1:modeCount
     plot(linspace(0,L,n+1),[0, modeShapes(:,i)'],"Marker",".")
end
title(sprintf('First %d Mode Shapes of a Basic Kalimba', modeCount));
xlabel('Beam Position [m]');
ylabel('Displacement [m]');
hold off

% Cleanup
clear row col i j 

% % Create a visual of the FEM matrix
% fig2 = visual_sparseMatrix(FEM);
% fig2.Title = "Basic Kalimba FEM Matrix";