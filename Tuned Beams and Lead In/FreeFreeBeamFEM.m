clc; clear all;
addpath('Functions\');
%% FreeFreeBeamFEM.m
%  Joseph Anthony
%
% Created:         6/3/25
% Last Modified:   6/3/25
%
% Description: Numerically solves the Euler-Bernoulli beam equation 
%   for a free-freebeam using FEM and determines the solution 
%   structure.
%
% Required:
%   Symbolic Math Toolbox               
%% Building Polynomial Test Functions
% Now, let's try to solve our model with FEM.
%   %   ρ * utt = Γ * uxxx                                      (eq. 1)
% with u''(0, t) = u''(L, t) = 0, u'''(0, t) = u'''(L, t) = 0.  (BCs)

rho         = 7.85;     % Mass per unit length [kg/m]
gamma       = 5.643;    % Flexural rigidity [N·m²]
L           = 1;        % Beam length [m]
n           = 50;       % Mesh size
modeCount   = 5;        % Number of displayed modes

% To satisfy the BCs, we will ensure that each function contains a product
%   of (x⁴+x+1) and (1/24*x⁴-L/6x³+L²/4*x²+x+1) so that the second and third derivatives at
%   the endpoints are zero, but that the displacement and first derivatives
%   are nonzero (no loading conditions).

% Create vector of test function polynomials and their derivatives
syms phi_i(x)
phi = sym(1:n);
Dphi = sym(1:n);
for i = 1:n
    fprintf('Creating test function %d\n', i);
    phi_i = (x^4+x+1)*(1/24*x^4-L/6*x^3+L^2/4*x^2+x+1);
    for j = 1:i
        phi_i = times(phi_i, (x - L*j/i));      % Multiply in equally-spaced zeros
    end
    phi(i) = phi_i;
    Dphi(i) = diff(phi_i, x);                   % Find each function's derivative
end

%% Building the Mass and Stiffness Matrices
% We want (i, j) of K to have value ∫ Φi' * Φj', and (i, j) of M to have
%   value ∫ Φi * Φj.

% Computing upper triangles of M and K (since they're symmetric)
M = zeros(n, n);
K = M;

for col = 1:n
    for row = 1:col
        M(row, col) = double(int(phi(row)  * phi(col),  [0, L]));
        K(row, col) = double(int(Dphi(row) * Dphi(col), [0, L]));
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

% Determine the effect of each basis function on the string
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
     plot(linspace(0,L,n+2),[0, modeShapes(:,i)', 0],"Marker",".")
end
title(sprintf('First %d Mode Shapes of a Free-Free Beam', modeCount));
xlabel('Beam Position [m]');
ylabel('Displacement [m]');
hold off

% Cleanup
clear row col i j 

% Create a visual of the FEM matrix
% fig2 = visual_sparseMatrix(FEM);
% fig2.Title = "Free-Free FEM Matrix";