clc; clear all;
addpath('Functions\');
%% FFFFBeamSymbolicFEM.m
%  Joseph Anthony
%
% Created:         6/3/25
% Last Modified:   7/7/25
%
% Description: Numerically solves the Euler-Bernoulli beam equation 
%   for a free-free free-free beam using FEM and determines the solution 
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

% Note that both of our boundary conditions are natural boundary
%   conditions, so they affect the "leftover terms" when converting E-B to
%   the weak form.

deltax = L/(n+1);

% Create vector of test function polynomials and their derivatives
syms phi_i(x)
phi = sym(1:n);
phi_i = 0;
D2phi = sym(1:n);
for i = 1:n
    fprintf('Creating test function %d\n', i);
    phi_i = add(phi_i, (x-deltax)^2)
    phi(i) = phi_i;
    D2phi(i) = diff(diff(phi_i, x), x);     % Find each function's 2nd derivative
end

%% Building the Mass and Stiffness Matrices
% We want (i, j) of K to have value ∫ Φi' * Φj', and (i, j) of M to have
%   value ∫ Φi * Φj.

% Computing upper triangles of M and K (since they're symmetric)
K = zeros(n, n);

for col = 1:n
    for row = 1:col
        K(row, col) = double(int(D2phi(row) * D2phi(col), [0, L]));
    end
    fprintf('Column %d of K calculated. \n', col);
end

% Reflect upper triangles across the main diagonal and correct duplicate
%   entries on the main diagonal
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
     plot(linspace(deltax,L-deltax,n),modeShapes(:,i)',"Marker",".")
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

return
%% Alternate Attempt
clc; clear all;
addpath('Functions\');

rho         = 7.85;     % Mass per unit length [kg/m]
gamma       = 5.643;    % Flexural rigidity [N·m²]
L           = 1;        % Beam length [m]
n           = 50;       % Mesh size
modeCount   = 5;        % Number of displayed modes
distance = L/10;               % Set total length of test parabola

deltax = L/(n+1);
a = -4/distance^2;

% Create vector of test functions
syms phi_i(x)
phi = sym(1:n);
phi_i = 0;
D2phi = sym(1:n);
for i = 1:n
    fprintf('Creating test function %d\n', i);
    phi_i = piecewise((i*deltax - distance/2 < x) & (x < i*deltax + distance/2), a*(x-i*deltax)^2+1,0);
    phi(i) = phi_i;
    D2phi(i) = piecewise((i*deltax - distance/2 < x) & (x < i*deltax + distance/2), 2*a,0);
end

K = zeros(n, n);

for col = 1:n
    for row = 1:col
        K(row, col) = double(int(D2phi(row) * D2phi(col), [0, L]));
    end
    fprintf('Column %d of K calculated. \n', col);
end

% Reflect upper triangles across the main diagonal and correct duplicate
%   entries on the main diagonal
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

step = distance/2/deltax;
step = floor(step);

% Determine the effect of each basis function on the beam
beamBasis = zeros(n,n);
for i = 1:n
    j = max(1, i-step);
    while j < min(n, i+step)
        beamBasis(i,j) = (a*(deltax^2)*(j-i)^2+1)
        j=j+1;
    end
end
disp('Basis function string effects calculated.');
beamBasis(:,n) = zeros(n,1);

% Normalize the basis functions
% for i = 1:n
%     beamBasis(:,i) = beamBasis(:,i)/max(abs(beamBasis(:,i)));
% end

% Create the mode shapes by scaling the basis functions by each eigenvector
modeShapes = beamBasis*eigVecs;
disp('Mode shapes calculated.');

% Plot modecount eigenvectors
hold on;
for i = 1:modeCount
     plot(linspace(deltax,L-deltax,n),modeShapes(:,i)',"Marker",".")
end
title(sprintf('First %d Mode Shapes of a Free-Free Beam', modeCount));
xlabel('Beam Position [m]');
ylabel('Displacement [m]');
hold off

% Cleanup
clear row col i j 