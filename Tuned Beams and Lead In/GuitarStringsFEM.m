clc; clear all;
addpath('Functions\');
%% GuitarStringsFEM.m
%  Joseph Anthony
%
% Created:         5/9/25
% Last Modified:   5/29/25
%
% Description: Numerically solves the wave equation for a guitar string of
%  an abitrary length using FEM and determines the solution structure.
%
% Required:
%   Symbolic Math Toolbox                                   (Add-on)
%% Building Polynomial Test Functions
% Now, let's try to solve our model with FEM.
%   utt = T / ρ * uxx                                          (eq. 1)
% with u(x, t) = u(L, t) = 0.                                  (BC)

rho = 1;                          % Parameter, density (ρ) [kg/m^2]
T   = 1;                          % Parameter, tension (T) [N]
L   = 1;                          % Parameter, length  (L) [m]
n   = 50;                         % Parameter, mesh size
modeCount = 5;                    % Parameter, # of displayed modes

% This time, we will approximate our solution in the form of the following:
%   u(x, t) = Σ αi(t) * Φi(x).                                (eq. 2)
%  
% We want αi to describe our solution's time-domain response, and Φi, or 
%   our test function, to describe our soltuions' spatial response. Each
%   test function needs to be indepenent.
%
% One easy way to do this is to build polynomials with zeroes equally
%   spaced zeroes between x = 0 and x = L, with the ith test function Φ
%   having i-1 zeroes in between the boundaries.

% Create vector of test function polynomials and their derivatives
syms phi_i(x)
phi = sym(1:n);
Dphi = sym(1:n);
for i = 1:n
    fprintf('Creating test function %d\n', i);
    phi_i = x;
    for j = 1:i
        phi_i = times(phi_i, (x - L*j/i));      % Multiply in equally-spaced zeros
    end
    phi(i) = phi_i;
    Dphi(i) = diff(phi_i, x);                   % Find each function's derivative
end

%% Procedure for the FDM Matrix
% Let T / ρ = k. Then utt = k * uxx.
%   Now, for any "test function" v = Φi, which acts as a solution component 
%   in eq. 3, our equation will be unaltered. That is,
%   v * utt = k * v * uxx                                       (eq. 3)
%   
%   ∫ v*utt dx = k * ∫ v*uxx dx                   (integrating across L)
%   ∫ Φi * utt dx = k * ∫ Φi * uxx dx                     (v is some Φi)
%   ∫ Φi * Σ (αj'' * Φj) dx = k∫ Φi * Σ (αj * Φj'') dx           (eq. 2)
%   Σ aj'' * (∫ Φi * Φj)dx = k*Σ(αj * ∫(Φi * Φj'')dx)
%   Σ aj'' * (∫ Φi * Φj)dx = k*Σ(αj * ∫(Φi' * Φj')dx) (int. by parts)
%
% Notice that the integrals on the LHS and RHS are both able to easily
%   determined by integrating across L, since all Φi and Φj are 
%   polynomials. We can describe this through the matrix-vector equation:
%   M.a'' = K.a                                   (M from LHS, K from RHS)
%   
%   where K represents the stiffness matrix, and M represents the mass
%   matrix.
%
%% Building the Mass and Stiffness Matrices
% We want (i, j) of K to have value ∫ Φi' * Φj', and (i, j) of M to have
%   value ∫ Φi * Φj. Note that both matrices are symmetric!

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
%% 

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
     plot(modeShapes(:,i),"Marker",".")
end
title(sprintf('First %d Mode Shapes of a Guitar String', modeCount));
xlabel('String Position (i)');
ylabel('Displacement [m]');

% Calculate angular frequencies
angFreq = diag(sqrt(T*abs(eigVals)/(rho*deltax^2)));
% Add respective angular frequencies
legends = zeros(1:modeCount);
for i = 1:modeCount
    legends(i) = string(sprintf('%.2f [rad/s]', angFreq(i))); 
end
legend(legends);
hold off;

% Cleanup
clear row col i j 

% Create a visual of the FEM matrix
fig2 = visual_sparseMatrix(FEM);
fig2.Title = "Guitar String FEM Matrix";