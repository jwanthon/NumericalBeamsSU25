clc; clear all; hold off;
%% GuitarStringsFEM.m
%  Joseph Anthony
%
% Created:         5/9/25
% Last Modified:   5/9/25
%
% Description: Numerically solves the wave equation for a guitar string of
%  an abitrary length using FEM and determines the solution structure.
%
% Required:
%   Symbolic Math Toolbox                                   (Add-on)
%% TODO:
% - Build matrices A and B
% - Build FDM matrix
% - Analyze FDM matrix modes
% - Explain integration by parts more
%% Building Polynomial Test Functions
% Now, let's try to solve our model with FEM.
%   utt = T / ρ * uxx                                          (eq. 1)
% with u(x, t) = u(L, t) = 0.                                  (BC)

rho = 1;                          % Parameter, density (ρ) [kg/m^2]
T   = 1;                          % Parameter, tension (T) [N]
L   = 1;                          % Parameter, length  (L) [m]
n   = 5;                          % Parameter, mesh size
modeCount = 5;                    % Parameter, # of displayed modes


% This time, we will approximate our solution in the form of the following:
%   u(x, t) = Σ αi(t) * Φi(x).                                (eq. 2)
%  
% We want αi to describe our solution's time-domain response, and Φi, or 
%   our test function, to describe our soltuions' spatial response. Each
%   test function needs to be indepenent.
%
% One easy way to do this is to build polynomials with zeroes at x = 0 and
%   x = L. Let Φ1 = x(x-L), Φ2 = x(x-L)^2, and so on.

% Build coefficient matrix for test functions Φ
disp('Building coefficient matrix...');
phi = zeros(n,n);          % ith row is ith test function, w/ coeffs x, x^2, ..., x^n
for i = 1:n
    for j = i:(2*i)
        phi(i,j) = (-L)^i * nchoosek(i,j-i);
    end
end
disp('Coefficient matrix complete.');
%% Procedure for the FDM Matrix
% Let T / ρ = k. Then utt = k * uxx.
%   Now, for any "test function" v = Φi, which acts as a solution component 
%   in eq. 3, our equation will be unaltered. That is,
%   v * utt = k * v * uxx                                       (eq. 3)
%   
%   ∫ v*utt dx = k * ∫ v*uxx dx                   (integrating across L)
%   ∫ Φi * utt dx = k * ∫ Φi * uxx dx                     (v is some Φi)
%   ∫ Φi * Σ (αj'' * Φj) dx = k∫ Φi * Σ (αj * Φj'') dx           (eq. 2)
%   Σ aj'' * (∫ Φi'' * Φj'')dx = k*Σ(αj * ∫(Φi * Φj'')dx)
%   Σ aj'' * (∫ Φi'' * Φj'')dx = k*Σ(αj * ∫(Φi' * Φj')dx) (int. by parts)
%
% Notice that the integrals on the LHS and RHS are both able to easily
%   determined by integrating across L, since all Φi and Φj are 
%   polynomials. We can describe this through the matrix-vector equation:
%   B.a'' = A.a                                   (A from LHS, B from RHS)
%   a'' = K.a                                     (K = inv(A)B)  (eq. 3)
%
%% Building the FDM Matrix
% We want (i, j) of B to have value ∫ Φi * Φj, and (i, j) of A to have
%   value ∫ Φi' * Φj'. Note that A and B are both symmetrix.

% Computing lower triangles of A and B
% Define Φi, Φj and their derivatives

A = zeros(n, n);
B = A;

syms phi_i(x) phi_j(x) Dphi_i(x) Dphi_j(x) intA(x) intB(x)
for i = 1:n
    phi_i = 0;
    phi_j = 0;
    for degree = i:(2*i)
        phi_i = plus(phi_i, phi(i, degree) * x^degree);
    end
    Dphi_i = diff(phi_i, x);
    for j = 1:i
        for degree = j:(2*j)
            phi_j = plus(phi_j, phi(j, degree) * x^degree);
        end
        Dphi_j = diff(phi_j, x);
        A(i,j) = double(int(Dphi_i * Dphi_j, [0, L]));
        B(i,j) = double(int(phi_i * phi_j, [0, L]));
    end
    disp(sprintf('Line %d of triangular matrices complete.', i));
end

% Manipulate and build into FEM matrix
A = A + A' - diag(A);
B = B + B' - diag(B);
FEM = B * inv(A);
disp('FEM matrix built.');

[eigVecs, eigVals] = eig(FEM);

hold on
for i = 1:modeCount
    plot(eigVecs(:,i),"Marker",".", "LineStyle","none")
end
title(sprintf('First %d Mode Shapes of a Guitar String', modeCount));
xlabel('String Position (i)');
ylabel('Displacement [m]');
hold off