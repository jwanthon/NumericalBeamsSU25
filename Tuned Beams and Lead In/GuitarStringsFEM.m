clc; clear all;
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
n   = 10;                        % Parameter, mesh size

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
phi = zeros(n,n);          % ith row is ith test function, w/ coeffs x, x^2, ..., x^n
for i = 1:n
    for j = i:(2*i)
        phi(i,j) = (-L)^i * nchoosek(i,j-i);
    end
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

% Computing upper triangles of A and B
% Define Φi, Φj and their derivatives
syms phi_i(x) phi_j(x) Dphi_i(x) Dphi_j(x)

