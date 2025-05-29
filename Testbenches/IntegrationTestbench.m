clc; clear all;
addpath('Functions\');
%% IntegrationTestbench.m
%  Joseph Anthony
%
% Created:         5/29/25
% Last Modified:   5/29/25
%
% Description: Testbench for different integration techniques for 
%   calculating FEM stiffness matrices with large quantities of basis 
%   functions. This uses the basis functions defined in GuitarStringsFEM.m.
%   
%   All tic/toc tests use n = 50 and L = 1 and approximate the timing for 
%   rough estimates.
%
% Required:
%   Symbolic Math Toolbox           

n = 50; % Mesh size
L = 1;  % Length  (L) [m]

%% Generating Basis Functions
% Create vector of test function polynomials and their derivative
% tic/toc: ~ 0.85 sec
tic
syms phi_i(x)
phi = sym(1:n);
Dphi = sym(1:n);
for i = 1:n
    phi_i = x;
    for j = 1:i
        phi_i = times(phi_i, (x - L*j/i));      % Multiply in equally-spaced zeros
    end
    phi(i) = phi_i;
    Dphi(i) = diff(phi_i, x);                   % Find each function's derivative
end
toc

%% Generating Coefficient Matrix
% Find coefficient matrices for Dphi
% tic/toc: ~ 2.4 sec
tic
    coeff_Dphi = zeros(n,n+1);
for i = 1:n
    temp = coeffs(Dphi(i));
    coeff_Dphi(i,:) = [temp, zeros(1,n-length(temp)+1)];
end
toc

%% Time multiplyCoeffs.m
A = rand(1,50)*10000;
B = rand(1,50)*10000;
f = @() multiplyCoeffs(A,B);
timeit(f)

%% Generate Stiffness Matrix
% Generate stiffness matrix by calculating coefficient vectors and 
%   doing dot products with polynomial integration coefficients
% tic/toc: ~ 0.161444 sec
tic

% Generate polynomial integration coefficients
coeff_integration = zeros(1,n+1);
for i = 1:2*n+2
    coeff_integration(i) = L^i/i;
end

% Build stiffness matrix
K = zeros(n);

for col = 1:n
    for row = 1:col
        % Figure out the element's new coefficient vector
        temp = multiplyCoeffs(coeff_Dphi(row,:), coeff_Dphi(col,:));
        K(row, col) = dot(temp, coeff_integration);
    end
    fprintf('Column %d of stiffness matrix calculated\n', col);
end
K = K + K' - diag(diag(K));
toc

