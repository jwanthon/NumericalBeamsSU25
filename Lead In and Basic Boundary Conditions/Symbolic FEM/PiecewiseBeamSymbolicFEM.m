clc; clear all;
addpath('Functions\');
%% PiecewiseBeamNumericFEM.m
%  Joseph Anthony
%
% Created:         7/18/25
% Last Modified:   7/18/25
%
% Description: Analytically solves the Euler-Bernoulli beam equation 
%   for a built in-built in beam using symbolic FEM and with different types
%   of piecewise functions.
% 
% Required:
%   Symbolic Math Toolbox      

shapes = 20;   % Number of shape functions used
n = shapes;
L = 1;         % Beam length
modeCount = 5; % Number of displayed modes

%% Generate Different Bases
syms phi_i(x)
phi = sym(1:shapes*2);
D2phi = sym(1:shapes*2);

h = L/(shapes+1);

for i = 1:shapes
    phi(i) = piecewise(((i-1)*h < x) & ((i+1)*h > x), (x-(i-1)*L/(n+1))^2*((x-(i-1)*L/(n+1))-2*h)^2/h^4,0);
    D2phi(i) = diff(diff(phi(i), x),x);   
end
fprintf('Humps generated.\n')

for i = 1:shapes
    phi(shapes+i) = piecewise(((i-1)*h < x) & ((i+1)*h > x),(x-(i-1)*L/(n+1))^2*((x-(i-1)*L/(n+1))-2*h)^2*(-1/h^3+(x-(i-1)*L/(n+1))/h^4),0);
    D2phi(i) = diff(diff(phi(i), x),x);   
end
fprintf('Wiggles generated.\n');

% Computing K
K = zeros(shapes*2);

for col = 1:(shapes*2)
    for row = 1:col
        K(row, col) = double(int(D2phi(row) * D2phi(col), [0, L]));
    end
end

K = K + K' - diag(diag(K));

SYMBOLIC_K = K
%% Working with FEM Matrix
disp('Stiffness matrix built.');

% Create a visual of the FEM matrix
fig2 = visual_sparseMatrix(K);
fig2.Title = "Stiffness Matrix for a BI-BI Beam";

