clc; clear all;
addpath('Functions\');
%% SymbolicFEMTestbench.m
%  Joseph Anthony
%
% Created:         5/29/25
% Last Modified:   6/22/25
%
% Description: Testbench for different properties of and ways to generate 
%   FEM stiffness matrices using symbolic methods. This includes generating
%   symbolic basis functions and derivatives, symbolic integration, and
%   ways to optimize symbolic steps by taking pseudo-numerical routes.
%
% Required:
%   Symbolic Math Toolbox (Specified under section)           

n = 50; % Mesh size
L = 10;  % Length  (L) [m]
modeCount = 5;

%% Generating Basis Functions
% Create vector of test function polynomials and their derivative for the
%   guitar string problem
% tic/toc: ~ 0.85 sec (n = shapes = 50)
tic
syms phi_i(x)
phi = sym(1:n);
D2phi = sym(1:n);
for i = 1:n
    phi_i = x;
    for j = 1:i
        phi_i = times(phi_i, (x - L*j/i));      % Multiply in equally-spaced zeros
    end
    phi(i) = phi_i;
    D2phi(i) = diff(diff(phi_i, x),x);          % Find each function's 2nd derivative
end
toc

%% Generating Coefficient Matrix
% Find coefficient matrices for Dphi
% tic/toc: ~ 2.4 sec (n = shapes = 50)
tic
    coeff_D2phi = zeros(n,n+1);
for i = 1:n
    temp = coeffs(D2phi(i));
    coeff_D2phi(i,:) = [temp, zeros(1,n-length(temp)+1)];
end
toc
%% Determine Derivative Basis Effects
% tic/toc: ~166.19 sec (n = shapes = 50)
tic
% Turn derivative functions into usable form for MATLAB
D2Basis = zeros(n,n);
deltax = L/(n+1); 
for i = 1:n
    D2fn = matlabFunction(D2phi(i),x);
    for j = 1:n
        D2Basis(i,j) = D2fn(deltax*j);
    end
end
toc
%% Generate Stiffness Matrix - Poly Method
% Generate stiffness matrix by calculating coefficient vectors and 
%   doing dot products with polynomial integration coefficients
% tic/toc: ~ 0.161444 sec (n = shapes = 50)
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
        temp = multiplyCoeffs(coeff_D2phi(row,:), coeff_D2phi(col,:));
        K(row, col) = dot(temp, coeff_integration);
    end
    fprintf('Column %d of stiffness matrix calculated\n', col);
end
K = K + K' - diag(diag(K));
toc
%% Generate Stiffness Matrix - Trap Sum
% Generate stiffness matrix by integrating trap rules across each product
%   of derivative basis functions
%
% tic/toc: >0.01 sec (n = shapes = 50)
tic

xvals = linspace(deltax, n*deltax, n);

% Build stiffness matrix
K = zeros(n);

for col = 1:n
    for row = 1:col
        % Determine the current product function
        temp = DbeamBasis(row,:).*DbeamBasis(col,:);
        % Find the trap rule sum of this beam
        K(row, col) = trapz(xvals, temp);
    end
end
K = K + K' - diag(diag(K));

toc
%% Graph FEM Matrix
% Find and sort eigenvectors
[eigVecs, eigVals] = eig(K);
[d, index] = sort(diag(eigVals));
eigVals = eigVals(index,index);
eigVecs = eigVecs(:, index);
eigVals = flip(flip(eigVals, 1), 2);
eigVecs = flip(eigVecs, 1);

% Turn basis functions into usable form for MATLAB
ht = zeros(1,n);
for i = 1:n
    ht = matlabFunction(phi(i));
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
title(sprintf('First %d Mode Shapes of a BI-BI Beam', modeCount));
xlabel('Beam Position [m]');
ylabel('Displacement [m]');
hold off
%% Integration vs Trapz Values
syms testfn(x)
testfn = 1;
resolution = 10;
for i = 1:3
    testfn = plus(testfn, rand*x^i);
end
disp(testfn);
val1 = int(testfn, [0,10]);
ht = matlabFunction(testfn);
fnout = zeros(1,resolution);
for i = 1:resolution
    fnout(i) = ht(10/resolution*i);
end
val2 = trapz(fnout)/(resolution/10);
diff = abs(val1-val2);
fprintf("Val1 %d, Val2 %d, Diff %d", val1, val2, diff);
return
