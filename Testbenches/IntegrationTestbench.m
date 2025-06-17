clc; clear all;
addpath('Functions\');
%% IntegrationTestbench.m
%  Joseph Anthony
%
% Created:         5/29/25
% Last Modified:   6/11/25
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
L = 10;  % Length  (L) [m]
modeCount = 5;

%% Generating Basis Functions
% Create vector of test function polynomials and their derivative
% tic/toc: ~ 0.85 sec
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
% tic/toc: ~ 2.4 sec
tic
    coeff_D2phi = zeros(n,n+1);
for i = 1:n
    temp = coeffs(D2phi(i));
    coeff_D2phi(i,:) = [temp, zeros(1,n-length(temp)+1)];
end
toc

%% Time multiplyCoeffs.m
A = rand(1,50)*10000;
B = rand(1,50)*10000;
f = @() multiplyCoeffs(A,B);
timeit(f)
%% Determine Derivative Basis Effects
% tic/toc: ~166.19 sec
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
xvals = linspace(deltax, n*deltax, n);

toc
%% Generate Stiffness Matrix - Poly Method
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
% tic/toc: >0.01 sec
tic

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
%% Numerical FEM Approach
tic
clc; clear all;
n = 100;         % Mesh size, number of interior points
shapes = 50;     % Number of shape functions, typically equal to n
L = 0.1;          % Beam length
modeCount = 3;  % Number of displayed modes

xvals       = linspace(0,L,n+2);
beamBasis   = zeros(shapes, n+2);
DbeamBasis  = zeros(shapes,n+2);
D2beamBasis = zeros(shapes,n+2);

% Generate shape functions
%   Current shape: equally spaced nodes x(x-L)(x-L/2) ...
for i = 1:shapes
    beamBasis(i,:) = xvals;
    for j = 1:i
        beamBasis(i,:) = beamBasis(i,:).*(xvals-j*L/i*ones(1,n+2));
    end
end


% % See shape functions
% hold on
% for i = 1:modeCount
%     plot(xvals,beamBasis(i,:));
% end
% return

% Generate second derivatives
for i = 1:shapes
    DbeamBasis(i,:) = gradient(beamBasis(i,:),xvals);
    D2beamBasis(i,:) = gradient(DbeamBasis(i,:),xvals);
end

% Build stiffness matrix
K = zeros(shapes);
for row = 1:shapes
    for col = 1:row
        product = D2beamBasis(row,:).*D2beamBasis(col,:);
        K(row,col) = trapz(product(2:n+1),xvals(2:n+1));
    end
end
K = K + K' - diag(diag(K));

% Find and sort eigenbasis from lowest mode to highest
[eigVecs, eigVals] = eig(K);
[d, index] = sort(diag(abs(eigVals))); % Sorts from highest mag to lowest
% [d, index] = sort(diag(eigVals));
eigVals = eigVals(index,index);
eigVecs = eigVecs(:, index);
eigVals = flip(flip(eigVals, 1), 2);% Used to sort from small to large
eigVecs = flip(eigVecs, 1);

% Create the mode shapes by scaling the basis functions by each eigenvector
modeShapes = beamBasis'*eigVecs;
modeShapes = modeShapes';
disp('Mode shapes calculated.');

% Normalize eigenfunctions
for i = 1:shapes
    modeShapes(i,:) = modeShapes(i,:)/max(abs(modeShapes(i,:)));
end

% Plot modecount eigenvectors
hold on;
for i = 1:modeCount
     plot(xvals,modeShapes(i,:),"Marker",".");
end
title(sprintf('First %d Mode Shapes of a Guitar String', modeCount));
xlabel('String Position [m]');
ylabel('Relative Displacement');
hold off

toc