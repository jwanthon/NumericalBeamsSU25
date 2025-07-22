clc; clear all;
addpath('Functions\');
%% FEMSolverTestbench.m
%  Joseph Anthony
%
% Created:          7/8/25
% Last Modified:    7/22/25
%
% Description: Testbench for solving the string and beam equation PDEs in
%   the time-domain using linear ODE solver tools and by implementing
%   linear damping to reduce ODE stiffness.
%
% Assumptions:
%  - Linear damping
%  - Zero initial velocity
%  - Linear interpolation for point displacement
%% Guitar String PDE Solver
% Parameters

n = 1000;         % Mesh size, number of interior points
shapes = 100;     % Number of shape functions

L = 1;            % String length
T = 1;            % String tension
rho = 1;          % String mass per unit length

beta = 10;        % Linear damping coefficient
IC_disp = .1;     % IC: Displacement [m]
IC_pos = .6;      % IC: Where displacement is located [m]
timespan = 10;     % Maximum time value, [s]

xvals       = linspace(0,L,n+2);
basis       = zeros(shapes, n+2);
Dbasis      = zeros(shapes,n+1);
D2basis     = zeros(shapes, n);
deltax      = L/(n+1);

% Create FEM Matrix
% Generate sinusoidal shape functions
for i = 1:shapes
    basis(i,:) = sin(pi*i*xvals/L);
    Dbasis(i,:) = diff(basis(i,:))/deltax;
end

% Build stiffness matrix
K = zeros(shapes);
M = zeros(shapes);
for row = 1:shapes
    for col = 1:row
        product = Dbasis(row,:).*Dbasis(col,:);
        K(row,col) = trapz(product);
        M(row, col) = trapz(basis(row,:).*basis(col,:));
    end
end
K = K + K' - diag(diag(K));
M = M + M' - diag(diag(M));

% Find and sort eigenbasis from lowest mode to highest
[eigVecs, eigVals] = eig(K);
[d, index] = sort(diag(abs(eigVals))); % Sort based on magnitude
eigVals = eigVals(index,index);
eigVecs = eigVecs(:, index);

% Create the mode shapes by scaling the basis functions by each eigenvector
modeShapes = basis'*eigVecs;
modeShapes = modeShapes';

% Solve time-domain solution
% Install ICs
K_scaled = -K * T / rho;
IC_string = xvals;
for ind = 1:length(xvals)
    if xvals(ind) < IC_pos
        IC_string(ind) = IC_disp/IC_pos*xvals(ind);
    else
        IC_string(ind) = -IC_disp/(L-IC_pos)*(xvals(ind)-IC_pos)+IC_disp;
    end
end

% Find closest shape vector
IC_shapes = [lsqr(basis', IC_string')', zeros(1,shapes)];

[t,y] = ode45(@(t,y) odefcn_2orderFEM(t, y, K_scaled, M, beta), [0 timespan], IC_shapes);


% Convert from shape vector to displacement vector
solution_string = y(:,1:shapes)*basis;

% Plot solution over time
surf(xvals,t,solution_string, "LineStyle","none");
camlight('headlight');
title('String Displacement over Time');
xlabel('String Position [m]');
ylabel('Time [sec]');
zlabel('Displacement [m]');

return
%% Troubleshooting w/ String FDM Matrix
clc; clear all;
addpath('Functions\');

n = 50;
rho = 1;                         
T   = 1;                     
L   = 1; 
beta = 0.2;
IC_pos = 0.1;
IC_disp = 0.1;
timespan = 10;

factor = 1000;

% Generate FDM matrix
deltax  = L/(n+1);
xvals = deltax:deltax:(L-deltax);
% temp1   = repmat(2, n, 1);
% temp2   = repmat(-1, n-1, 1);
% FDM       = diag(temp1) + diag(temp2, 1) + diag(temp2, -1);
% scaledFDM = -T / (rho * deltax^2) * FDM;


FDM = diag(repmat(-4, n-1, 1), 1) + diag(ones(n-2, 1), 2);
FDM = FDM + FDM';
FDM = FDM + diag(repmat(6, n, 1));
FDM(1,1) = 7;
FDM(n,n) = 7;
scaledFDM = -FDM * factor;

disp('FDM matrix generated.')

% Solve ODE K.u = K.u'' + βu'
IC_string = zeros(1,2*n);
for i = 1:n
    if xvals(i) < IC_pos
        IC_string(i) = IC_disp/IC_pos*xvals(i);
    else
        IC_string(i) = -IC_disp/(L-IC_pos)*(xvals(i)-IC_pos)+IC_disp;
    end
end

disp('ICs generated.')

[t,y] = ode45(@(t,y) odefcn_2orderFDM(t, y, scaledFDM, beta), [0 timespan], IC_string);

disp('System solved.')

surf(xvals,t,y(:,1:n), "LineStyle","none");
title('Sring Displacement over Time');
xlabel('String Position [m]');
ylabel('Time [sec]');
zlabel('Displacement [m]');