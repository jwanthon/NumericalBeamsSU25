clc; clear all;
addpath('Functions\');
%% FEMSolverTestbench.m
%  Joseph Anthony
%
% Created:          7/8/25
% Last Modified:    7/8/25
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

beta = 1000;        % Linear damping coefficient
IC_disp = .1;     % IC: Displacement [m]
IC_pos = .6;      % IC: Where displacement is located [m]
timespan = 1;     % Maximum time value, [s]

xvals       = linspace(0,L,n+2);
beamBasis   = zeros(shapes, n+2);
DbeamBasis  = zeros(shapes,n+1);
deltax      = L/(n+1);

% Create FEM Matrix
% Generate sinusoidal shape functions
for i = 1:shapes
    beamBasis(i,:) = sin(pi*i*xvals/L);
    DbeamBasis(i,:) = diff(beamBasis(i,:))/deltax;
end

% Build stiffness matrix
K = zeros(shapes);
M = zeros(shapes);
for row = 1:shapes
    for col = 1:row
        product = DbeamBasis(row,:).*DbeamBasis(col,:);
        K(row,col) = trapz(product);
        M(row, col) = trapz(beamBasis(row,:).*beamBasis(col,:));
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
modeShapes = beamBasis'*eigVecs;
modeShapes = modeShapes';

% Solve time-domain solution
% Install ICs
K_scaled = K * T / rho;
IC_string = xvals;
for ind = 1:length(xvals)
    if xvals(ind) < IC_pos
        IC_string(ind) = IC_disp/IC_pos*xvals(ind);
    else
        IC_string(ind) = -IC_disp/(L-IC_pos)*(xvals(ind)-IC_pos)+IC_disp;
    end
end
% Find closest shape vector
IC_shapes = [lsqr(beamBasis', IC_string')', zeros(1,shapes)];

[t,y] = ode89(@(t,y) odefcn_2orderFEM(t, y, K_scaled, M, beta), [0 timespan], IC_shapes);

% Convert from shape vector to displacement vector
solution = y(:,1:shapes)*beamBasis;

% Plot solution over time
surf(xvals,t,solution, "LineStyle","none");
camlight('headlight');
title('String Displacement over Time');
xlabel('String Position [m]');
ylabel('Time [sec]');
zlabel('Displacement [m]');
