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
M_ode = [M, zeros(shapes);
         zeros(shapes), zeros(shapes)];

F = ode;
F.InitialValue = IC_shapes;
F.ODEFcn = @(t, y) [ y(shapes+1:end);
                     K_scaled*y(1:shapes)];
F.MassMatrix = odeMassMatrix(MassMatrix = M, Singular = "yes");

solve(F, 0, timespan);


% Convert from shape vector to displacement vector
solution_string = solution_shape(:,1:shapes)*beamBasis;

% Plot solution over time
surf(xvals,t,solution_string, "LineStyle","none");
camlight('headlight');
title('String Displacement over Time');
xlabel('String Position [m]');
ylabel('Time [sec]');
zlabel('Displacement [m]');


%% Troubleshooting w/ String FDM Matrix
clc; clear all;
addpath('Functions\');

n = 100;
rho = 1;                         
T   = 1;                     
L   = 1; 
beta = 10000;
IC_pos = 0.1;
IC_disp = 0.1;
timespan = 100;

factor = 1;

% Generate FDM matrix
deltax  = L/(n+1);
xvals = deltax:deltax:(L-deltax);
% temp1   = repmat(2, n, 1);
% temp2   = repmat(-1, n-1, 1);
% FDM       = diag(temp1) + diag(temp2, 1) + diag(temp2, -1);
% scaledFDM = -T / (rho * deltax^2) * FDM;

%%
%%
%% CHANGE IC TO THE EXPECTED ONE FOR THE BI-BI AND RETEST
%%
%%

FDM = diag(repmat(-4, n-1, 1), 1) + diag(ones(n-2, 1), 2);
FDM = FDM + FDM';
FDM = FDM + diag(repmat(6, n, 1));
FDM(1,1) = 7;
FDM(n,n) = 7;
scaledFDM = FDM / factor;

% Solve ODE K.u = K.u'' + βu'
IC_string = zeros(1,2*n);
for i = 1:n
    if xvals(i) < IC_pos
        IC_string(i) = IC_disp/IC_pos*xvals(i);
    else
        IC_string(i) = -IC_disp/(L-IC_pos)*(xvals(i)-IC_pos)+IC_disp;
    end
end

[t,y] = ode45(@(t,y) odefcn_2orderFDM(t, y, scaledFDM, beta), [0 timespan], IC_string);

surf(xvals,t,y(:,1:100), "LineStyle","none");
title('String Displacement over Time');
xlabel('String Position [m]');
ylabel('Time [sec]');
zlabel('Displacement [m]');