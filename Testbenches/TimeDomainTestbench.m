clc; clear all;
addpath('Functions\');
%% TimeDomainTestbench.m
%  Joseph Anthony
%
% Created:         9/17/25
% Last Modified:   9/17/25
%
% Description: Numerical solver analyzing an FEM solution for a BI-BI
%   beam, with the goal of optimizing its runtime
%
%% Parameters
n = 100;         % Mesh size, number of interior points
shapes = 50;     % Number of shape functions

L = pi;           % Beam length
modeCount = 5;    % Number of displayed modes

beta = 10;        % Linear damping coefficient
IC_disp = .1;     % IC: Displacement [m]
IC_pos = .6;      % IC: Where displacement is located [m]
timespan = 1;     % Maximum time value, [s]

%% Flags
DISPLAY_ICS = false;

%% Build K and M
xvals = linspace(0,L,n+2);
deltax = mean(diff(xvals));
basis = zeros(shapes, n+2);
Dbasis = zeros(shapes, n+1);
D2basis = zeros(shapes, n);

% Generate shape functions and their respective derivatives
for i = 1:shapes
    basis(i,:) = xvals.*(xvals-L*ones(1,n+2)).*sin(i*xvals);
    Dbasis(i,:) = diff(basis(i,:))/deltax;
    D2basis(i,:) = diff(Dbasis(i,:))/deltax;
end

% Generate K and M
K = zeros(shapes);
M = zeros(shapes);

for row = 1:shapes
    for col = 1:row
        K(row,col) = trapz(xvals(2:end-1), D2basis(row,:).*D2basis(col,:));
        M(row, col) = trapz(xvals, basis(row,:).*basis(col,:));
    end
end
K = K + K' - diag(diag(K));
M = M + M' - diag(diag(M));

%% Time-Domain Solver
K_scaled = -K;

% Generate initial conditions: overall triangle-shaped displacement
IC_beam = xvals;

for ind = 1:length(xvals)
    if xvals(ind) < IC_pos
        IC_beam(ind) = IC_disp/IC_pos*xvals(ind);
    else
        IC_beam(ind) = -IC_disp/(L-IC_pos)*(xvals(ind)-IC_pos)+IC_disp;
    end
end
if DISPLAY_ICS  
    plot(IC_beam); 
    title("BI-BI Beam Initial Conditions");
    xlabel("Position Index");
    ylabel("Displacement");
end

% Since M and K deal with shape vectors, we need to translate our ICs to a
%   shape vector, run it through the ODE solver, and then convert back to
%   displacement vectors
% Find closest shape vector
IC_shapes = [lsqr(basis', IC_beam')', zeros(1,shapes)];

% Solve shape vector and translate to displacement
tic
[t,y] = ode45(@(t,y) odefcn_2orderFEM(t, y, K_scaled, M, beta), [0 timespan], IC_shapes);
toc
solution_beam = y(:,1:shapes)*basis;