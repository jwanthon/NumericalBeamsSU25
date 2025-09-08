clc; clear all;
addpath('Functions\');
%% BISSBeamNumericFEM.m
%  Joseph Anthony
%
% Created:         9/8/25
% Last Modified:   9/8/25
%
% Description: Numerically solves the Euler-Bernoulli beam equation 
%   for a built in-simply supported beam using numeric FEM, determines the solution 
%   structure, and installs a time-domain solution.
%
%% Flags
FLAG.DISPLAY_MODES = true;
FLAG.DISPLAY_SPECTRUM = true;
FLAG.CALCULATE_TIMESOLVER = false;

%% Parameters
n = 100;         % Mesh size, number of interior points
shapes = 50;     % Number of shape functions

L = pi;           % Beam length
modeCount = 5;    % Number of displayed modes

beta = 10;        % Linear damping coefficient
IC_disp = .1;     % IC: Displacement [m]
IC_pos = .6;      % IC: Where displacement is located [m]
timespan = 1;     % Maximum time value, [s]
factor = 1;       % Scaling factor for K

%% Generate Modes
xvals = linspace(0,L,n+2);
deltax = mean(diff(xvals));
basis = zeros(shapes, n+2);
Dbasis = zeros(shapes, n+1);
D2basis = zeros(shapes, n);

% Generate shape functions and their respective derivatives
for i = 1:shapes
    basis(i,:) = xvals.^2.*cos(i*xvals).*(xvals-L*ones(1,n+2));
    Dbasis(i,:) = diff(basis(i,:))/deltax;
    D2basis(i,:) = diff(Dbasis(i,:))/deltax;
end

% Generate symmetric stiffness matrices
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

% Determine sorted eigenbasis
[evecs, evals] = eig(K);
[~, index] = sort(diag(abs(evals))); 
evals = evals(index,index);
evecs = evecs(:, index);

% Display mode shapes
modeshapes = basis'*evecs;
modeshapes = modeshapes';

if FLAG.DISPLAY_MODES
    figure();
    hold on
    for i = 1:modeCount
        plot(xvals, modeshapes(i,:));
    end
    xlabel("Beam Position [m]");
    ylabel("Relative Displacement");
    title("Mode Shapes of a Built In-Simply Supported Beam");
end


% Plot eigenvalue spectrum
espectrum = abs(diag(evals))/abs(evals(1,1));
scaled_espectrum = espectrum.^(1/4);

if FLAG.DISPLAY_SPECTRUM
    figure();
    hold on
    plot(scaled_espectrum, "LineStyle", "none", "Marker", ".");
    title("Adjusted BI-SS Eigenvalue Magnitudes");
    ylabel("Fourth Root of λ");
end

%% Time-Domain Solver
if FLAG.CALCULATE_TIMESOLVER
    K_scaled = -K * factor;
    IC_beam = xvals;
    for ind = 1:length(xvals)
        if xvals(ind) < IC_pos
            IC_beam(ind) = IC_disp/IC_pos*xvals(ind);
        else
            IC_beam(ind) = -IC_disp/(L-IC_pos)*(xvals(ind)-IC_pos)+IC_disp;
        end
    end
    
    % Find closest shape vector
    IC_shapes = [lsqr(basis', IC_beam')', zeros(1,shapes)];
    
    % Solve shape vector and translate to displacement
    [t,y] = ode45(@(t,y) odefcn_2orderFEM(t, y, K_scaled, M, beta), [0 timespan], IC_shapes);
    solution_beam = y(:,1:shapes)*basis;

    % Plot solution over time
    figure();
    surf(xvals,t,solution_beam, "LineStyle","none");
    camlight('headlight');
    title('Beam Displacement over Time');
    xlabel('Beam Position [m]');
    ylabel('Time [sec]');
    zlabel('Displacement [m]');
end