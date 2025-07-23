clc; clear all;
addpath('Functions\');
%% KalimbaNumericFEM.m
%  Joseph Anthony
%
% Created:         7/16/25
% Last Modified:   7/23/25
%
% Description: Numerically solves the Euler-Bernoulli beam equation 
%   for a kalimba using numeric FEM and with piecewise functions, 
%   determines the solution structure, and installs a time-domain solution.
%
%% Flags
FLAG.DISPLAY_MODES = true;
FLAG.DISPLAY_SPECTRUM = true;
FLAG.CALCULATE_TIMESOLVER = true;
FLAG.CALCULATE_ASPECTRA = false;

%% Parameters
n = 100;          % Mesh size, number of interior points
shapes = 50;      % Number of shape functions

L = 1;           % Beam length
a = L/5;          % Position of the simple support
modeCount = 5;    % Number of displayed modes

spectra = 100;    % Number of a values to calculate for spectra

beta = 0;        % Linear damping coefficient
IC_disp = .1;     % IC: Displacement [m]
timespan = 1;     % Maximum time value, [s]
factor = 100;     % Scaling factor for K

%% Generate Modes
xvals = linspace(0,L,n+2);
deltax = mean(diff(xvals));
h = L/(shapes+1);

basis = zeros(shapes,n+2);
Dbasis = zeros(shapes, n+1);
D2basis = zeros(shapes,n);

% Basis = sin(i²x)(x-a)
for i = 1:shapes
    basis(i,:) = sin(i^2*xvals).*(xvals-a*ones(1,n+2));
end


% Generate derivatives
for i = 1:shapes
    Dbasis(i,:) = diff(basis(i,:))/deltax;
    D2basis(i,:) = diff(Dbasis(i,:))/deltax;
end  

% Build matrices
K = zeros(shapes);
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

for i = 1:modeCount
    modeshapes(i,:) = modeshapes(i,:)/max(abs(modeshapes(i,:)));
end

if FLAG.DISPLAY_MODES
    figure();
    hold on
    for i = 1:modeCount
        plot(xvals, modeshapes(i,:));
    end
    xlabel("Kalimba Position [m]");
    ylabel("Relative Displacement");
    title("Mode Shapes of a Kalimba");
end

% Plot eigenvalue spectrum
espectrum = abs(diag(evals))/abs(evals(1,1));
scaled_espectrum = espectrum.^(1/4);

if FLAG.DISPLAY_SPECTRUM
    figure();
    hold on
    plot(scaled_espectrum, "LineStyle", "none", "Marker", ".");
    title("Adjusted Kalimba Eigenvalue Magnitudes");
    ylabel("Fourth Root of λ");
end

%% Time-Domain Solver
% IC Based on cantilever beam
% (https://mechanicalc.com/reference/beam-deflection-tables)
if FLAG.CALCULATE_TIMESOLVER
    K_scaled = -K * factor;

    % Install ICs
    lever_length = (L-a);
    lever_constant = IC_disp/(2*lever_length^3);
    
    IC_beam = xvals;
    for ind = 1:length(xvals)
        if xvals(ind) < a
            IC_beam(ind) = 0;
        else
            IC_beam(ind) = lever_constant*(xvals(ind)-a)^2*(3*lever_length-(xvals(ind)-a));
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

%% Recompute FEM for different values a
if FLAG.CALCULATE_ASPECTRA
    a_vec = linspace(L/(spectra+1), L-L/(spectra+1), spectra); % All positions of the simple support
    for iteration = 1:length(a_vec)
        a = a_vec(iteration);

        % Basis = sin(i²x)(x-a)
        for i = 1:shapes
            basis(i,:) = sin(i^2.*xvals).*(xvals-a*ones(1,n+2));
            Dbasis(i,:) = diff(basis(i,:))/deltax;
            D2basis(i,:) = diff(Dbasis(i,:))/deltax;
        end  

        % Build K
        K = zeros(shapes);
        for row = 1:shapes
            for col = 1:row
                K(row,col) = trapz(xvals(2:end-1), D2basis(row,:).*D2basis(col,:));
            end
        end
        K = K + K' - diag(diag(K));
        
        % Determine sorted eigenbasis
        [evecs, evals] = eig(K);
        [~, index] = sort(diag(abs(evals))); 
        evals = evals(index,index);
        evecs = evecs(:, index);
        
        % Find eigenvalue spectrum
        espectrum = abs(diag(evals))/abs(evals(1,1));
        scaled_espectrum = espectrum.^(1/4);
    
        full_aspectra(iteration,:) = scaled_espectrum;
    end
    figure();
    surf(full_aspectra);
end

