clc; clear all;
addpath('Functions\');

%% BeamExamplesFEM.m
%  Joseph Anthony
%
% Created:         11/08/25
% Last Modified:   11/08/25
%
% Description: Produces several diagrams of BI-BI, BI-SS, SS-SS, and BI-FF beams
%   with their eigenspectra.
%% Flags
FLAG.DISPLAY_MODES = true;
FLAG.DISPLAY_SPECTRUM = true;
FLAG.CALCULATE_TIMESOLVER = true;

%% Parameters
n = 100;         % Mesh size, number of interior points
shapes = n;     % Number of shape functions

L = 1;           % Beam length
modeCount = 5;    % Number of displayed modes

%% BI-BI Beam
xvals = linspace(0,L,n+2);
deltax = mean(diff(xvals));
basis = zeros(shapes, n+2);
Dbasis = zeros(shapes, n+1);
D2basis = zeros(shapes, n);

% Generate shape functions and their respective derivatives
for i = 1:shapes
    basis(i,:) = xvals.*(xvals-L*ones(1,n+2)).*sin(pi*i*xvals);
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
        modeshapes(i,:) = modeshapes(i,:)/max(abs(modeshapes(i,:)));
        plot(xvals, modeshapes(i,:));
    end
    xlabel("Beam Position");
    ylabel("Relative Displacement");
    legend('1','2','3','4','5');
end


% Plot eigenvalue spectrum
espectrum = abs(diag(evals))/abs(evals(1,1));
scaled_espectrum = espectrum.^(1/4);

if FLAG.DISPLAY_SPECTRUM
    figure();
    hold on
    plot(scaled_espectrum(1:50), "LineStyle", "none", "Marker", ".");
end

%% BI-SS Beam
xvals = linspace(0,L,n+2);
deltax = mean(diff(xvals));
basis = zeros(shapes, n+2);
Dbasis = zeros(shapes, n+1);
D2basis = zeros(shapes, n);

% Generate shape functions and their respective derivatives
for i = 1:shapes
    basis(i,:) = xvals.*sin(i*pi*xvals);
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
        modeshapes(i,:) = modeshapes(i,:)/max(abs(modeshapes(i,:)));
        plot(xvals, modeshapes(i,:));
    end
    xlabel("Beam Position");
    ylabel("Relative Displacement");
    legend('1','2','3','4','5');
end


% Plot eigenvalue spectrum
espectrum = abs(diag(evals))/abs(evals(1,1));
scaled_espectrum = espectrum.^(1/4);

if FLAG.DISPLAY_SPECTRUM
    figure();
    hold on
    plot(scaled_espectrum(1:50), "LineStyle", "none", "Marker", ".");
end

%% BI-FF Beam
xvals = linspace(0,L,n+2);
deltax = mean(diff(xvals));
basis = zeros(shapes, n+2);
Dbasis = zeros(shapes, n+1);
D2basis = zeros(shapes, n);

% Generate shape functions and their respective derivatives
for i = 1:shapes
    basis(i,:) = xvals.*sin((2*i+1)*pi/2*xvals);
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
        modeshapes(i,:) = modeshapes(i,:)/max(abs(modeshapes(i,:)));
        plot(xvals, modeshapes(i,:));
    end
    xlabel("Beam Position");
    ylabel("Relative Displacement");
    legend('1','2','3','4','5');
end


% Plot eigenvalue spectrum
espectrum = abs(diag(evals))/abs(evals(1,1));
scaled_espectrum = espectrum.^(1/4);

if FLAG.DISPLAY_SPECTRUM
    figure();
    hold on
    plot(scaled_espectrum(1:50), "LineStyle", "none", "Marker", ".");
end

%% SS-SS Beam
xvals = linspace(0,L,n+2);
deltax = mean(diff(xvals));
basis = zeros(shapes, n+2);
Dbasis = zeros(shapes, n+1);
D2basis = zeros(shapes, n);

% Generate shape functions and their respective derivatives
for i = 1:shapes
    basis(i,:) = sin(i*xvals*pi);
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
        modeshapes(i,:) = modeshapes(i,:)/max(abs(modeshapes(i,:)));
        plot(xvals, modeshapes(i,:));
    end
    xlabel("Beam Position");
    ylabel("Relative Displacement");
    legend('1','2','3','4','5');
end


% Plot eigenvalue spectrum
espectrum = abs(diag(evals))/abs(evals(1,1));
scaled_espectrum = espectrum.^(1/4);

if FLAG.DISPLAY_SPECTRUM
    figure();
    hold on
    plot(scaled_espectrum(1:50), "LineStyle", "none", "Marker", ".");
end