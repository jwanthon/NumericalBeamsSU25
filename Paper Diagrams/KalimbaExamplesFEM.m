clc; clear all;
addpath('Functions\');

%% KalimbaExamplesFEM.m
%  Joseph Anthony
%
% Created:         11/11/25
% Last Modified:   11/23/25
%
% Description: Produces several diagrams related to the whammy-bar kalimba
%% Parameters
n = 100;         % Mesh size, number of interior points
shapes = n;     % Number of shape functions

L = 1;           % Beam length
modeCount = 5;    % Number of displayed modes

%% Single Mode
b = L/4;

xvals = linspace(0,L,n+2);
deltax = mean(diff(xvals));
basis = zeros(shapes, n+2);
Dbasis = zeros(shapes, n+1);
D2basis = zeros(shapes, n);

% Generate shape functions and their respective derivatives
for i = 1:shapes
    basis(i,:) = xvals.*(xvals-b*ones(1,n+2)).*sin((2*i-1)*pi/2*xvals);
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


figure();
hold on
for i = 1:modeCount
    modeshapes(i,:) = modeshapes(i,:)/max(abs(modeshapes(i,:)));
    plot(xvals, modeshapes(i,:));
end
xlabel("Beam Position");
ylabel("Relative Displacement");
legend('1','2','3','4','5');


% Plot eigenvalue spectrum
espectrum = abs(diag(evals))/abs(evals(1,1));
scaled_espectrum = espectrum.^(1/4);


figure();
hold on
plot(scaled_espectrum(1:50), "LineStyle", "none", "Marker", ".");

%% Find frequency spectrum

% Build ICs
ics = zeros(length(xvals),1);
for i = 1:length(ics)
    if xvals(i) > b
        ics(i) = -(xvals(i)-b)^2*(3*(L-b)-xvals(i));
    end
end

% Best fit of modes to ICs
ics_shapes = lsqr(basis', ics);
loglog(abs(diag(evals)), abs(ics_shapes), 'Marker', '.', 'LineStyle', 'none', 'MarkerSize', 15);
grid on