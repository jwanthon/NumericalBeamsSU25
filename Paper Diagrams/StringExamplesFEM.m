clc; clear all;
addpath('Functions\');
%% TODO:
% - Fix modeshapes using waves

%% StringExamplesFEM.m
%  Joseph Anthony
%
% Created:         10/21/25
% Last Modified:   10/21/25
%
% Description: Produces several diagrams of DD and DN string mode shapes
%   and basis functions.

n = 100;       % Mesh size, number of interior points
shapes = n;   % Number of shape functions, typically equal to n
L = pi;          % Beam length
modecount = 10;  % Number of displayed modes

xvals = linspace(0,L,n+2);
deltax = xvals(2)-xvals(1);
Dxvals = zeros(1, n+1); % Generate a mesh of xvals of derivative
for i = 1:length(xvals) - 1
    Dxvals(i) = (xvals(i) + xvals(i+1))/2;
end

% Initialize variables
basis = zeros(shapes, n+2);
Dbasis = zeros(shapes,n+1);
basis_DDmodes = basis;
basis_DNmodes = basis;
basis_tents = basis;
basis_waves = basis;
Dbasis_DDmodes = Dbasis;
Dbasis_DNmodes = Dbasis;
Dbasis_tents = Dbasis;
Dbasis_waves = Dbasis;

% Generate shape functions and their respective derivatives
for i = 1:shapes
    basis_DDmodes(i,:) = sin(i*xvals);
    Dbasis_DDmodes(i,:) = diff(basis_DDmodes(i,:))/deltax;
    
    basis_DNmodes(i,:) = sin((2*i-1)/2*xvals);
    Dbasis_DNmodes(i,:) = diff(basis_DNmodes(i,:))/deltax;

    basis_tents(i,i+1) = 1;
    Dbasis_tents(i,:) = diff(basis_tents(i,:))/deltax;

    basis_waves(i,i:end) = linspace(0,1,shapes-i+3); 
    Dbasis_waves(i,:) = diff(basis_waves(i,:))/deltax;
end

% Generate stiffness matrices
K_DDmodes = zeros(shapes);
K_DNmodes = zeros(shapes);
K_tents = zeros(shapes);
K_waves =  zeros(shapes);

for row = 1:shapes
    for col = 1:row
        K_DDmodes(row,col) = trapz(Dxvals, Dbasis_DDmodes(row,:).*Dbasis_DDmodes(col,:));
        K_DNmodes(row,col) = trapz(Dxvals, Dbasis_DNmodes(row,:).*Dbasis_DNmodes(col,:));
        K_tents(row,col) = trapz(Dxvals, Dbasis_tents(row,:).*Dbasis_tents(col,:));
        K_waves(row,col) = trapz(Dxvals, Dbasis_waves(row,:).*Dbasis_waves(col,:));
    end
end

K_DDmodes = K_DDmodes + K_DDmodes' - diag(diag(K_DDmodes));
K_DNmodes = K_DNmodes + K_DNmodes' - diag(diag(K_DNmodes));
K_tents = K_tents + K_tents' - diag(diag(K_tents));
K_waves = K_waves + K_waves' - diag(diag(K_waves));

% Determine sorted eigenbasis
[evecs_DDmodes, evals_DDmodes] = eig(K_DDmodes);
[~, index] = sort(diag(abs(evals_DDmodes))); 
evals_DDmodes = evals_DDmodes(index,index);
evecs_DDmodes = evecs_DDmodes(:, index);

[evecs_DNmodes, evals_DNmodes] = eig(K_DNmodes);
[~, index] = sort(diag(abs(evals_DNmodes))); 
evals_DNmodes = evals_DNmodes(index,index);
evecs_DNmodes = evecs_DNmodes(:, index);

[evecs_tents, evals_tents] = eig(K_tents);
[~, index] = sort(diag(abs(evals_tents))); 
evals_tents = evals_tents(index,index);
evecs_tents = evecs_tents(:, index);

[evecs_waves, evals_waves] = eig(K_waves);
[~, index] = sort(diag(abs(evecs_waves))); 
evecs_waves = evecs_waves(index,index);
evals_waves = evals_waves(:, index);


% Compute mode shapes
modes_DDmodes = basis_DDmodes'*evecs_DDmodes;
modes_DDmodes = modes_DDmodes';

modes_DNmodes = basis_DNmodes'*evecs_DNmodes;
modes_DNmodes = modes_DNmodes';

modes_tents = basis_tents'*evecs_tents;
modes_tents = modes_tents';

modes_waves = basis_waves'*evecs_waves;
modes_waves = modes_waves';

figure();
hold on
for i = 1:modecount
    plot(xvals, modes_tents(i,:));
end
xlabel("String Position [m]");
ylabel("Relative Displacement");
title("Mode Shapes of a Dirichlet-Dirichley String");

% 
% % Plot eigenvalue spectrum
% espectrum = abs(diag(evals))/abs(evals(1,1));
% rootvals = sqrt(espectrum);
% figure();
% plot(rootvals, "LineStyle", "none", "Marker", ".");
% title("Squareroots of Eigenvalue Magnitudes");
