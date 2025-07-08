clc; clear all;
addpath('Functions\');
%% FFFFBeamNumericFEM.m
%  Joseph Anthony
%
% Created:         7/7/25
% Last Modified:   6/12/25
%
% Description: Numerically solves the Euler-Bernoulli beam equation 
%   for a free-free free-free beam using FEM and determines the solution 
%   structure.
%
%% Compute Stiffness Matrix

rho         = 7.85;     % Mass per unit length [kg/m]
gamma       = 5.643;    % Flexural rigidity [N·m²]
L           = 1;        % Beam length [m]
n           = 1000;     % Mesh size, number of interior points
shapes      = 100;      % Number of shape functions
modeCount   = 5;        % Number of displayed modes

xvals = linspace(0,L,n+2);
deltax = mean(diff(xvals));
beambasis = zeros(shapes, n+2);
Dbeambasis = zeros(shapes, n+1);
D2beambasis = zeros(shapes, n);

% Note that both of our boundary conditions are natural boundary
%   conditions, so they affect the "leftover terms" when converting E-B to
%   the weak form.

% Generate shape functions and their respective derivatives 
%   Piecewise parabs (problem: 2nd derivs are nonzero at endpoints)
% for i = 1:shapes
%     beambasis(i,:) = ones(1,n+2) - 100*(xvals - ones(1,n+2)*i*L/(shapes+1)).^2;
% end
% beambasis(beambasis < 0) = 0;

%   Polys * sines (problem: shape functions are flat)
% beambasis = repmat(1/30*xvals.^6-1/10*L*xvals.^5+1/6*L^2*xvals.^4, shapes, 1);
% for i = 1:shapes
%     beambasis(i,:) = beambasis(i,:).*sin(i*(xvals+0.1*ones(1,n+2)));
% end
% 
% beambasis = [flip(beambasis(:,(n/2+2):end),2), beambasis(:,(n/2+2):end)];

for i = 1:shapes
    beambasis(i,:) = (xvals-L/2*ones(1,n+2)).^(i-1);
end

% % Display shape functions
% hold on
% for i = 1:shapes
%     plot(beambasis(i,:));
% end
% return

for i = 1:shapes
    Dbeambasis(i,:) = diff(beambasis(i,:))/deltax;
    D2beambasis(i,:) = diff(Dbeambasis(i,:))/deltax;
end

% Generate symmetric stiffness matrices
K = zeros(shapes);

for row = 1:shapes
    for col = 1:row
        K(row,col) = trapz(xvals(2:end-1), D2beambasis(row,:).*D2beambasis(col,:));
    end
end
K = K + K' - diag(diag(K));

[evecs, evals] = eig(K);
[~, index] = sort(diag(abs(evals))); 
evals = evals(index,index);
evecs = evecs(:, index);

% Display mode shapes
modeshapes = beambasis'*evecs;
modeshapes = modeshapes';

figure();
hold on
for i = 1:modeCount
    plot(xvals, modeshapes(i,:));
end
xlabel("Beam Position [m]");
ylabel("Relative Displacement");
title("Mode Shapes of a Free Free-Free Free Beam");


% Plot eigenvalue spectrum
figure();
espectrum = abs(diag(evals))/abs(evals(1,1));

scaled_espectrum = espectrum.^(1/4);

hold on
plot(scaled_espectrum, "LineStyle", "none", "Marker", ".");
title("Adjusted FF-FF Eigenvalue Magnitudes");
ylabel("Fourth Root of λ");