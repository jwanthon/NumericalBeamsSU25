% clc; clear all;
addpath('Functions\');
%% PiecewiseBeamNumericFEM.m
%  Joseph Anthony
%
% Created:         7/14/25
% Last Modified:   7/18/25
%
% Description: Numerically solves the Euler-Bernoulli beam equation 
%   for a built in-built in beam using numeric FEM and with different types
%   of piecewise functions.

n = 20;       % Mesh size, number of interior points
shapes = 20;     % Number of shape functions used
L = 1;         % Beam length
modeCount = 5;  % Number of displayed modes

%% Generate Different Bases
xvals = linspace(0,L,n+2);
deltax = mean(diff(xvals));
h = L/(shapes+1);

% Type 1, Hump: x²(x-2h)²(1/h⁴)
hump_basis = zeros(shapes,n+2);
hump_Dbasis = zeros(shapes, n+1);
hump_D2basis = zeros(shapes,n);

for i = 1:shapes
    x0 = (i-1)*L/(shapes+1);
    for j = 1:length(xvals)
        x = xvals(j);
        if (i-1)*h < x && x < (i+1)*h
            hump_basis(i,j) = (x-x0)^2*(x-x0-2*h)^2/h^4;
        end
    end
end

% Type 2, Wiggle:  x²(x-2h)²(x/h⁴-1/h³)
wigg_basis = zeros(shapes,n+2);
wigg_Dbasis = zeros(shapes, n+1);
wigg_D2basis = zeros(shapes,n);

for i = 1:shapes
    x0 = (i-1)*L/(shapes+1);
    for j = 1:length(xvals)
        x = xvals(j);
        if (i-1)*h < x && x < (i+1)*h
            wigg_basis(i,j) = (x-x0)^2*(x-x0-2*h)^2*(1/h^4*(x-x0)-1/h^3);
        end
    end
end
hold on
% Generate derivatives and normalize basis functions
for i = 1:shapes
    hump_basis(i,:) = hump_basis(i,:)/max(abs(hump_basis(i,:)));
    wigg_basis(i,:) = wigg_basis(i,:)/max(abs(wigg_basis(i,:)));
end

hump_basis(isnan(hump_basis)) = 0;
wigg_basis(isnan(wigg_basis)) = 0;

for i = 1:shapes
    hump_Dbasis(i,:) = diff(hump_basis(i,:))/deltax;
    hump_D2basis(i,:) = diff(hump_Dbasis(i,:))/deltax;
    wigg_Dbasis(i,:) = diff(wigg_basis(i,:))/deltax;
    wigg_D2basis(i,:) = diff(wigg_Dbasis(i,:))/deltax;
end  

beambasis = [hump_basis ; wigg_basis];
D2basis = [hump_D2basis ; wigg_D2basis];

%% Build Stiffness Matrix

% Hump submatrix
% hump_K = zeros(shapes);
% for row = 1:shapes
%     for col = 1:row
%         hump_K(row,col) = trapz(xvals(2:end-1), hump_D2basis(row,:).*hump_D2basis(col,:));
% 
%     end
% end
% 
% % Wiggle submatrix
% wigg_K = zeros(shapes);
% for row = 1:shapes
%     for col = 1:row
%         wigg_K(row,col) = trapz(xvals(2:end-1), hump_D2basis(row,:).*hump_D2basis(col,:));
% 
%     end
% end
% 
% % Build the full matrix
% K = [hump_K,        zeros(shapes);
%      zeros(shapes), wigg_K];
% 
% K = K + K' - diag(diag(K));

% Build the full matrix
K = zeros(2*shapes);
for row = 1:(2*shapes)
    for col = 1:(2*shapes)
        K(row,col) = trapz(xvals(2:end-1), D2basis(row,:).*D2basis(col,:));
    end
end
% K = K + K' - diag(diag(K));

%%
K(shapes+1:2*shapes, 1:shapes) = zeros(shapes);
K(1:shapes, shapes+1:2*shapes) = zeros(shapes);

%%
K = SYMBOLIC_K
% Determine sorted eigenbasis
[evecs, evals] = eig(K);
[~, index] = sort(diag(abs(evals))); 
evals = evals(index,index);
evecs = evecs(:, index);

% Display mode shapes
modeshapes = beambasis'*evecs;
modeshapes = modeshapes';

figure();
hold on
for i = [1, 4, 6, 7, 9]
    modeshapes(i,:) = modeshapes(i,:)/max(abs(modeshapes(i,:)));
    plot(xvals, modeshapes(i,:));
end
xlabel("Beam Position [m]");
ylabel("Relative Displacement");
title("Mode Shapes of a Built In-Built In Beam");


% Plot eigenvalue spectrum
figure();
espectrum = abs(diag(evals))/abs(evals(1,1));

scaled_espectrum = espectrum.^(1/4);

hold on
plot(scaled_espectrum, "LineStyle", "none", "Marker", ".");
title("Adjusted BI-BI Eigenvalue Magnitudes");
ylabel("Fourth Root of λ");

% [1, 4, 6, 7, 9]
% [2, 3, 5, 8, 10]

figure();
fig2 = visual_sparseMatrix(K);

%%
